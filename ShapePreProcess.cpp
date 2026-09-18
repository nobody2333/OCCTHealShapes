#include "ShapePreProcess.h"

#include "GeometryUtils.h"
#include "TopologyUtils.h"
#include "RepairDiagnostics.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRep_Tool.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Wire.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepCheck_ListIteratorOfListOfStatus.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <BRepTools.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <GC_MakeSegment.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomProjLib.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>

namespace {
constexpr int kMaxAutomaticMergeIterations = 1000;
constexpr double kBilateralAnchorLengthDifferenceTolerance = 1.0e-5;

struct BRepCheckFailureInfo {
    TopoDS_Shape shape;
    TopAbs_ShapeEnum type = TopAbs_SHAPE;
    std::vector<BRepCheck_Status> statuses;
};

std::string failureText(const Standard_Failure& failure)
{
    const Standard_CString message = failure.GetMessageString();
    return message != nullptr ? std::string(message) : std::string("OCCT failure");
}

bool containsStatus(const std::vector<BRepCheck_Status>& statuses, BRepCheck_Status status)
{
    return std::find(statuses.begin(), statuses.end(), status) != statuses.end();
}

bool statusesAreSubset(
    const std::vector<BRepCheck_Status>& statuses,
    const std::vector<BRepCheck_Status>& baselineStatuses)
{
    for (const BRepCheck_Status status : statuses) {
        if (!containsStatus(baselineStatuses, status)) {
            return false;
        }
    }
    return true;
}

std::vector<BRepCheckFailureInfo> collectBRepCheckFailures(
    const TopoDS_Shape& shape,
    std::string* errorMessage = nullptr)
{
    std::vector<BRepCheckFailureInfo> failures;
    if (shape.IsNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = "shape is null";
        }
        return failures;
    }

    try {
        BRepCheck_Analyzer analyzer(shape, Standard_True);
        TopTools_IndexedMapOfShape allShapes;
        TopExp::MapShapes(shape, allShapes);

        for (int i = 1; i <= allShapes.Extent(); ++i) {
            const TopoDS_Shape& subShape = allShapes.FindKey(i);
            Handle(BRepCheck_Result) result = analyzer.Result(subShape);
            if (result.IsNull()) {
                continue;
            }

            BRepCheckFailureInfo failure;
            failure.shape = subShape;
            failure.type = subShape.ShapeType();

            const BRepCheck_ListOfStatus& statuses = result->Status();
            for (BRepCheck_ListIteratorOfListOfStatus it(statuses); it.More(); it.Next()) {
                const BRepCheck_Status status = it.Value();
                if (status == BRepCheck_NoError || containsStatus(failure.statuses, status)) {
                    continue;
                }
                failure.statuses.push_back(status);
            }

            if (!failure.statuses.empty()) {
                failures.push_back(failure);
            }
        }
    }
    catch (Standard_Failure& failure) {
        if (errorMessage != nullptr) {
            *errorMessage = "BRepCheck failure collection failed: " + failureText(failure);
        }
    }

    return failures;
}

bool isExistingBRepCheckFailure(
    const BRepCheckFailureInfo& failure,
    const std::vector<BRepCheckFailureInfo>& baselineFailures)
{
    for (const BRepCheckFailureInfo& baseline : baselineFailures) {
        if (failure.type != baseline.type || !failure.shape.IsSame(baseline.shape)) {
            continue;
        }

        if (statusesAreSubset(failure.statuses, baseline.statuses)) {
            return true;
        }
    }

    return false;
}

std::string describeBRepCheckFailures(const std::vector<BRepCheckFailureInfo>& failures)
{
    std::ostringstream os;
    int reported = 0;
    constexpr int maxReported = 30;
    for (const BRepCheckFailureInfo& failure : failures) {
        if (reported >= maxReported) {
            os << " ...";
            break;
        }

        os << " [type=" << static_cast<int>(failure.type) << " status=";
        for (std::size_t i = 0; i < failure.statuses.size(); ++i) {
            if (i > 0) {
                os << ",";
            }
            os << static_cast<int>(failure.statuses[i]);
        }
        os << "]";
        ++reported;
    }

    return os.str();
}

bool listContainsSameShape(const TopTools_ListOfShape& shapes, const TopoDS_Shape& target)
{
    for (TopTools_ListIteratorOfListOfShape it(shapes); it.More(); it.Next()) {
        if (it.Value().IsSame(target)) {
            return true;
        }
    }
    return false;
}

bool listsShareShape(const TopTools_ListOfShape& lhs, const TopTools_ListOfShape& rhs)
{
    for (TopTools_ListIteratorOfListOfShape it(lhs); it.More(); it.Next()) {
        if (listContainsSameShape(rhs, it.Value())) {
            return true;
        }
    }
    return false;
}


bool edgesShareFaceWire(
    const TopoDS_Face& face, const std::vector<TopoDS_Edge>& required)
{
    for (TopExp_Explorer wireIt(face, TopAbs_WIRE); wireIt.More(); wireIt.Next()) {
        std::size_t found = 0;
        for (const TopoDS_Edge& requiredEdge : required) {
            for (TopExp_Explorer edgeIt(wireIt.Current(), TopAbs_EDGE);
                 edgeIt.More(); edgeIt.Next()) {
                if (edgeIt.Current().IsSame(requiredEdge)) {
                    ++found;
                    break;
                }
            }
        }
        if (found == required.size()) {
            return true;
        }
    }
    return false;
}

double pointToEdgeDistance(const gp_Pnt& point, const TopoDS_Edge& edge)
{
    Standard_Real first = 0.0;
    Standard_Real last = 0.0;
    const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);
    if (curve.IsNull() || !std::isfinite(first) || !std::isfinite(last)
        || last - first <= Precision::PConfusion()) {
        return std::numeric_limits<double>::infinity();
    }
    // The closest point can be an endpoint. OCCT's interior projector may
    // report no solution for that case, which must not become infinite gap.
    double nearest = std::min(
        point.Distance(curve->Value(first)),
        point.Distance(curve->Value(last)));
    GeomAPI_ProjectPointOnCurve projector(point, curve, first, last);
    if (projector.NbPoints() > 0) {
        nearest = std::min(nearest, projector.LowerDistance());
    }
    return nearest;
}


TopoDS_Edge orientEdgeFromTo(
    const TopoDS_Edge& source,
    const TopoDS_Vertex& from,
    const TopoDS_Vertex& to)
{
    TopoDS_Edge edge = source;
    TopoDS_Vertex first;
    TopoDS_Vertex last;
    TopExp::Vertices(edge, first, last, Standard_True);
    if (first.IsSame(from) && last.IsSame(to)) {
        return edge;
    }
    if (first.IsSame(to) && last.IsSame(from)) {
        edge.Reverse();
        return edge;
    }
    return TopoDS_Edge();
}

double pointToPathDistance(
    const gp_Pnt& point,
    const std::vector<TopoDS_Edge>& path)
{
    double nearest = std::numeric_limits<double>::infinity();
    for (const TopoDS_Edge& edge : path) {
        nearest = std::min(nearest, pointToEdgeDistance(point, edge));
    }
    return nearest;
}

double maximumPathGap(
    const std::vector<TopoDS_Edge>& firstPath,
    const std::vector<TopoDS_Edge>& secondPath,
    double maxAllowed)
{
    if (firstPath.empty() || secondPath.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    double maximum = 0.0;
    for (const auto* direction : {&firstPath, &secondPath}) {
        const auto& target = direction == &firstPath ? secondPath : firstPath;
        for (const TopoDS_Edge& edge : *direction) {
            Standard_Real first = 0.0;
            Standard_Real last = 0.0;
            const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);
            if (curve.IsNull() || !std::isfinite(first) || !std::isfinite(last)
                || last - first <= Precision::PConfusion()) {
                return std::numeric_limits<double>::infinity();
            }
            const double length = GeometryUtils::ComputeEdgeLength(edge);
            if (!std::isfinite(length)) {
                return std::numeric_limits<double>::infinity();
            }
            const int samples = static_cast<int>(std::min(64.0,
                std::max(8.0, std::ceil(length /
                    std::max(4.0 * maxAllowed, Precision::Confusion())))));
            for (int sample = 0; sample <= samples; ++sample) {
                const gp_Pnt point = curve->Value(
                    first + (last - first) * sample / samples);
                maximum = std::max(maximum, pointToPathDistance(point, target));
                if (maximum > maxAllowed) {
                    return maximum;
                }
            }
        }
    }
    return maximum;
}

bool pathIsSimple(const std::vector<TopoDS_Vertex>& vertices)
{
    if (vertices.size() < 2) {
        return false;
    }
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (vertices[i].IsNull()) {
            return false;
        }
        for (std::size_t j = i + 1; j < vertices.size(); ++j) {
            if (vertices[i].IsSame(vertices[j])) {
                return false;
            }
        }
    }
    return true;
}

bool uniqueEdgeUseOrientation(
    const TopoDS_Face& face,
    const TopoDS_Edge& edge,
    TopAbs_Orientation& orientation)
{
    int useCount = 0;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next()) {
        if (!it.Current().IsSame(edge)) {
            continue;
        }
        orientation = it.Current().Orientation();
        ++useCount;
    }
    return useCount == 1 &&
        (orientation == TopAbs_FORWARD || orientation == TopAbs_REVERSED);
}
namespace bilateral {
struct Params {
    double maxGapWidth = 0.5;
    double maxEndpointCollapseDistance = 0.8;
    double endpointLengthDifferenceTolerance = 1.0e-5;
    double pcurveTolerance = 1.0e-1;
    double maxRebuiltTolerance = 1.0e-1;
    double maxRelativeAreaChange = 0.15;
    int maxBoundaryEdges = 8;
};
using Edges = std::vector<TopoDS_Edge>;
void requireValidRepairCondition(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
double normalizePositiveTolerance(double value, double fallback) {
    return std::isfinite(value) && value >= Precision::Confusion() ? value : fallback;
}
TopoDS_Vertex getEdgeStartVertex(const TopoDS_Edge& e) { return TopExp::FirstVertex(e, Standard_True); }
TopoDS_Vertex getEdgeEndVertex(const TopoDS_Edge& e) { return TopExp::LastVertex(e, Standard_True); }
TopoDS_Edge orientEdgeFromVertex(TopoDS_Edge e, const TopoDS_Vertex& v) {
    if (!getEdgeStartVertex(e).IsSame(v)) e.Reverse();
    requireValidRepairCondition(getEdgeStartVertex(e).IsSame(v), "disconnected edge use");
    return e;
}
Handle(Geom_Curve) getEdgeCurveInWorldCoordinates(const TopoDS_Edge& e, double& first, double& last) {
    auto c = BRep_Tool::Curve(e, first, last); // 带 Location 的世界坐标曲线。
    requireValidRepairCondition(!c.IsNull() && std::isfinite(first) && std::isfinite(last)
        && last > first, "missing or unbounded 3D curve");
    return c;
}
double computeBoundaryLength(const Edges& edges) {
    double sum = 0;
    for (const auto& e : edges) sum += GeometryUtils::ComputeEdgeLength(e);
    return sum;
}
double computePointToBoundaryDistance(const gp_Pnt& point, const Edges& edges) {
    double best = std::numeric_limits<double>::infinity();
    for (const auto& e : edges) {
        double f, l;
        const auto c = getEdgeCurveInWorldCoordinates(e, f, l);
        best = std::min(best, std::min(point.Distance(c->Value(f)), point.Distance(c->Value(l))));
        GeomAPI_ProjectPointOnCurve projection(point, c, f, l);
        if (projection.NbPoints()) best = std::min(best, projection.LowerDistance());
    }
    return best;
}
template<class F> void forEachBoundarySample(const Edges& edges, F action) {
    for (const auto& e : edges) {
        double f, l;
        const auto c = getEdgeCurveInWorldCoordinates(e, f, l);
        // 每段独立覆盖，不能让较短的边界片段漏采。只是筛选，最终另做精度检查。
        for (int i = 0; i <= 24; ++i) action(c->Value(f + (l - f) * i / 24.0));
    }
}
bool containsSameEdge(const Edges& edges, const TopoDS_Shape& shape) {
    return std::any_of(edges.begin(), edges.end(), [&](const auto& e) { return e.IsSame(shape); });
}
struct Candidate { TopoDS_Vertex p, q; TopoDS_Face a, b; Edges sideA, sideB; };
struct ProjectionQuality {
    double meanSideImbalance = std::numeric_limits<double>::infinity();
    double maxSurfaceDeviation = std::numeric_limits<double>::infinity();
    double maxEndpointDeviation = std::numeric_limits<double>::infinity();
};
struct Index {
    TopTools_IndexedDataMapOfShapeListOfShape ef, ve;
    TopTools_IndexedMapOfShape anchors;
    std::vector<Edges> anchorPairs;
    explicit Index(const TopoDS_Shape& shape, const Params& params, std::vector<std::string>& messages) {
        TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_EDGE, TopAbs_FACE, ef);
        TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, ve);
        for (int i = 1; i <= ve.Extent(); ++i) {
            Edges boundary;
            std::vector<double> lengths;
            for (const auto& shapeEdge : ve(i)) {
                if (!ef.Contains(shapeEdge) || ef.FindFromKey(shapeEdge).Extent() != 1) continue;
                const auto e = TopoDS::Edge(shapeEdge);
                if (BRep_Tool::Degenerated(e) || getEdgeStartVertex(e).IsNull() || getEdgeEndVertex(e).IsNull()
                    || getEdgeStartVertex(e).IsSame(getEdgeEndVertex(e))) continue;
                const double len = GeometryUtils::ComputeEdgeLength(e);
                if (!std::isfinite(len) || len <= Precision::Confusion()
                    || len > params.maxEndpointCollapseDistance) continue;
                boundary.push_back(e); lengths.push_back(len);
            }
            Edges pair;
            int matches = 0;
            for (std::size_t a = 0; a < boundary.size(); ++a) for (std::size_t b = a + 1; b < boundary.size(); ++b) {
                if (std::abs(lengths[a] - lengths[b]) > params.endpointLengthDifferenceTolerance) continue;
                // 本模块修复两面之间的缺陷，排除同一面上的普通等长边角。
                if (ef.FindFromKey(boundary[a]).First().IsSame(ef.FindFromKey(boundary[b]).First())) continue;
                ++matches; pair = {boundary[a], boundary[b]};
                const auto pt = BRep_Tool::Pnt(TopoDS::Vertex(ve.FindKey(i)));
                std::ostringstream os;
                os.precision(12);
                os << "Symmetric boundary pair: vertex=" << i << ", position=(" << pt.X() << ',' << pt.Y() << ',' << pt.Z()
                   << "), edges=" << ef.FindIndex(boundary[a]) << '/' << ef.FindIndex(boundary[b])
                   << ", lengths=" << lengths[a] << '/' << lengths[b]
                   << ", difference=" << std::abs(lengths[a] - lengths[b]);
                messages.push_back(os.str());
            }
            if (matches == 1) { anchors.Add(ve.FindKey(i)); anchorPairs.push_back(pair); }
            else if (matches > 1) messages.push_back("Ambiguous equal-length pairs at vertex " + std::to_string(i));
        }
    }
    bool traceBoundaryToAnchor(const TopoDS_Vertex& p, TopoDS_Edge e, int limit,
              Edges& chain, TopoDS_Face& face, TopoDS_Vertex& q) const {
        TopoDS_Vertex v = p;
        while (static_cast<int>(chain.size()) < limit) {
            if (!ef.Contains(e) || ef.FindFromKey(e).Extent() != 1 || containsSameEdge(chain, e)) return false;
            const auto owner = TopoDS::Face(ef.FindFromKey(e).First());
            if (face.IsNull()) face = owner;
            if (!owner.IsSame(face) || BRep_Tool::Degenerated(e) || BRep_Tool::IsClosed(e, face)) return false;
            e = orientEdgeFromVertex(e, v);
            if (getEdgeEndVertex(e).IsNull() || getEdgeStartVertex(e).IsSame(getEdgeEndVertex(e))) return false;
            chain.push_back(e);
            v = getEdgeEndVertex(e);
            if (anchors.Contains(v)) {
                if (!containsSameEdge(anchorPairs[anchors.FindIndex(v) - 1], e)) return false;
                q = v; return !q.IsSame(p);
            }
            // 内部待删除顶点不能被第三条边使用，避免损伤其它 face。
            if (!ve.Contains(v) || ve.FindFromKey(v).Extent() != 2) return false;
            bool found = false;
            for (const auto& next : ve.FindFromKey(v)) {
                if (!next.IsSame(e)) { e = TopoDS::Edge(next); found = true; break; }
            }
            if (!found) return false;
        }
        return false;
    }
};
void validateEndpointCollapseCandidate(const Candidate& c, const Params& p) {
    // 每侧只有一条主体长边，允许在其两端各有一条错误短边。
    auto tips = [&](const Edges& side) {
        requireValidRepairCondition(!side.empty() && side.size() <= 3, "unsupported boundary: expected short/main/short chain");
        auto longest = std::max_element(side.begin(), side.end(), [](const auto& x, const auto& y) {
            return GeometryUtils::ComputeEdgeLength(x) < GeometryUtils::ComputeEdgeLength(y);
        });
        const auto k = static_cast<std::size_t>(longest - side.begin());
        requireValidRepairCondition(k <= 1 && side.size() - k <= 2, "more than one end-cap edge");
        for (std::size_t i = 0; i < side.size(); ++i)
            if (i != k) {
                const double capLength = GeometryUtils::ComputeEdgeLength(side[i]);
                requireValidRepairCondition(capLength <= p.maxEndpointCollapseDistance,
                    "end-cap edge exceeds collapse budget: length=" + std::to_string(capLength)
                    + ", budget=" + std::to_string(p.maxEndpointCollapseDistance));
            }
        return std::make_pair(getEdgeStartVertex(*longest), getEdgeEndVertex(*longest));
    };
    const auto ta = tips(c.sideA), tb = tips(c.sideB);
    for (int end = 0; end < 2; ++end) {
        const gp_Pnt a = BRep_Tool::Pnt(end ? ta.second : ta.first);
        const gp_Pnt b = BRep_Tool::Pnt(end ? tb.second : tb.first);
        const gp_Pnt target = BRep_Tool::Pnt(end ? c.q : c.p);
        const gp_Pnt midpoint((a.XYZ() + b.XYZ()) * 0.5);
        requireValidRepairCondition(target.Distance(a) <= p.maxEndpointCollapseDistance
            && target.Distance(b) <= p.maxEndpointCollapseDistance
            && target.Distance(midpoint) <= p.maxEndpointCollapseDistance * 0.5,
            "existing anchor is not a local midpoint candidate");
    }
    requireValidRepairCondition(c.sideA.size() > 1 || c.sideB.size() > 1, "no incorrect end-cap edge to remove");
}
std::vector<Candidate> findBilateralCuspCandidates(const TopoDS_Shape& shape, const Params& params,
                            std::vector<std::string>& messages) {
    Index index(shape, params, messages);
    std::vector<Candidate> result;
    std::set<std::pair<int, int>> seen;
    messages.push_back("Anchor search: equal-length boundary pairs=" + std::to_string(index.anchors.Extent())
        + ", lengthDifferenceTolerance=" + std::to_string(params.endpointLengthDifferenceTolerance));
    for (int i = 1; i <= index.anchors.Extent(); ++i) {
        const auto p = TopoDS::Vertex(index.anchors(i));
        const auto& first = index.anchorPairs[i - 1];
        for (int direction = 0; direction < 2; ++direction) {
            Candidate c; c.p = p;
            if (!index.traceBoundaryToAnchor(p, first[direction], params.maxBoundaryEdges, c.sideA, c.a, c.q)) {
                messages.push_back("Anchor " + std::to_string(index.ve.FindIndex(p))
                    + ": free-boundary search failed, startingEdge=" + std::to_string(index.ef.FindIndex(first[direction])));
                continue;
            }
            const int pId = index.ve.FindIndex(p), qId = index.ve.FindIndex(c.q);
            const std::pair<int, int> key{std::min(pId, qId), std::max(pId, qId)};
            if (seen.count(key)) continue;
            const auto point = BRep_Tool::Pnt(c.q);
            std::ostringstream os;
            os << "Candidate anchors " << pId << " -> " << qId
               << ", Q=(" << point.X() << ',' << point.Y() << ',' << point.Z() << ")";
            TopoDS_Vertex otherEnd;
            if (!index.traceBoundaryToAnchor(p, first[1 - direction], params.maxBoundaryEdges, c.sideB, c.b, otherEnd)
                || !otherEnd.IsSame(c.q)) {
                messages.push_back(os.str() + ": rejected: opposite free boundary does not reach the same symmetric endpoint"); continue;
            }
            seen.insert(key);
            try {
                requireValidRepairCondition(!c.a.IsSame(c.b), "both paths belong to the same face");
                validateEndpointCollapseCandidate(c, params);
                // 扩展候选后不能把一条属于第三个面的共享边误当作可删除短边。
                for (const auto& side : {c.sideA, c.sideB}) for (const auto& e : side) {
                    for (const auto& owner : index.ef.FindFromKey(e))
                        requireValidRepairCondition(owner.IsSame(c.a) || owner.IsSame(c.b),
                            "boundary interval uses third-face edge " + std::to_string(index.ef.FindIndex(e))
                            + "; adjacent-face retrimming required");
                }
                messages.push_back(os.str() + ": selected, chainEdges=" + std::to_string(c.sideA.size())
                    + "/" + std::to_string(c.sideB.size()));
                result.push_back(c);
            }
            catch (const std::exception& e) { messages.push_back(os.str() + ": rejected: " + e.what()); }
        }
    }
    return result;
}
// 将 P-Q 连线投影到当前试算 A 的原支撑面；调用方会交换 A/B 做独立试算。
Edges buildSharedEdgeByProjectingPQ(const Candidate& c, const Params& params,
                          std::vector<std::string>& messages, ProjectionQuality& quality) {
    const auto surface = BRep_Tool::Surface(c.a); // 世界坐标，包含原 Location。
    requireValidRepairCondition(!surface.IsNull(), "projection surface A is null");
    requireValidRepairCondition(!surface->IsUPeriodic() && !surface->IsVPeriodic(), "periodic projection surface unsupported");
    const auto p = BRep_Tool::Pnt(c.p), q = BRep_Tool::Pnt(c.q);
    requireValidRepairCondition(p.Distance(q) > Precision::Confusion(), "P-Q segment has zero length");
    GC_MakeSegment segment(p, q);
    requireValidRepairCondition(segment.IsDone(), "cannot construct P-Q line segment");
    messages.push_back("Shared boundary method: project P-Q segment onto trial surface A");
    // Project 对整条有限线段作法向投影并生成空间曲线，不使用原 P-Q 直线作共享曲线。
    const Handle(Geom_Curve) projected = GeomProjLib::Project(segment.Value(), surface);
    requireValidRepairCondition(!projected.IsNull(), "P-Q normal projection returned no curve");
    const double first = projected->FirstParameter(), last = projected->LastParameter();
    requireValidRepairCondition(std::isfinite(first) && std::isfinite(last)
        && !Precision::IsInfinite(first) && !Precision::IsInfinite(last)
        && last - first > Precision::PConfusion(), "projected curve has invalid finite range");
    const auto head = projected->Value(first), tail = projected->Value(last);
    const bool reversed = p.Distance(tail) + q.Distance(head) < p.Distance(head) + q.Distance(tail);
    const auto firstVertex = reversed ? c.q : c.p;
    const auto lastVertex = reversed ? c.p : c.q;
    const double pDeviation = p.Distance(reversed ? tail : head);
    const double qDeviation = q.Distance(reversed ? head : tail);
    quality.maxEndpointDeviation = std::max(pDeviation, qDeviation);
    std::ostringstream diagnostics; diagnostics.precision(12);
    diagnostics << "P-Q projection endpoint deviations: P=" << pDeviation << ", Q=" << qDeviation
                << ", maxRebuiltTolerance=" << params.maxRebuiltTolerance;
    messages.push_back(diagnostics.str());
    // 顶点仅在隔离 trial 中按实测值增加容差；不移动 P/Q，也不截掉投影的缺失部分。
    BRep_Builder builder;
    for (const auto& vertex : {c.p, c.q}) {
        const double deviation = vertex.IsSame(c.p) ? pDeviation : qDeviation;
        const double existing = BRep_Tool::Tolerance(vertex);
        requireValidRepairCondition(std::isfinite(existing) && existing >= 0.0
            && std::isfinite(deviation), "invalid projected endpoint tolerance or deviation");
        const double needed = std::max(existing, deviation + Precision::Confusion());
        if (needed > existing && needed > params.maxRebuiltTolerance) {
            std::ostringstream error; error.precision(12);
            error << "projected endpoint tolerance increase exceeds budget: vertex="
                  << (vertex.IsSame(c.p) ? "P" : "Q")
                  << ", existing=" << existing << ", deviation=" << deviation
                  << ", required=" << needed << ", maximum=" << params.maxRebuiltTolerance;
            requireValidRepairCondition(false, error.str());
        }
        if (needed > existing) builder.UpdateVertex(vertex, needed);
    }
    BRepBuilderAPI_MakeEdge maker(projected, firstVertex, lastVertex, first, last);
    requireValidRepairCondition(maker.IsDone(), "cannot construct shared projected edge with P/Q vertices");
    const Edges shared{orientEdgeFromVertex(maker.Edge(), c.p)};
    const double projectedLength = computeBoundaryLength(shared);
    requireValidRepairCondition(std::isfinite(projectedLength) && projectedLength > Precision::Confusion()
        && projectedLength <= 1.5 * std::max(computeBoundaryLength(c.sideA), computeBoundaryLength(c.sideB)), "projected curve detours or degenerates");
    double imbalanceSum = 0.0;
    int imbalanceSamples = 0;
    forEachBoundarySample(shared, [&](const gp_Pnt& point) {
        const double distanceA = computePointToBoundaryDistance(point, c.sideA);
        const double distanceB = computePointToBoundaryDistance(point, c.sideB);
        requireValidRepairCondition(std::isfinite(distanceA) && std::isfinite(distanceB)
            && distanceA <= params.maxGapWidth && distanceB <= params.maxGapWidth,
            "projected edge leaves defect band");
        imbalanceSum += std::abs(distanceA - distanceB);
        ++imbalanceSamples;
    });
    requireValidRepairCondition(imbalanceSamples > 0, "projected edge has no valid samples");
    quality.meanSideImbalance = imbalanceSum / imbalanceSamples;
    forEachBoundarySample(c.sideA, [&](const gp_Pnt& point) {
        requireValidRepairCondition(computePointToBoundaryDistance(point, shared) <= params.maxGapWidth, "projected edge misses side A");
    });
    forEachBoundarySample(c.sideB, [&](const gp_Pnt& point) {
        requireValidRepairCondition(computePointToBoundaryDistance(point, shared) <= params.maxGapWidth, "projected edge misses side B");
    });
    messages.push_back("Shared projected edge constructed: length=" + std::to_string(projectedLength));
    return shared;
}
TopoDS_Face rebuildFaceWithSharedBoundary(const TopoDS_Face& input, const Edges& removed, const Edges& shared,
                   const Candidate& c, const Params& params) {
    const auto face = TopoDS::Face(input.Oriented(TopAbs_FORWARD));
    TopTools_IndexedMapOfShape wires;
    TopExp::MapShapes(face, TopAbs_WIRE, wires);
    requireValidRepairCondition(wires.Extent() == 1, "inner/multiple wires unsupported");
    const auto oldWire = BRepTools::OuterWire(face);
    Edges ordered;
    for (BRepTools_WireExplorer it(oldWire, face); it.More(); it.Next()) ordered.push_back(it.Current());
    TopTools_IndexedMapOfShape oldEdges;
    TopExp::MapShapes(oldWire, TopAbs_EDGE, oldEdges);
    requireValidRepairCondition(ordered.size() == static_cast<std::size_t>(oldEdges.Extent()), "unordered/repeated wire edges");
    std::size_t begin = ordered.size();
    for (std::size_t i = 0; i < ordered.size(); ++i)
        if (containsSameEdge(removed, ordered[i]) && !containsSameEdge(removed, ordered[(i + ordered.size() - 1) % ordered.size()])) {
            requireValidRepairCondition(begin == ordered.size(), "removed boundary is not a single wire interval"); begin = i;
        }
    requireValidRepairCondition(begin < ordered.size(), "cannot locate removed interval");
    std::rotate(ordered.begin(), ordered.begin() + begin, ordered.end());
    std::size_t count = 0;
    while (count < ordered.size() && containsSameEdge(removed, ordered[count])) ++count;
    requireValidRepairCondition(count == removed.size() && count < ordered.size(), "boundary interval mismatch");
    Edges replacement = shared;
    const auto first = getEdgeStartVertex(ordered[0]);
    if (!first.IsSame(getEdgeStartVertex(replacement.front()))) {
        std::reverse(replacement.begin(), replacement.end());
        for (auto& e : replacement) e.Reverse();
    }
    requireValidRepairCondition(first.IsSame(getEdgeStartVertex(replacement.front()))
        && getEdgeEndVertex(ordered[count - 1]).IsSame(getEdgeEndVertex(replacement.back())), "replacement anchors mismatch");
    replacement.insert(replacement.end(), ordered.begin() + count, ordered.end());
    BRep_Builder builder;
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    for (std::size_t i = 0; i < replacement.size(); ++i) {
        requireValidRepairCondition(getEdgeEndVertex(replacement[i]).IsSame(getEdgeStartVertex(replacement[(i + 1) % replacement.size()])), "wire is not closed");
        requireValidRepairCondition(!BRep_Tool::Degenerated(replacement[i]) && !BRep_Tool::IsClosed(replacement[i], face),
            "seam/degenerate wire edge unsupported");
        requireValidRepairCondition(GeometryUtils::ComputeEdgeLength(replacement[i]) > Precision::Confusion(),
            "near-zero wire edge after replacement");
        builder.Add(wire, replacement[i]);
    }
    wire.Closed(Standard_True);
    TopLoc_Location loc;
    const auto surface = BRep_Tool::Surface(face, loc);
    TopoDS_Face rebuilt;
    builder.MakeFace(rebuilt, surface, loc, BRep_Tool::Tolerance(face));
    builder.Add(rebuilt, wire);
    for (auto edge : shared) {
        std::string error;
        const bool repaired = GeometryUtils::RepairPCurveForEdgeOnFace(edge, rebuilt, params.pcurveTolerance, &error);
        requireValidRepairCondition(repaired, "cannot create shared pcurve: " + error);
        BRepLib_CheckCurveOnSurface check(edge, rebuilt);
        check.Perform();
        requireValidRepairCondition(check.IsDone() && std::isfinite(check.MaxDistance())
            && check.MaxDistance() <= params.pcurveTolerance, "shared curve deviates from original surface");
        const double needed = std::max(BRep_Tool::Tolerance(edge), check.MaxDistance() + Precision::Confusion());
        if (needed > params.maxRebuiltTolerance) {
            std::ostringstream error;
            error.precision(12);
            error << "shared edge tolerance exceeds budget: edge="
                  << BRep_Tool::Tolerance(edge) << ", deviation="
                  << check.MaxDistance() << ", limit=" << params.maxRebuiltTolerance;
            requireValidRepairCondition(false, error.str());
        }
        builder.UpdateEdge(edge, needed);
    }
    BRepCheck_Wire wireCheck(wire);
    TopoDS_Edge e1, e2;
    requireValidRepairCondition(wireCheck.SelfIntersect(rebuilt, e1, e2) == BRepCheck_NoError, "rebuilt UV wire self-intersects");
    BRepCheck_Analyzer check(rebuilt, Standard_True);
    requireValidRepairCondition(check.IsValid(), "rebuilt face invalid: " + RepairDiagnostics::DescribeBRepCheckFailures(check, rebuilt));
    GProp_GProps oldArea, newArea;
    BRepGProp::SurfaceProperties(face, oldArea);
    BRepGProp::SurfaceProperties(rebuilt, newArea);
    requireValidRepairCondition(oldArea.Mass() > Precision::Confusion() && std::isfinite(newArea.Mass())
        && newArea.Mass() > Precision::Confusion()
        && std::abs(newArea.Mass() - oldArea.Mass()) <= oldArea.Mass() * params.maxRelativeAreaChange,
        "excessive trimmed area change");
    // 远离缺陷的原面内部点必须仍在重裁后的区域内；不简单选择最大环。
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    int preserved = 0;
    const auto worldSurface = BRep_Tool::Surface(face);
    for (int u = 1; u < 8; ++u) for (int v = 1; v < 8; ++v) {
        const gp_Pnt2d uv(u0 + (u1 - u0) * u / 8.0, v0 + (v1 - v0) * v / 8.0);
        BRepClass_FaceClassifier before(face, uv, Precision::PConfusion());
        const auto pt = worldSurface->Value(uv.X(), uv.Y());
        if (before.State() != TopAbs_IN || computePointToBoundaryDistance(pt, c.sideA) <= params.maxGapWidth
            || computePointToBoundaryDistance(pt, c.sideB) <= params.maxGapWidth) continue;
        BRepClass_FaceClassifier after(rebuilt, uv, Precision::PConfusion());
        requireValidRepairCondition(after.State() == TopAbs_IN || after.State() == TopAbs_ON, "lost trusted face interior");
        ++preserved;
    }
    requireValidRepairCondition(preserved > 0, "no trusted interior sample outside defect band");
    rebuilt.Orientation(input.Orientation());
    return rebuilt;
}

// 方向调整作用于整个 face use，不能单独反转闭合 wire 中的一条 edge。
// 仅枚举两个新面的四种朝向；外围相邻面保持不变。
TopoDS_Shape replaceFacesWithConsistentOrientation(
    const TopoDS_Shape& trial, const Candidate& c,
    const TopoDS_Face& a, const TopoDS_Face& b,
    const Edges& shared, std::vector<std::string>& messages) {
    const auto edgeUse = [](const TopoDS_Face& face, const TopoDS_Shape& edge,
                            TopAbs_Orientation& orientation) {
        int count = 0;
        for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next()) {
            if (!it.Current().IsSame(edge)) continue;
            orientation = it.Current().Orientation();
            ++count;
        }
        return count == 1 && (orientation == TopAbs_FORWARD || orientation == TopAbs_REVERSED);
    };
    TopTools_IndexedMapOfShape boundary;
    TopExp::MapShapes(a, TopAbs_EDGE, boundary);
    TopExp::MapShapes(b, TopAbs_EDGE, boundary);
    std::string lastError;
    for (int mask : {0, 2, 1, 3}) {
        auto orientedA = a, orientedB = b;
        if (mask & 1) orientedA.Reverse();
        if (mask & 2) orientedB.Reverse();
        BRepTools_ReShape reshape;
        const auto replace = [&](const TopoDS_Face& oldFace, TopoDS_Face newFace) {
            requireValidRepairCondition(oldFace.Orientation() == TopAbs_FORWARD || oldFace.Orientation() == TopAbs_REVERSED,
                "unsupported original face orientation");
            // FORWARD key 的 replacement 存储相对方向；ReShape 会合成原 face use 方向。
            if (oldFace.Orientation() == TopAbs_REVERSED) newFace.Reverse();
            reshape.Replace(oldFace.Oriented(TopAbs_FORWARD), newFace);
        };
        replace(c.a, orientedA);
        replace(c.b, orientedB);
        const auto result = reshape.Apply(trial);
        requireValidRepairCondition(!result.IsNull() && result.ShapeType() == trial.ShapeType(), "orientation replacement changed root type");
        TopTools_IndexedDataMapOfShapeListOfShape ef;
        TopExp::MapShapesAndUniqueAncestors(result, TopAbs_EDGE, TopAbs_FACE, ef);
        bool valid = true;
        for (int i = 1; i <= boundary.Extent(); ++i) {
            const auto& edge = boundary(i);
            if (!ef.Contains(edge)) { valid = false; lastError = "rebuilt boundary disappeared"; break; }
            const auto& owners = ef.FindFromKey(edge);
            if (owners.Extent() == 1 && !containsSameEdge(shared, edge)) continue;
            if (owners.Extent() != 2) { valid = false; lastError = "unsupported boundary face count"; break; }
            TopAbs_Orientation uses[2] = {TopAbs_INTERNAL, TopAbs_INTERNAL};
            int j = 0;
            for (const auto& owner : owners) {
                if (!edgeUse(TopoDS::Face(owner), edge, uses[j++])) valid = false;
            }
            if (!valid || uses[0] == uses[1]) {
                valid = false;
                lastError = containsSameEdge(shared, edge)
                    ? "shared edge directions agree" : "rebuilt face conflicts with an outer neighbor";
                break;
            }
        }
        if (!valid) {
            messages.push_back("Face orientation trial: flipA=" + std::to_string(mask & 1)
                + ", flipB=" + std::to_string((mask >> 1) & 1) + ", rejected: " + lastError);
            continue;
        }
        messages.push_back("Face orientation accepted: flipA=" + std::to_string(mask & 1)
            + ", flipB=" + std::to_string((mask >> 1) & 1)
            + "; shared and outer boundary directions verified in replacement shape");
        return result;
    }
    requireValidRepairCondition(false, "cannot orient both rebuilt faces consistently with shared edge and outer neighbors: " + lastError);
    return {};
}

TopoDS_Shape repairCandidateInIsolatedShape(const TopoDS_Shape& current, const Candidate& source, const Params& params,
                   std::vector<std::string>& messages, ProjectionQuality& quality) {
    BRepBuilderAPI_Copy copy(current, Standard_True, Standard_False);
    requireValidRepairCondition(!copy.Shape().IsNull(), "null trial copy");
    TopTools_IndexedMapOfShape originalVertices;
    TopExp::MapShapes(copy.Shape(), TopAbs_VERTEX, originalVertices);
    std::vector<double> vertexTolerances;
    for (int i = 1; i <= originalVertices.Extent(); ++i)
        vertexTolerances.push_back(BRep_Tool::Tolerance(TopoDS::Vertex(originalVertices(i))));
    Candidate c;
    c.p = TopoDS::Vertex(copy.ModifiedShape(source.p)); c.q = TopoDS::Vertex(copy.ModifiedShape(source.q));
    c.a = TopoDS::Face(copy.ModifiedShape(source.a)); c.b = TopoDS::Face(copy.ModifiedShape(source.b));
    // ModifiedShape 的映射不应取代原 face 在 shell 中的 use 方向。
    c.a.Orientation(source.a.Orientation()); c.b.Orientation(source.b.Orientation());
    for (const auto& e : source.sideA) c.sideA.push_back(TopoDS::Edge(copy.ModifiedShape(e).Oriented(e.Orientation())));
    for (const auto& e : source.sideB) c.sideB.push_back(TopoDS::Edge(copy.ModifiedShape(e).Oriented(e.Orientation())));
    {
        TopLoc_Location la, lb;
        const auto sa = BRep_Tool::Surface(c.a, la), sb = BRep_Tool::Surface(c.b, lb);
        std::ostringstream os;
        os << "Candidate supports: A=" << (sa.IsNull() ? "null" : sa->DynamicType()->Name())
           << ", B=" << (sb.IsNull() ? "null" : sb->DynamicType()->Name())
           << ", sameSurfaceHandle=" << (sa == sb) << ", sameLocation=" << la.IsEqual(lb);
        messages.push_back(os.str());
    }
    auto shared = buildSharedEdgeByProjectingPQ(c, params, messages, quality);
    const auto a = rebuildFaceWithSharedBoundary(c.a, c.sideA, shared, c, params);
    const auto b = rebuildFaceWithSharedBoundary(c.b, c.sideB, shared, c, params);
    quality.maxSurfaceDeviation = 0.0;
    for (const auto& face : {a, b}) for (const auto& edge : shared) {
        BRepLib_CheckCurveOnSurface check(edge, face);
        check.Perform();
        requireValidRepairCondition(check.IsDone() && std::isfinite(check.MaxDistance()),
            "cannot measure rebuilt shared-edge surface deviation");
        quality.maxSurfaceDeviation = std::max(quality.maxSurfaceDeviation, check.MaxDistance());
    }
    // 第二张面补 pcurve 后，重新检查两张面，防止共享边被后续修改破坏。
    requireValidRepairCondition(BRepCheck_Analyzer(a, Standard_True).IsValid() && BRepCheck_Analyzer(b, Standard_True).IsValid(),
        "shared edge update invalidated rebuilt face");
    const auto result = replaceFacesWithConsistentOrientation(copy.Shape(), c, a, b, shared, messages);
    requireValidRepairCondition(!result.IsNull() && result.ShapeType() == current.ShapeType(), "root type changed");
    TopTools_IndexedMapOfShape oldFaces, newFaces, outputEdges;
    TopExp::MapShapes(current, TopAbs_FACE, oldFaces); TopExp::MapShapes(result, TopAbs_FACE, newFaces);
    TopExp::MapShapes(result, TopAbs_EDGE, outputEdges);
    requireValidRepairCondition(oldFaces.Extent() == newFaces.Extent(), "face count changed");
    for (const auto& side : {c.sideA, c.sideB}) for (const auto& e : side)
        requireValidRepairCondition(!outputEdges.Contains(e), "incorrect boundary edge still present");
    TopTools_IndexedDataMapOfShapeListOfShape ef;
    TopExp::MapShapesAndUniqueAncestors(result, TopAbs_EDGE, TopAbs_FACE, ef);
    for (const auto& e : shared) {
        requireValidRepairCondition(ef.Contains(e) && ef.FindFromKey(e).Extent() == 2, "projected edge not shared by two faces");
        std::vector<TopAbs_Orientation> uses;
        for (const auto& f : ef.FindFromKey(e)) for (TopExp_Explorer it(f, TopAbs_EDGE); it.More(); it.Next())
            if (it.Current().IsSame(e)) uses.push_back(it.Current().Orientation());
        requireValidRepairCondition(uses.size() == 2 && uses[0] != uses[1], "shared edge use directions disagree");
        requireValidRepairCondition(std::isfinite(BRep_Tool::Tolerance(e)) && BRep_Tool::Tolerance(e) <= params.maxRebuiltTolerance,
            "final shared edge exceeds tolerance budget");
        for (const auto& f : {a, b}) {
            BRepLib_CheckCurveOnSurface check(e, f);
            check.Perform();
            requireValidRepairCondition(check.IsDone() && std::isfinite(check.MaxDistance())
                && check.MaxDistance() <= params.pcurveTolerance, "final shared pcurve mismatch");
        }
    }
    TopTools_IndexedMapOfShape resultVertices;
    TopExp::MapShapes(result, TopAbs_VERTEX, resultVertices);
    for (int i = 1; i <= resultVertices.Extent(); ++i) {
        const auto v = TopoDS::Vertex(resultVertices(i));
        const int old = originalVertices.FindIndex(v);
        const double tolerance = BRep_Tool::Tolerance(v);
        requireValidRepairCondition(std::isfinite(tolerance) && (tolerance <= params.maxRebuiltTolerance
            || (old > 0 && tolerance <= vertexTolerances[old - 1])),
            "new/increased vertex tolerance exceeds budget");
    }
    std::string error;
    const auto before = RepairDiagnostics::CollectBRepCheckFailures(current, &error);
    requireValidRepairCondition(error.empty(), "baseline check failed: " + error);
    const auto after = RepairDiagnostics::CollectBRepCheckFailures(result, &error);
    requireValidRepairCondition(error.empty(), "result check failed: " + error);
    const auto added = RepairDiagnostics::FindBRepCheckFailureCountIncrease(before, after);
    requireValidRepairCondition(added.empty(), "new BRepCheck errors: " + RepairDiagnostics::DescribeBRepCheckFailures(added));
    requireValidRepairCondition(TopologyUtils::CountFreeEdges(result) < TopologyUtils::CountFreeEdges(current), "free edges did not decrease");
    return result;
}
struct Result {
    TopoDS_Shape shape;
    int committed = 0;
    std::vector<std::string> messages;
};

Result perform(const TopoDS_Shape& shape, Params params) {
    params.maxGapWidth = normalizePositiveTolerance(params.maxGapWidth, 0.5);
    params.maxEndpointCollapseDistance = normalizePositiveTolerance(params.maxEndpointCollapseDistance, 0.8);
    params.endpointLengthDifferenceTolerance = normalizePositiveTolerance(params.endpointLengthDifferenceTolerance, 1.0e-5);
    params.pcurveTolerance = normalizePositiveTolerance(params.pcurveTolerance, 1.0e-5);
    params.maxRebuiltTolerance = normalizePositiveTolerance(params.maxRebuiltTolerance, 1.0e-4);
    params.maxRelativeAreaChange = std::min(0.5, normalizePositiveTolerance(params.maxRelativeAreaChange, 0.15));
    params.maxBoundaryEdges = std::clamp(params.maxBoundaryEdges, 3, 32);
    Result output; output.shape = shape;
    if (shape.IsNull()) { output.messages.push_back("Null input"); return output; }
    int committed = 0;
    try {
        // 成功后重新提取，避免引用过期的 face。每次消除一对端点间的自由边链，数量单调减少。
        while (true) {
            const auto candidates = findBilateralCuspCandidates(output.shape, params, output.messages);
            output.messages.push_back("Bilateral cusp candidates=" + std::to_string(candidates.size()));
            bool changed = false;
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                struct TrialOutcome {
                    TopoDS_Shape shape;
                    ProjectionQuality quality;
                    std::vector<std::string> messages;
                    std::string error;
                };
                std::array<TrialOutcome, 2> trials;
                for (int surfaceIndex = 0; surfaceIndex < 2; ++surfaceIndex) {
                    Candidate candidate = candidates[i];
                    if (surfaceIndex == 1) {
                        std::swap(candidate.a, candidate.b);
                        std::swap(candidate.sideA, candidate.sideB);
                    }
                    try {
                        trials[surfaceIndex].shape = repairCandidateInIsolatedShape(
                            output.shape, candidate, params, trials[surfaceIndex].messages,
                            trials[surfaceIndex].quality);
                        const ProjectionQuality& quality = trials[surfaceIndex].quality;
                        requireValidRepairCondition(std::isfinite(quality.meanSideImbalance)
                            && std::isfinite(quality.maxSurfaceDeviation)
                            && std::isfinite(quality.maxEndpointDeviation),
                            "projection quality is not finite");
                    } catch (const Standard_Failure& e) {
                        trials[surfaceIndex].shape.Nullify();
                        trials[surfaceIndex].error = e.GetMessageString()
                            ? e.GetMessageString() : "OCCT failure";
                    } catch (const std::exception& e) {
                        trials[surfaceIndex].shape.Nullify();
                        trials[surfaceIndex].error = e.what();
                    }
                }
                if (trials[0].shape.IsNull() && trials[1].shape.IsNull()) {
                    output.messages.push_back("Candidate rejected: projection on A: "
                        + trials[0].error + "; projection on B: " + trials[1].error);
                    continue;
                }
                const auto score = [&params](const ProjectionQuality& quality) {
                    // Compare dimensionless geometric risks; neither trial is
                    // accepted without full face and whole-shape validation.
                    return quality.meanSideImbalance / params.maxGapWidth
                        + quality.maxSurfaceDeviation / params.maxRebuiltTolerance
                        + quality.maxEndpointDeviation / params.maxRebuiltTolerance;
                };
                const int selected = trials[0].shape.IsNull() ? 1
                    : trials[1].shape.IsNull() ? 0
                    : (score(trials[1].quality) < score(trials[0].quality) ? 1 : 0);
                output.shape = trials[selected].shape;
                ++committed;
                changed = true;
                output.messages.insert(output.messages.end(), trials[selected].messages.begin(),
                    trials[selected].messages.end());
                std::ostringstream selection;
                selection.precision(9);
                const int validTrialCount =
                    static_cast<int>(!trials[0].shape.IsNull())
                    + static_cast<int>(!trials[1].shape.IsNull());
                selection << "Projection surface selected: " << (selected == 0 ? 'A' : 'B')
                          << ", valid trials=" << validTrialCount
                          << ", mean side imbalance="
                          << trials[selected].quality.meanSideImbalance
                          << ", max surface deviation="
                          << trials[selected].quality.maxSurfaceDeviation;
                output.messages.push_back(selection.str());
                output.messages.push_back("Committed candidate " + std::to_string(i + 1)
                    + ": original surfaces and existing anchor vertices preserved");
                break;
            }
            if (!changed) break;
        }
    } catch (const Standard_Failure& e) {
        output.messages.push_back(std::string("Detection stopped (OCCT): ") + (e.GetMessageString() ? e.GetMessageString() : "unknown"));
    } catch (const std::exception& e) { output.messages.push_back(std::string("Detection stopped: ") + e.what()); }
    output.messages.push_back("Bilateral cusp committed=" + std::to_string(committed));
    output.committed = committed;
    return output;
}
} // namespace bilateral

} // namespace

ShapePreProcess::ShapePreProcess(const Params& params)
    : params_(params)
{
    const auto positive = [](double value, double fallback) {
        return std::isfinite(value) && value >= Precision::Confusion()
            ? value : fallback;
    };
    params_.zeroDistanceTolerance = positive(params_.zeroDistanceTolerance, 1.0e-4);
    params_.maxShortBoundaryEdgeLength = positive(params_.maxShortBoundaryEdgeLength, 1.1e-1);
    params_.minLongToShortLengthRatio =
        std::isfinite(params_.minLongToShortLengthRatio) &&
        params_.minLongToShortLengthRatio >= 2.0
        ? params_.minLongToShortLengthRatio : 3.0;
    params_.maxBoundaryChange = positive(params_.maxBoundaryChange, 1.1e-1);
    params_.maxEdgeTolerance = positive(params_.maxEdgeTolerance, 1.0e-2);
    params_.maxBilateralSharedEdgeTolerance = positive(
        params_.maxBilateralSharedEdgeTolerance, 1.0e-1);
}

void ShapePreProcess::setInputShape(const TopoDS_Shape& shape)
{
    inputShape_ = shape;
    resultShape_ = inputShape_;
    modified_ = false;
    mergedCount_ = 0;
    collapsedBoundaryCount_ = 0;
    messages_.clear();
}

int ShapePreProcess::perform()
{
    modified_ = false;
    mergedCount_ = 0;
    collapsedBoundaryCount_ = 0;
    messages_.clear();
    resultShape_ = inputShape_;

    if (inputShape_.IsNull()) {
        appendMessage("Input shape is null");
        return 0;
    }

    TopoDS_Shape currentShape = inputShape_;
    if (params_.enableDuplicateVertexMerge) {
        const int iterationLimit = computeAutomaticIterationLimit(currentShape);
        for (int iter = 0; iter < iterationLimit; ++iter) {
            std::vector<Candidate> candidates = findDuplicateVertexCandidates(currentShape);
            std::sort(candidates.begin(), candidates.end(),
                [](const Candidate& lhs, const Candidate& rhs) {
                    return lhs.confidence > rhs.confidence;
                });
            bool accepted = false;
            for (const Candidate& candidate : candidates) {
                TopoDS_Shape mergedShape;
                std::string errorMessage;
                if (!tryMergeCandidate(currentShape, candidate, mergedShape, errorMessage)) {
                    if (params_.verbose && messages_.size() < 12) {
                        appendMessage("Vertex pair skipped: " + errorMessage);
                    }
                    continue;
                }
                currentShape = mergedShape;
                modified_ = true;
                ++mergedCount_;
                accepted = true;
                break;
            }
            if (!accepted) {
                break;
            }
        }
    }

    if (params_.enableShortBoundaryEdgeCollapse) {
        // Two symmetric short caps at each anchor mean neither face boundary
        // is a trustworthy one-sided replacement. Rebuild both faces around
        // one projected edge before scanning unilateral Zero/One Cusp rings.
        bilateral::Params bilateralParams;
        bilateralParams.maxGapWidth = params_.maxBoundaryChange;
        bilateralParams.maxEndpointCollapseDistance =
            params_.maxShortBoundaryEdgeLength;
        bilateralParams.endpointLengthDifferenceTolerance =
            kBilateralAnchorLengthDifferenceTolerance;
        bilateralParams.pcurveTolerance =
            params_.maxBilateralSharedEdgeTolerance;
        bilateralParams.maxRebuiltTolerance =
            params_.maxBilateralSharedEdgeTolerance;
        const bilateral::Result bilateralRepair =
            bilateral::perform(currentShape, bilateralParams);
        const int bilateralCount = bilateralRepair.committed;
        if (bilateralCount > 0) {
            std::string validationError;
            const TopoDS_Shape bilateralShape = bilateralRepair.shape;
            if (validateMergeResult(currentShape, bilateralShape,
                    params_.zeroDistanceTolerance, validationError)) {
                currentShape = bilateralShape;
                modified_ = true;
                collapsedBoundaryCount_ += bilateralCount;
                std::string projectionChoice;
                if (bilateralCount == 1) {
                    for (const std::string& message : bilateralRepair.messages) {
                        if (message.find("Projection surface selected:") == 0) {
                            projectionChoice = message;
                        }
                    }
                }
                appendMessage("Bilateral cusp: repaired cases="
                    + std::to_string(bilateralCount)
                    + (projectionChoice.empty() ? std::string()
                        : ", " + projectionChoice));
            }
            else if (params_.verbose) {
                appendMessage("Bilateral cusp skipped: " + validationError);
            }
        }
        else if (params_.verbose) {
            const auto& bilateralMessages = bilateralRepair.messages;
            const auto failure = std::find_if(bilateralMessages.rbegin(),
                bilateralMessages.rend(), [](const std::string& message) {
                    return message.find("Candidate rejected:") == 0;
                });
            if (failure != bilateralMessages.rend()) {
                appendMessage("Bilateral cusp skipped: " + *failure);
            }
        }
        const int iterationLimit = computeAutomaticIterationLimit(currentShape);
        for (int iter = 0; iter < iterationLimit; ++iter) {
            std::vector<BoundaryGapCandidate> candidates;
            try {
                candidates = findBoundaryGapCandidates(currentShape);
            }
            catch (const Standard_Failure& failure) {
                appendMessage("Boundary scan stopped: " + failureText(failure));
                break;
            }
            bool accepted = false;
            for (const BoundaryGapCandidate& candidate : candidates) {
                TopoDS_Shape repairedShape;
                std::string errorMessage;
                if (!tryRepairBoundaryGap(currentShape, candidate, repairedShape,
                        errorMessage)) {
                    if (params_.verbose && messages_.size() < 12) {
                        appendMessage("Boundary gap skipped: " + errorMessage);
                    }
                    continue;
                }
                currentShape = repairedShape;
                modified_ = true;
                ++collapsedBoundaryCount_;
                accepted = true;
                std::ostringstream os;
                os << (candidate.shortEdges.size() == 2 ? "Zero Cusp" : "One Cusp")
                   << ": repaired, side edges=" << candidate.keptPath.edges.size()
                   << "/" << candidate.movedPath.edges.size()
                   << ", width=" << candidate.maxWidth;
                appendMessage(os.str());
                break;
            }
            if (!accepted) {
                break;
            }
        }
    }

    resultShape_ = currentShape;
    std::ostringstream summary;
    summary << "PreProcess: merged vertex pairs=" << mergedCount_
            << ", short-edge cases=" << collapsedBoundaryCount_;
    appendMessage(summary.str());
    return mergedCount_ + collapsedBoundaryCount_;
}

int ShapePreProcess::perform(const TopoDS_Shape& shape)
{
    setInputShape(shape);
    return perform();
}

std::vector<ShapePreProcess::Candidate> ShapePreProcess::findDuplicateVertexCandidates(
    const TopoDS_Shape& shape)
{
    std::vector<Candidate> candidates;
    if (shape.IsNull()) {
        return candidates;
    }

    TopTools_IndexedMapOfShape vertexMap;
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertexMap);
    if (vertexMap.Extent() < 2) {
        return candidates;
    }

    TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap;
    TopologyUtils::BuildVertexEdgeMap(shape, vertexEdgeMap);

    struct RawPair {
        int firstIndex = 0;
        int secondIndex = 0;
        double distance = 0.0;
    };

    struct IndexedPoint {
        int vertexIndex = 0;
        gp_Pnt point;
    };

    std::vector<int> nearZeroPartnerCount(static_cast<std::size_t>(vertexMap.Extent()) + 1, 0);
    std::vector<RawPair> rawPairs;
    int multiVertexPairCount = 0;
    std::vector<IndexedPoint> sortedPoints;
    sortedPoints.reserve(static_cast<std::size_t>(vertexMap.Extent()));
    for (int i = 1; i <= vertexMap.Extent(); ++i) {
        const TopoDS_Vertex firstVertex = TopoDS::Vertex(vertexMap.FindKey(i));
        if (firstVertex.IsNull()) {
            continue;
        }
        const gp_Pnt point = BRep_Tool::Pnt(firstVertex);
        if (std::isfinite(point.X()) && std::isfinite(point.Y())
            && std::isfinite(point.Z())) {
            sortedPoints.push_back({ i, point });
        }
    }
    std::sort(sortedPoints.begin(), sortedPoints.end(),
        [](const IndexedPoint& lhs, const IndexedPoint& rhs) {
            return lhs.point.X() < rhs.point.X();
        });
    for (std::size_t i = 0; i < sortedPoints.size(); ++i) {
        for (std::size_t j = i + 1; j < sortedPoints.size()
             && sortedPoints[j].point.X() - sortedPoints[i].point.X()
                 <= params_.zeroDistanceTolerance; ++j) {
            if (std::abs(sortedPoints[j].point.Y() - sortedPoints[i].point.Y())
                    > params_.zeroDistanceTolerance
                || std::abs(sortedPoints[j].point.Z() - sortedPoints[i].point.Z())
                    > params_.zeroDistanceTolerance) {
                continue;
            }
            const double distance = sortedPoints[i].point.Distance(sortedPoints[j].point);
            if (distance > params_.zeroDistanceTolerance) {
                continue;
            }
            const int firstIndex = sortedPoints[i].vertexIndex;
            const int secondIndex = sortedPoints[j].vertexIndex;
            rawPairs.push_back({ firstIndex, secondIndex, distance });
            ++nearZeroPartnerCount[static_cast<std::size_t>(firstIndex)];
            ++nearZeroPartnerCount[static_cast<std::size_t>(secondIndex)];
        }
    }

    for (const RawPair& rawPair : rawPairs) {
        if (nearZeroPartnerCount[static_cast<std::size_t>(rawPair.firstIndex)] != 1 ||
            nearZeroPartnerCount[static_cast<std::size_t>(rawPair.secondIndex)] != 1) {
            ++multiVertexPairCount;
            continue;
        }

        const TopoDS_Vertex firstVertex = TopoDS::Vertex(vertexMap.FindKey(rawPair.firstIndex));
        const TopoDS_Vertex secondVertex = TopoDS::Vertex(vertexMap.FindKey(rawPair.secondIndex));
        const TopoDS_Vertex keptVertex = chooseKeptVertex(firstVertex, secondVertex, vertexEdgeMap);
        const TopoDS_Vertex duplicateVertex =
            keptVertex.IsSame(firstVertex) ? secondVertex : firstVertex;

        Candidate candidate;
        candidate.keptVertex = keptVertex;
        candidate.duplicateVertex = duplicateVertex;
        candidate.distance = rawPair.distance;
        candidate.confidence = 1.0 - std::min(1.0, rawPair.distance / params_.zeroDistanceTolerance);
        candidates.push_back(candidate);
    }

    if (multiVertexPairCount > 0 && params_.verbose) {
        std::ostringstream os;
        os << "Skipped near-zero pairs involved in 3+ vertex groups: " << multiVertexPairCount;
        appendMessage(os.str());
    }

    return candidates;
}


std::vector<ShapePreProcess::BoundaryGapCandidate>
ShapePreProcess::findBoundaryGapCandidates(const TopoDS_Shape& shape) const
{
    std::vector<BoundaryGapCandidate> candidates;
    if (shape.IsNull()) {
        return candidates;
    }

    struct BoundaryUse {
        TopoDS_Edge edge;
        int first = 0;
        int last = 0;
        int face = 0;
        double length = 0.0;
    };
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap;
    TopTools_IndexedDataMapOfShapeListOfShape faceShellMap;
    TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);
    TopologyUtils::BuildVertexEdgeMap(shape, vertexEdgeMap);
    TopologyUtils::BuildFaceShellMap(shape, faceShellMap);

    TopTools_IndexedMapOfShape vertices;
    TopTools_IndexedMapOfShape faces;
    std::vector<BoundaryUse> uses;
    for (int index = 1; index <= edgeFaceMap.Extent(); ++index) {
        if (edgeFaceMap.FindFromIndex(index).Extent() != 1) {
            continue;
        }
        const TopoDS_Edge edge = TopoDS::Edge(edgeFaceMap.FindKey(index));
        const TopoDS_Face face =
            TopoDS::Face(edgeFaceMap.FindFromIndex(index).First());
        const TopoDS_Vertex first = TopExp::FirstVertex(edge);
        const TopoDS_Vertex last = TopExp::LastVertex(edge);
        Standard_Real curveFirst = 0.0;
        Standard_Real curveLast = 0.0;
        if (first.IsNull() || last.IsNull() || first.IsSame(last)
            || BRep_Tool::Degenerated(edge) || BRep_Tool::IsClosed(edge, face)
            || BRep_Tool::Curve(edge, curveFirst, curveLast).IsNull()
            || !std::isfinite(curveFirst) || !std::isfinite(curveLast)
            || curveLast - curveFirst <= Precision::PConfusion()) {
            continue;
        }
        const double length = GeometryUtils::ComputeEdgeLength(edge);
        if (!std::isfinite(length) || length <= Precision::Confusion()) {
            continue;
        }
        uses.push_back({edge, vertices.Add(first), vertices.Add(last),
                        faces.Add(face), length});
    }
    if (uses.size() < 3) {
        return candidates;
    }

    std::vector<std::vector<int>> incident(
        static_cast<std::size_t>(vertices.Extent()) + 1);
    std::vector<int> shortUses;
    for (std::size_t index = 0; index < uses.size(); ++index) {
        const BoundaryUse& use = uses[index];
        incident[static_cast<std::size_t>(use.first)].push_back(
            static_cast<int>(index));
        incident[static_cast<std::size_t>(use.last)].push_back(
            static_cast<int>(index));
        if (use.length <= params_.maxShortBoundaryEdgeLength) {
            shortUses.push_back(static_cast<int>(index));
        }
    }

    // A pair of similar short free edges from different faces at each anchor
    // describes a bilateral defect, not a safe one-sided kept boundary.
    std::vector<bool> bilateralAnchor(
        static_cast<std::size_t>(vertices.Extent()) + 1, false);
    for (int vertex = 1; vertex <= vertices.Extent(); ++vertex) {
        const auto& adjacent = incident[static_cast<std::size_t>(vertex)];
        for (std::size_t i = 0; i < adjacent.size() && !bilateralAnchor[vertex]; ++i) {
            const BoundaryUse& first = uses[static_cast<std::size_t>(adjacent[i])];
            if (first.length > params_.maxShortBoundaryEdgeLength) {
                continue;
            }
            for (std::size_t j = i + 1; j < adjacent.size(); ++j) {
                const BoundaryUse& second = uses[static_cast<std::size_t>(adjacent[j])];
                if (first.face != second.face
                    && second.length <= params_.maxShortBoundaryEdgeLength
                    && std::abs(first.length - second.length)
                        <= kBilateralAnchorLengthDifferenceTolerance) {
                    bilateralAnchor[vertex] = true;
                    break;
                }
            }
        }
    }

    const auto vertexAt = [&vertices](int index) {
        return TopoDS::Vertex(vertices.FindKey(index));
    };
    const auto faceAt = [&faces](int index) {
        return TopoDS::Face(faces.FindKey(index));
    };
    const auto hasFaceEdge = [&incident, &uses](int vertex, int face) {
        for (const int edgeIndex : incident[static_cast<std::size_t>(vertex)]) {
            if (uses[static_cast<std::size_t>(edgeIndex)].face == face) {
                return true;
            }
        }
        return false;
    };
    const auto sameShell = [&faceShellMap, &faceAt](int first, int second) {
        const TopoDS_Face firstFace = faceAt(first);
        const TopoDS_Face secondFace = faceAt(second);
        return first != second && faceShellMap.Contains(firstFace)
            && faceShellMap.Contains(secondFace)
            && faceShellMap.FindFromKey(firstFace).Extent() == 1
            && faceShellMap.FindFromKey(secondFace).Extent() == 1
            && faceShellMap.FindFromKey(firstFace).First().IsSame(
                faceShellMap.FindFromKey(secondFace).First());
    };
    const auto preferredFixed = [&vertexAt, &vertexEdgeMap](
        int fixed, int moved) {
        const TopoDS_Vertex fixedVertex = vertexAt(fixed);
        const TopoDS_Vertex movedVertex = vertexAt(moved);
        const int fixedDegree = TopologyUtils::GetVertexDegree(
            fixedVertex, vertexEdgeMap);
        const int movedDegree = TopologyUtils::GetVertexDegree(
            movedVertex, vertexEdgeMap);
        return fixedDegree > movedDegree
            || (fixedDegree == movedDegree
                && BRep_Tool::Tolerance(fixedVertex) <=
                    BRep_Tool::Tolerance(movedVertex) + Precision::Confusion());
    };

    const auto shortestPath = [&uses, &incident, &vertexAt, &vertices](
        int face, int start, int finish, int blockedA, int blockedB) {
        BoundaryPath path;
        if (start == finish) {
            return path;
        }
        const double infinity = std::numeric_limits<double>::infinity();
        std::vector<double> distance(
            static_cast<std::size_t>(vertices.Extent()) + 1, infinity);
        std::vector<int> previousEdge(distance.size(), -1);
        std::vector<int> previousVertex(distance.size(), -1);
        using QueueItem = std::pair<double, int>;
        std::priority_queue<QueueItem, std::vector<QueueItem>,
                            std::greater<QueueItem>> queue;
        distance[static_cast<std::size_t>(start)] = 0.0;
        queue.push({0.0, start});
        while (!queue.empty()) {
            const auto [currentDistance, current] = queue.top();
            queue.pop();
            if (currentDistance > distance[static_cast<std::size_t>(current)]) {
                continue;
            }
            if (current == finish) {
                break;
            }
            for (const int edgeIndex : incident[static_cast<std::size_t>(current)]) {
                if (edgeIndex == blockedA || edgeIndex == blockedB) {
                    continue;
                }
                const BoundaryUse& use = uses[static_cast<std::size_t>(edgeIndex)];
                if (use.face != face) {
                    continue;
                }
                const int next = use.first == current ? use.last : use.first;
                const double nextDistance = currentDistance + use.length;
                if (nextDistance >= distance[static_cast<std::size_t>(next)]) {
                    continue;
                }
                distance[static_cast<std::size_t>(next)] = nextDistance;
                previousEdge[static_cast<std::size_t>(next)] = edgeIndex;
                previousVertex[static_cast<std::size_t>(next)] = current;
                queue.push({nextDistance, next});
            }
        }
        if (!std::isfinite(distance[static_cast<std::size_t>(finish)])) {
            return path;
        }
        std::vector<int> reversedEdges;
        std::vector<int> reversedVertices{finish};
        for (int current = finish; current != start;) {
            const int edgeIndex = previousEdge[static_cast<std::size_t>(current)];
            const int previous = previousVertex[static_cast<std::size_t>(current)];
            if (edgeIndex < 0 || previous < 0 || reversedEdges.size() >= 128) {
                return BoundaryPath();
            }
            reversedEdges.push_back(edgeIndex);
            reversedVertices.push_back(previous);
            current = previous;
        }
        path.length = distance[static_cast<std::size_t>(finish)];
        for (auto it = reversedVertices.rbegin(); it != reversedVertices.rend(); ++it) {
            path.vertices.push_back(vertexAt(*it));
        }
        for (std::size_t index = 0; index < reversedEdges.size(); ++index) {
            const int edgeIndex = reversedEdges[reversedEdges.size() - index - 1];
            const TopoDS_Edge oriented = orientEdgeFromTo(
                uses[static_cast<std::size_t>(edgeIndex)].edge,
                path.vertices[index], path.vertices[index + 1]);
            if (oriented.IsNull()) {
                return BoundaryPath();
            }
            path.edges.push_back(oriented);
        }
        return path;
    };

    const auto appendCandidate = [this, &candidates, &faceAt,
                                  &vertices, &bilateralAnchor](
        const BoundaryPath& kept,
        const BoundaryPath& moved,
        const std::vector<TopoDS_Edge>& caps,
        int fixedFaceIndex,
        int defectiveFaceIndex) {
        if (kept.edges.empty() || moved.edges.empty()
            || kept.vertices.size() != kept.edges.size() + 1
            || moved.vertices.size() != moved.edges.size() + 1
            || !pathIsSimple(kept.vertices) || !pathIsSimple(moved.vertices)) {
            return;
        }
        const int firstAnchor = vertices.FindIndex(kept.vertices.front());
        const int lastAnchor = vertices.FindIndex(kept.vertices.back());
        if (firstAnchor > 0 && lastAnchor > 0
            && bilateralAnchor[static_cast<std::size_t>(firstAnchor)]
            && bilateralAnchor[static_cast<std::size_t>(lastAnchor)]) {
            return;
        }
        const bool zeroCusp = caps.size() == 2;
        for (std::size_t i = 0; i < kept.vertices.size(); ++i) {
            for (std::size_t j = 0; j < moved.vertices.size(); ++j) {
                if (kept.vertices[i].IsSame(moved.vertices[j])
                    && !( !zeroCusp && i + 1 == kept.vertices.size()
                          && j + 1 == moved.vertices.size())) {
                    return;
                }
            }
        }
        const double maxCap = caps.size() == 2
            ? std::max(GeometryUtils::ComputeEdgeLength(caps[0]),
                       GeometryUtils::ComputeEdgeLength(caps[1]))
            : GeometryUtils::ComputeEdgeLength(caps[0]);
        const double longer = std::max(kept.length, moved.length);
        if (!std::isfinite(maxCap) || maxCap > params_.maxShortBoundaryEdgeLength
            || std::min(kept.length, moved.length) <
                params_.minLongToShortLengthRatio * maxCap
            || std::min(kept.length, moved.length) < 0.5 * longer) {
            return;
        }
        std::vector<TopoDS_Edge> removed = caps;
        removed.insert(removed.end(), moved.edges.begin(), moved.edges.end());
        const TopoDS_Face defectiveFace = faceAt(defectiveFaceIndex);
        const TopoDS_Face fixedFace = faceAt(fixedFaceIndex);
        if (!edgesShareFaceWire(defectiveFace, removed)
            || !edgesShareFaceWire(fixedFace, kept.edges)) {
            return;
        }
        const double width = maximumPathGap(
            kept.edges, moved.edges, params_.maxBoundaryChange);
        if (!std::isfinite(width) || width > params_.maxBoundaryChange) {
            return;
        }
        BoundaryGapCandidate candidate;
        candidate.keptPath = kept;
        candidate.movedPath = moved;
        candidate.shortEdges = caps;
        candidate.keptFace = fixedFace;
        candidate.defectiveFace = defectiveFace;
        candidate.maxWidth = width;
        candidates.push_back(std::move(candidate));
    };

    for (const int firstCapIndex : shortUses) {
        const BoundaryUse& firstCap = uses[static_cast<std::size_t>(firstCapIndex)];
        for (int direction = 0; direction < 2; ++direction) {
            const int fixedStart = direction == 0 ? firstCap.first : firstCap.last;
            const int movedStart = direction == 0 ? firstCap.last : firstCap.first;
            if (!preferredFixed(fixedStart, movedStart)
                || !hasFaceEdge(movedStart, firstCap.face)) {
                continue;
            }
            std::vector<int> fixedFaces;
            for (const int edgeIndex : incident[static_cast<std::size_t>(fixedStart)]) {
                const int faceIndex = uses[static_cast<std::size_t>(edgeIndex)].face;
                if (faceIndex == firstCap.face
                    || std::find(fixedFaces.begin(), fixedFaces.end(), faceIndex)
                        != fixedFaces.end()) {
                    continue;
                }
                fixedFaces.push_back(faceIndex);
            }
            for (const int fixedFaceIndex : fixedFaces) {
                if (!sameShell(fixedFaceIndex, firstCap.face)) {
                    continue;
                }
                for (const int secondCapIndex : shortUses) {
                    if (secondCapIndex <= firstCapIndex) {
                        continue;
                    }
                    const BoundaryUse& secondCap =
                        uses[static_cast<std::size_t>(secondCapIndex)];
                    if (secondCap.face != firstCap.face) {
                        continue;
                    }
                    for (int secondDirection = 0; secondDirection < 2;
                         ++secondDirection) {
                        const int fixedEnd = secondDirection == 0
                            ? secondCap.first : secondCap.last;
                        const int movedEnd = secondDirection == 0
                            ? secondCap.last : secondCap.first;
                        if (fixedEnd == fixedStart || movedEnd == movedStart
                            || !preferredFixed(fixedEnd, movedEnd)
                            || !hasFaceEdge(fixedEnd, fixedFaceIndex)
                            || !hasFaceEdge(movedEnd, firstCap.face)) {
                            continue;
                        }
                        const BoundaryPath kept = shortestPath(
                            fixedFaceIndex, fixedStart, fixedEnd,
                            firstCapIndex, secondCapIndex);
                        const BoundaryPath moved = shortestPath(
                            firstCap.face, movedStart, movedEnd,
                            firstCapIndex, secondCapIndex);
                        appendCandidate(kept, moved,
                            {firstCap.edge, secondCap.edge},
                            fixedFaceIndex, firstCap.face);
                    }
                }
                for (int tip = 1; tip <= vertices.Extent(); ++tip) {
                    if (tip == fixedStart || tip == movedStart
                        || !hasFaceEdge(tip, fixedFaceIndex)
                        || !hasFaceEdge(tip, firstCap.face)) {
                        continue;
                    }
                    const BoundaryPath kept = shortestPath(
                        fixedFaceIndex, fixedStart, tip, firstCapIndex, -1);
                    const BoundaryPath moved = shortestPath(
                        firstCap.face, movedStart, tip, firstCapIndex, -1);
                    appendCandidate(kept, moved, {firstCap.edge},
                        fixedFaceIndex, firstCap.face);
                }
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(),
        [](const BoundaryGapCandidate& lhs, const BoundaryGapCandidate& rhs) {
            if (lhs.shortEdges.size() != rhs.shortEdges.size()) {
                return lhs.shortEdges.size() > rhs.shortEdges.size();
            }
            return lhs.maxWidth < rhs.maxWidth;
        });
    return candidates;
}

bool ShapePreProcess::isZeroDistanceVertexPair(
    const TopoDS_Vertex& firstVertex,
    const TopoDS_Vertex& secondVertex,
    double& distance) const
{
    distance = std::numeric_limits<double>::max();
    if (firstVertex.IsNull() || secondVertex.IsNull() || firstVertex.IsSame(secondVertex)) {
        return false;
    }

    const gp_Pnt firstPoint = BRep_Tool::Pnt(firstVertex);
    const gp_Pnt secondPoint = BRep_Tool::Pnt(secondVertex);
    distance = firstPoint.Distance(secondPoint);
    return distance <= params_.zeroDistanceTolerance;
}

TopoDS_Vertex ShapePreProcess::chooseKeptVertex(
    const TopoDS_Vertex& firstVertex,
    const TopoDS_Vertex& secondVertex,
    const TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap) const
{
    const int firstDegree = TopologyUtils::GetVertexDegree(firstVertex, vertexEdgeMap);
    const int secondDegree = TopologyUtils::GetVertexDegree(secondVertex, vertexEdgeMap);
    if (firstDegree != secondDegree) {
        return firstDegree > secondDegree ? firstVertex : secondVertex;
    }

    const double firstTolerance = BRep_Tool::Tolerance(firstVertex);
    const double secondTolerance = BRep_Tool::Tolerance(secondVertex);
    if (std::abs(firstTolerance - secondTolerance) > Precision::Confusion()) {
        return firstTolerance <= secondTolerance ? firstVertex : secondVertex;
    }

    return firstVertex;
}

bool ShapePreProcess::isSafeToMergePair(
    const TopoDS_Shape& shape,
    const Candidate& candidate,
    const TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap,
    std::string& errorMessage) const
{
    if (shape.IsNull() || candidate.keptVertex.IsNull() || candidate.duplicateVertex.IsNull()) {
        errorMessage = "shape or vertex is null";
        return false;
    }

    if (candidate.keptVertex.IsSame(candidate.duplicateVertex)) {
        errorMessage = "candidate vertices are already the same topological vertex";
        return false;
    }

    double distance = 0.0;
    if (!isZeroDistanceVertexPair(candidate.keptVertex, candidate.duplicateVertex, distance)) {
        errorMessage = "candidate vertices are not within zero-distance tolerance";
        return false;
    }

    const auto checkAdjacentEdges = [&](const TopoDS_Vertex& sourceVertex) -> bool {
        if (!vertexEdgeMap.Contains(sourceVertex)) {
            return true;
        }

        const TopTools_ListOfShape& edges = vertexEdgeMap.FindFromKey(sourceVertex);
        for (TopTools_ListIteratorOfListOfShape it(edges); it.More(); it.Next()) {
            const TopoDS_Edge edge = TopoDS::Edge(it.Value());
            if (edge.IsNull()) {
                continue;
            }

            TopoDS_Vertex first = TopExp::FirstVertex(edge);
            TopoDS_Vertex last = TopExp::LastVertex(edge);
            if (first.IsSame(candidate.duplicateVertex)) {
                first = candidate.keptVertex;
            }
            if (last.IsSame(candidate.duplicateVertex)) {
                last = candidate.keptVertex;
            }

            if (first.IsSame(last)) {
                return false;
            }
        }

        return true;
    };

    if (!checkAdjacentEdges(candidate.keptVertex) ||
        !checkAdjacentEdges(candidate.duplicateVertex)) {
        errorMessage = "merge would collapse at least one adjacent edge";
        return false;
    }

    if (params_.forbidAcrossDifferentSolids) {
        TopTools_IndexedDataMapOfShapeListOfShape vertexSolidMap;
        TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_SOLID, vertexSolidMap);
        if (vertexSolidMap.Contains(candidate.keptVertex) &&
            vertexSolidMap.Contains(candidate.duplicateVertex)) {
            const TopTools_ListOfShape& keptSolids = vertexSolidMap.FindFromKey(candidate.keptVertex);
            const TopTools_ListOfShape& duplicateSolids = vertexSolidMap.FindFromKey(candidate.duplicateVertex);
            if (keptSolids.Extent() > 0 &&
                duplicateSolids.Extent() > 0 &&
                !listsShareShape(keptSolids, duplicateSolids)) {
                errorMessage = "candidate vertices belong to different solids";
                return false;
            }
        }
    }

    return true;
}

bool ShapePreProcess::tryMergeCandidate(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& mergedShape,
    std::string& errorMessage)
{
    mergedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "current shape is null";
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap;
    TopologyUtils::BuildVertexEdgeMap(currentShape, vertexEdgeMap);
    if (!isSafeToMergePair(currentShape, candidate, vertexEdgeMap, errorMessage)) {
        return false;
    }

    try {
        const gp_Pnt keptPoint = BRep_Tool::Pnt(candidate.keptVertex);
        const double mergedTolerance = std::max({
            Precision::Confusion(),
            params_.zeroDistanceTolerance,
            candidate.distance,
            BRep_Tool::Tolerance(candidate.keptVertex),
            BRep_Tool::Tolerance(candidate.duplicateVertex)
        });

        BRep_Builder builder;
        TopoDS_Vertex mergedVertex;
        builder.MakeVertex(mergedVertex, keptPoint, mergedTolerance);

        BRepTools_ReShape reshaper;
        reshaper.Replace(candidate.keptVertex, mergedVertex);
        reshaper.Replace(candidate.duplicateVertex, mergedVertex);

        TopoDS_Shape trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()) {
            errorMessage = "reshaper produced a null shape";
            return false;
        }

        if (!repairAffectedTopology(trialShape, candidate, errorMessage)) {
            return false;
        }

        if (!validateMergeResult(currentShape, trialShape, params_.zeroDistanceTolerance, errorMessage)) {
            return false;
        }

        mergedShape = trialShape;
        return true;
    }
    catch (Standard_Failure& failure) {
        errorMessage = "OCCT failure while merging duplicate vertices: " + failureText(failure);
        return false;
    }
}


bool ShapePreProcess::tryRepairBoundaryGap(
    const TopoDS_Shape& currentShape,
    const BoundaryGapCandidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    repairedShape.Nullify();
    if (currentShape.IsNull() || candidate.keptFace.IsNull()
        || candidate.defectiveFace.IsNull()
        || candidate.keptPath.edges.empty() || candidate.movedPath.edges.empty()
        || (candidate.shortEdges.size() != 1 && candidate.shortEdges.size() != 2)) {
        errorMessage = "incomplete boundary-gap candidate";
        return false;
    }

    try {
        BRepBuilderAPI_Copy copier(currentShape, Standard_True, Standard_False);
        TopoDS_Shape trialShape = copier.Shape();
        const auto mappedEdge = [&copier](const TopoDS_Edge& source) {
            const TopoDS_Shape mapped = copier.ModifiedShape(source);
            return mapped.IsNull() || mapped.ShapeType() != TopAbs_EDGE
                ? TopoDS_Edge() : TopoDS::Edge(mapped);
        };
        const auto mappedVertex = [&copier](const TopoDS_Vertex& source) {
            const TopoDS_Shape mapped = copier.ModifiedShape(source);
            return mapped.IsNull() || mapped.ShapeType() != TopAbs_VERTEX
                ? TopoDS_Vertex() : TopoDS::Vertex(mapped);
        };
        const auto mappedFace = [&copier](const TopoDS_Face& source) {
            const TopoDS_Shape mapped = copier.ModifiedShape(source);
            return mapped.IsNull() || mapped.ShapeType() != TopAbs_FACE
                ? TopoDS_Face() : TopoDS::Face(mapped);
        };
        const TopoDS_Face oldFace = mappedFace(candidate.defectiveFace);
        const TopoDS_Face keptFace = mappedFace(candidate.keptFace);
        if (trialShape.IsNull() || oldFace.IsNull() || keptFace.IsNull()) {
            errorMessage = "failed to isolate boundary-gap faces";
            return false;
        }

        std::vector<TopoDS_Edge> keptEdges;
        std::vector<TopoDS_Vertex> keptVertices;
        std::vector<TopoDS_Edge> removedEdges;
        for (const TopoDS_Edge& edge : candidate.keptPath.edges) {
            keptEdges.push_back(mappedEdge(edge));
        }
        for (const TopoDS_Vertex& vertex : candidate.keptPath.vertices) {
            keptVertices.push_back(mappedVertex(vertex));
        }
        for (const TopoDS_Edge& edge : candidate.shortEdges) {
            removedEdges.push_back(mappedEdge(edge));
        }
        for (const TopoDS_Edge& edge : candidate.movedPath.edges) {
            removedEdges.push_back(mappedEdge(edge));
        }
        if (keptVertices.size() != keptEdges.size() + 1
            || std::any_of(keptEdges.begin(), keptEdges.end(),
                    [](const TopoDS_Edge& edge) { return edge.IsNull(); })
            || std::any_of(keptVertices.begin(), keptVertices.end(),
                    [](const TopoDS_Vertex& vertex) { return vertex.IsNull(); })
            || std::any_of(removedEdges.begin(), removedEdges.end(),
                    [](const TopoDS_Edge& edge) { return edge.IsNull(); })) {
            errorMessage = "failed to map boundary-gap edges or vertices";
            return false;
        }

        TopoDS_Wire targetWire;
        int targetWireCount = 0;
        std::vector<int> removedUseCounts(removedEdges.size(), 0);
        TopoDS_Edge firstCapUse;
        for (TopoDS_Iterator wireIt(oldFace, Standard_False, Standard_False);
             wireIt.More(); wireIt.Next()) {
            if (wireIt.Value().ShapeType() != TopAbs_WIRE) {
                errorMessage = "defective face contains a non-wire child";
                return false;
            }
            const TopoDS_Wire wire = TopoDS::Wire(wireIt.Value());
            int foundInWire = 0;
            for (TopoDS_Iterator edgeIt(wire, Standard_False, Standard_False);
                 edgeIt.More(); edgeIt.Next()) {
                if (edgeIt.Value().ShapeType() != TopAbs_EDGE) {
                    errorMessage = "defective wire contains a non-edge child";
                    return false;
                }
                const TopoDS_Edge edgeUse = TopoDS::Edge(edgeIt.Value());
                for (std::size_t index = 0; index < removedEdges.size(); ++index) {
                    if (!edgeUse.IsSame(removedEdges[index])) {
                        continue;
                    }
                    ++removedUseCounts[index];
                    ++foundInWire;
                    if (index == 0) {
                        firstCapUse = edgeUse;
                    }
                }
            }
            if (foundInWire == static_cast<int>(removedEdges.size())) {
                targetWire = wire;
                ++targetWireCount;
            }
        }
        if (targetWireCount != 1 || firstCapUse.IsNull()
            || std::any_of(removedUseCounts.begin(), removedUseCounts.end(),
                    [](int count) { return count != 1; })) {
            errorMessage = "target boundary is not one unique face-wire interval";
            return false;
        }

        TopoDS_Vertex capFirst;
        TopoDS_Vertex capLast;
        TopExp::Vertices(firstCapUse, capFirst, capLast, Standard_True);
        const bool replacementForward = capFirst.IsSame(keptVertices.front());
        if (!replacementForward && !capLast.IsSame(keptVertices.front())) {
            errorMessage = "first short edge does not meet the fixed boundary";
            return false;
        }
        BRep_Builder builder;
        TopoDS_Wire rebuiltWire;
        builder.MakeWire(rebuiltWire);
        bool inserted = false;
        for (TopoDS_Iterator edgeIt(targetWire, Standard_False, Standard_False);
             edgeIt.More(); edgeIt.Next()) {
            const TopoDS_Edge edgeUse = TopoDS::Edge(edgeIt.Value());
            const bool removeUse = std::any_of(removedEdges.begin(), removedEdges.end(),
                [&edgeUse](const TopoDS_Edge& removed) {
                    return edgeUse.IsSame(removed);
                });
            if (!removeUse) {
                builder.Add(rebuiltWire, edgeUse);
                continue;
            }
            if (inserted) {
                continue;
            }
            inserted = true;
            for (std::size_t index = 0; index < keptEdges.size(); ++index) {
                const std::size_t pathIndex = replacementForward
                    ? index : keptEdges.size() - index - 1;
                const TopoDS_Vertex from = replacementForward
                    ? keptVertices[pathIndex] : keptVertices[pathIndex + 1];
                const TopoDS_Vertex to = replacementForward
                    ? keptVertices[pathIndex + 1] : keptVertices[pathIndex];
                const TopoDS_Edge edgeUse = orientEdgeFromTo(
                    keptEdges[pathIndex], from, to);
                if (edgeUse.IsNull()) {
                    errorMessage = "kept chain orientation is inconsistent";
                    return false;
                }
                builder.Add(rebuiltWire, edgeUse);
            }
        }
        if (!inserted) {
            errorMessage = "target wire has no removable edge use";
            return false;
        }
        rebuiltWire.Orientation(targetWire.Orientation());

        TopLoc_Location surfaceLocation;
        const Handle(Geom_Surface) surface = BRep_Tool::Surface(
            oldFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "defective face has no support surface";
            return false;
        }
        TopoDS_Face rebuiltFace;
        builder.MakeFace(rebuiltFace, surface, surfaceLocation,
            BRep_Tool::Tolerance(oldFace));
        builder.NaturalRestriction(rebuiltFace,
            BRep_Tool::NaturalRestriction(oldFace));
        for (TopoDS_Iterator wireIt(oldFace, Standard_False, Standard_False);
             wireIt.More(); wireIt.Next()) {
            const TopoDS_Wire oldWire = TopoDS::Wire(wireIt.Value());
            builder.Add(rebuiltFace,
                oldWire.IsSame(targetWire) ? rebuiltWire : oldWire);
        }
        rebuiltFace.Orientation(oldFace.Orientation());

        for (TopoDS_Edge& edge : keptEdges) {
            std::string pcurveError;
            if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                    edge, rebuiltFace, params_.maxEdgeTolerance,
                    &pcurveError)) {
                errorMessage = "kept chain pcurve failed: " + pcurveError;
                return false;
            }
            BRepLib_CheckCurveOnSurface checker(edge, rebuiltFace);
            checker.Perform();
            std::string keptPCurveError;
            if (!checker.IsDone() || !std::isfinite(checker.MaxDistance())
                || checker.MaxDistance() > params_.maxEdgeTolerance
                || BRep_Tool::Tolerance(edge) > params_.maxEdgeTolerance
                || !GeometryUtils::HasValidPCurveOnFace(
                    edge, keptFace, params_.maxEdgeTolerance,
                    &keptPCurveError)) {
                errorMessage = "shared edge exceeds surface tolerance or lost its original pcurve";
                return false;
            }
        }

        BRepTools_ReShape reshaper;
        reshaper.Replace(oldFace, rebuiltFace);
        trialShape = reshaper.Apply(trialShape);
        if (trialShape.IsNull() || trialShape.ShapeType() != currentShape.ShapeType()) {
            errorMessage = "face replacement changed the root shape";
            return false;
        }
        TopTools_IndexedDataMapOfShapeListOfShape trialEdgeFaces;
        TopologyUtils::BuildEdgeFaceMap(trialShape, trialEdgeFaces);
        for (const TopoDS_Edge& edge : removedEdges) {
            if (trialEdgeFaces.Contains(edge)) {
                errorMessage = "removed boundary edge remains in the trial";
                return false;
            }
        }
        for (const TopoDS_Edge& edge : keptEdges) {
            if (!trialEdgeFaces.Contains(edge)
                || trialEdgeFaces.FindFromKey(edge).Extent() != 2) {
                errorMessage = "kept chain is not shared by two faces";
                return false;
            }
            TopoDS_Face otherFace;
            bool foundKeptFace = false;
            for (TopTools_ListIteratorOfListOfShape it(
                     trialEdgeFaces.FindFromKey(edge)); it.More(); it.Next()) {
                if (it.Value().IsSame(keptFace)) {
                    foundKeptFace = true;
                }
                else {
                    otherFace = TopoDS::Face(it.Value());
                }
            }
            TopAbs_Orientation keptOrientation = TopAbs_EXTERNAL;
            TopAbs_Orientation rebuiltOrientation = TopAbs_EXTERNAL;
            const bool keptUse = uniqueEdgeUseOrientation(
                keptFace, edge, keptOrientation);
            const bool rebuiltUse = !otherFace.IsNull()
                && uniqueEdgeUseOrientation(otherFace, edge, rebuiltOrientation);
            if (!foundKeptFace || otherFace.IsNull()
                || !keptUse || !rebuiltUse) {
                errorMessage = "shared edge use is missing or invalid";
                return false;
            }
            // An imported shell can already be unorientable. Keep the existing
            // face orientations and let the whole-shape non-regression check
            // decide whether this shared use worsens that baseline.
        }

        if (!validateMergeResult(
                currentShape, trialShape, params_.zeroDistanceTolerance,
                errorMessage)) {
            return false;
        }
        repairedShape = trialShape;
        return true;
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure during boundary-gap repair: " + failureText(failure);
        return false;
    }
}

bool ShapePreProcess::repairAffectedTopology(
    TopoDS_Shape& shape,
    const Candidate&,
    std::string& errorMessage) const
{
    if (shape.IsNull()) {
        errorMessage = "shape is null after vertex replacement";
        return false;
    }

    if (!params_.runShapeFixAfterMerge) {
        return true;
    }

    try {
        ShapeFix_Shape fixer(shape);
        fixer.SetPrecision(std::max(Precision::Confusion(), params_.zeroDistanceTolerance));
        fixer.SetMaxTolerance(std::max(Precision::Confusion(), params_.zeroDistanceTolerance * 10.0));
        fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
        fixer.FixFaceTool()->FixMissingSeamMode() = 1;
        fixer.Perform();

        const TopoDS_Shape fixedShape = fixer.Shape();
        if (!fixedShape.IsNull()) {
            shape = fixedShape;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "ShapeFix failed after duplicate vertex merge: " + failureText(failure);
        return false;
    }

    return true;
}

bool ShapePreProcess::validateMergeResult(
    const TopoDS_Shape& beforeShape,
    const TopoDS_Shape& afterShape,
    double collapsedEdgeTolerance,
    std::string& errorMessage) const
{
    if (beforeShape.IsNull() || afterShape.IsNull()) {
        errorMessage = "before or after shape is null";
        return false;
    }

    if (beforeShape.ShapeType() != afterShape.ShapeType()) {
        errorMessage = "shape type changed after preprocessing";
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(beforeShape);
    const int freeEdgesAfter = TopologyUtils::CountFreeEdges(afterShape);
    if (freeEdgesAfter > freeEdgesBefore) {
        errorMessage = "free edge count increased after preprocessing";
        return false;
    }

    const int collapsedEdgesBefore = countCollapsedEdges(beforeShape, collapsedEdgeTolerance);
    const int collapsedEdgesAfter = countCollapsedEdges(afterShape, collapsedEdgeTolerance);
    if (collapsedEdgesAfter > collapsedEdgesBefore) {
        errorMessage = "degenerate edge count increased after preprocessing";
        return false;
    }

    std::string baselineError;
    const std::vector<BRepCheckFailureInfo> baselineFailures =
        collectBRepCheckFailures(beforeShape, &baselineError);
    if (!baselineError.empty()) {
        errorMessage = "baseline BRepCheck failed: " + baselineError;
        return false;
    }

    std::string afterError;
    const std::vector<BRepCheckFailureInfo> afterFailures =
        collectBRepCheckFailures(afterShape, &afterError);
    if (!afterError.empty()) {
        errorMessage = "result BRepCheck failed: " + afterError;
        return false;
    }

    std::vector<BRepCheckFailureInfo> newFailures;
    for (const BRepCheckFailureInfo& failure : afterFailures) {
        if (!isExistingBRepCheckFailure(failure, baselineFailures)) {
            newFailures.push_back(failure);
        }
    }

    if (!newFailures.empty()) {
        // 顶点替换会重建相邻 face 的 TShape，身份变化不能证明错误是新增的。
        // 保持与流水线一致，按拓扑类型/错误状态的数量检查是否恶化。
        const auto toDiagnostics = [](const std::vector<BRepCheckFailureInfo>& failures) {
            RepairDiagnostics::FailureList result;
            for (const auto& failure : failures) {
                result.push_back({ failure.shape, failure.type, failure.statuses });
            }
            return result;
        };
        const auto increased = RepairDiagnostics::FindBRepCheckFailureCountIncrease(
            toDiagnostics(baselineFailures), toDiagnostics(afterFailures));
        if (increased.empty()) {
            return true;
        }
        errorMessage = "new BRepCheck failures after preprocessing:";
        errorMessage += RepairDiagnostics::DescribeBRepCheckFailures(increased);
        if (!baselineFailures.empty()) {
            errorMessage += " baseline failures:";
            errorMessage += describeBRepCheckFailures(baselineFailures);
        }
        return false;
    }

    return true;
}

int ShapePreProcess::computeAutomaticIterationLimit(const TopoDS_Shape& shape) const
{
    if (shape.IsNull()) {
        return 1;
    }

    TopTools_IndexedMapOfShape vertexMap;
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertexMap);
    if (vertexMap.Extent() <= 0) {
        return 1;
    }

    return std::min(vertexMap.Extent(), kMaxAutomaticMergeIterations);
}

int ShapePreProcess::countCollapsedEdges(const TopoDS_Shape& shape, double tolerance) const
{
    if (shape.IsNull()) {
        return 0;
    }

    TopTools_IndexedMapOfShape edgeMap;
    TopExp::MapShapes(shape, TopAbs_EDGE, edgeMap);

    int collapsedCount = 0;
    const double safeTolerance = std::max(tolerance, Precision::Confusion());
    for (int i = 1; i <= edgeMap.Extent(); ++i) {
        const TopoDS_Edge edge = TopoDS::Edge(edgeMap.FindKey(i));
        if (edge.IsNull()) {
            continue;
        }

        const TopoDS_Vertex firstVertex = TopExp::FirstVertex(edge);
        const TopoDS_Vertex lastVertex = TopExp::LastVertex(edge);
        if (!firstVertex.IsNull() && !lastVertex.IsNull() && firstVertex.IsSame(lastVertex)) {
            ++collapsedCount;
            continue;
        }

        const double length = GeometryUtils::ComputeEdgeLength(edge);
        if (length <= safeTolerance) {
            ++collapsedCount;
        }
    }

    return collapsedCount;
}

void ShapePreProcess::appendMessage(const std::string& message)
{
    messages_.push_back(message);
}
