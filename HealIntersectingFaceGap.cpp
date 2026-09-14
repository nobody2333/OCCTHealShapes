#include "HealIntersectingFaceGap.h"
#include "GeometryUtils.h"
#include "TopologyUtils.h"
#include "RepairDiagnostics.h"
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Wire.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRepTools.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GC_MakeSegment.hxx>
#include <GeomProjLib.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <set>
#include <stdexcept>

namespace {
using Params = HealIntersectingFaceGap::Params;
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
std::vector<Candidate> findIntersectingFaceGapCandidates(const TopoDS_Shape& shape, const Params& params,
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
// 固定使用候选 A 的原始支撑面；只投影一次，不做双面择优或曲面求交。
Edges buildSharedEdgeByProjectingPQ(const Candidate& c, const Params& params,
                          std::vector<std::string>& messages) {
    const auto surface = BRep_Tool::Surface(c.a); // 世界坐标，包含原 Location。
    requireValidRepairCondition(!surface.IsNull(), "projection surface A is null");
    requireValidRepairCondition(!surface->IsUPeriodic() && !surface->IsVPeriodic(), "periodic projection surface unsupported");
    const auto p = BRep_Tool::Pnt(c.p), q = BRep_Tool::Pnt(c.q);
    requireValidRepairCondition(p.Distance(q) > Precision::Confusion(), "P-Q segment has zero length");
    GC_MakeSegment segment(p, q);
    requireValidRepairCondition(segment.IsDone(), "cannot construct P-Q line segment");
    messages.push_back("Shared boundary method: project P-Q segment onto original surface A only");
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
    forEachBoundarySample(shared, [&](const gp_Pnt& point) {
        requireValidRepairCondition(computePointToBoundaryDistance(point, c.sideA) <= params.maxGapWidth
            && computePointToBoundaryDistance(point, c.sideB) <= params.maxGapWidth, "projected edge leaves defect band");
    });
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
        requireValidRepairCondition(needed <= params.maxRebuiltTolerance, "shared edge tolerance exceeds budget");
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
                   std::vector<std::string>& messages) {
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
    auto shared = buildSharedEdgeByProjectingPQ(c, params, messages);
    const auto a = rebuildFaceWithSharedBoundary(c.a, c.sideA, shared, c, params);
    const auto b = rebuildFaceWithSharedBoundary(c.b, c.sideB, shared, c, params);
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
} // namespace

HealIntersectingFaceGap::HealIntersectingFaceGap() : HealIntersectingFaceGap(Params{}) {}
HealIntersectingFaceGap::HealIntersectingFaceGap(const Params& params) : params_(params) {
    params_.maxGapWidth = normalizePositiveTolerance(params_.maxGapWidth, 0.5);
    params_.maxEndpointCollapseDistance = normalizePositiveTolerance(params_.maxEndpointCollapseDistance, 0.8);
    params_.endpointLengthDifferenceTolerance = normalizePositiveTolerance(params_.endpointLengthDifferenceTolerance, 1.0e-5);
    params_.pcurveTolerance = normalizePositiveTolerance(params_.pcurveTolerance, 1.0e-5);
    params_.maxRebuiltTolerance = normalizePositiveTolerance(params_.maxRebuiltTolerance, 1.0e-4);
    params_.maxRelativeAreaChange = std::min(0.5, normalizePositiveTolerance(params_.maxRelativeAreaChange, 0.15));
    params_.maxBoundaryEdges = std::clamp(params_.maxBoundaryEdges, 3, 32);
}
int HealIntersectingFaceGap::perform(const TopoDS_Shape& shape) {
    result_ = shape; modified_ = false; messages_.clear();
    if (shape.IsNull()) { messages_.push_back("Null input"); return 0; }
    int committed = 0;
    try {
        // 成功后重新提取，避免引用过期的 face。每次消除一对端点间的自由边链，数量单调减少。
        while (true) {
            const auto candidates = findIntersectingFaceGapCandidates(result_, params_, messages_);
            messages_.push_back("HealIntersectingFaceGap candidates=" + std::to_string(candidates.size()));
            bool changed = false;
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                try {
                    const auto trial = repairCandidateInIsolatedShape(result_, candidates[i], params_, messages_);
                    result_ = trial; ++committed; modified_ = true; changed = true;
                    messages_.push_back("Committed candidate " + std::to_string(i + 1)
                        + ": original surfaces and existing anchor vertices preserved");
                    break;
                } catch (const Standard_Failure& e) {
                    messages_.push_back(std::string("Candidate rejected (OCCT): ") + (e.GetMessageString() ? e.GetMessageString() : "unknown"));
                } catch (const std::exception& e) { messages_.push_back(std::string("Candidate rejected: ") + e.what()); }
            }
            if (!changed) break;
        }
    } catch (const Standard_Failure& e) {
        messages_.push_back(std::string("Detection stopped (OCCT): ") + (e.GetMessageString() ? e.GetMessageString() : "unknown"));
    } catch (const std::exception& e) { messages_.push_back(std::string("Detection stopped: ") + e.what()); }
    messages_.push_back("HealIntersectingFaceGap committed=" + std::to_string(committed));
    return committed;
}
