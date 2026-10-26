#include "GeometryUtils.h"

#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRepBndLib.hxx>
#include <BRepLib.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <Extrema_ExtPC.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <Bnd_Box.hxx>
#include <Precision.hxx>
#include <ShapeFix_Edge.hxx>
#include <Standard_Failure.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <TColStd_HArray1OfReal.hxx>

#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <map>
#include <iterator>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {
// 仅在常规修复失败后使用；所有试算均在独立 TShape 上完成。
bool rebuildPCurveByProjection(TopoDS_Edge& edge,
                               const TopoDS_Face& face,
                               double tolerance,
                               std::string& error)
{
    try {
        if (BRep_Tool::Degenerated(edge) || BRep_Tool::IsClosed(edge, face)) {
            error = "sampled pcurve rebuild does not support degenerate or seam edges";
            return false;
        }
        // 此回退只重建当前 face 的表示，不为其他 face 的 pcurve 宣称同参数成立。
        if (!BRep_Tool::SameParameter(edge) || !BRep_Tool::SameRange(edge)) {
            error = "sampled pcurve rebuild requires synchronized edge parameters";
            return false;
        }
        Standard_Real first = 0.0, last = 0.0;
        // 两个不带 location 输出参数的接口均返回世界坐标系中的几何。
        const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);
        const Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
        if (curve.IsNull() || surface.IsNull() ||
            !std::isfinite(first) || !std::isfinite(last) ||
            last - first <= Precision::PConfusion()) {
            error = "sampled pcurve rebuild has invalid geometry or parameter range";
            return false;
        }
        // 周期接缝和闭曲面的 UV 分支需要专门展开，不在此回退中猜测。
        if (surface->IsUClosed() || surface->IsVClosed() ||
            surface->IsUPeriodic() || surface->IsVPeriodic()) {
            error = "sampled pcurve rebuild does not support closed or periodic surfaces";
            return false;
        }

        constexpr int initialIntervals = 128;
        constexpr std::size_t maxSamples = 2049;
        constexpr int maxIterations = 16;
        const double parameterEpsilon = std::max(
            Precision::PConfusion(), (last - first) * 1.0e-12);
        const double projectionPrecision = Precision::Confusion();
        double uFirst, uLast, vFirst, vLast;
        surface->Bounds(uFirst, uLast, vFirst, vLast);
        GeomAPI_ProjectPointOnSurf projector;
        projector.Init(surface, uFirst, uLast, vFirst, vLast, projectionPrecision);
        std::map<double, gp_Pnt2d> samples;
        double maxProjectionDistance = 0.0;
        auto addSample = [&](double parameter) {
            auto next = samples.lower_bound(parameter);
            if ((next != samples.end() && next->first - parameter <= parameterEpsilon) ||
                (next != samples.begin() && parameter - std::prev(next)->first <= parameterEpsilon)) {
                return true;
            }
            if (samples.size() >= maxSamples) {
                error = "sampled pcurve rebuild reached sample limit";
                return false;
            }
            const gp_Pnt point = curve->Value(parameter);
            projector.Perform(point);
            if (!projector.IsDone() || projector.NbPoints() == 0) {
                error = "sampled pcurve rebuild could not project a curve point";
                return false;
            }
            Standard_Real u, v;
            projector.LowerDistanceParameters(u, v);
            const double distance = point.Distance(surface->Value(u, v));
            if (!std::isfinite(u) || !std::isfinite(v) || !std::isfinite(distance) ||
                distance > tolerance) {
                error = "sampled pcurve rebuild exceeds curve-to-surface tolerance";
                return false;
            }
            maxProjectionDistance = std::max(maxProjectionDistance, distance);
            samples.emplace_hint(next, parameter, gp_Pnt2d(u, v));
            return true;
        };
        for (int i = 0; i <= initialIntervals; ++i) {
            if (!addSample(i == initialIntervals ? last :
                    first + (last - first) * i / initialIntervals)) {
                return false;
            }
        }

        // EmptyCopied 克隆 edge 的曲线表示及其范围，不复制/修改原曲线几何。
        TopoDS_Edge probe = TopoDS::Edge(edge.EmptyCopied());
        probe.Orientation(TopAbs_FORWARD);
        BRep_Builder builder;
        double maxDeviation = 0.0;
        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            const int count = static_cast<int>(samples.size());
            Handle(TColgp_HArray1OfPnt2d) points = new TColgp_HArray1OfPnt2d(1, count);
            Handle(TColStd_HArray1OfReal) parameters = new TColStd_HArray1OfReal(1, count);
            int index = 1;
            for (const auto& sample : samples) {
                points->SetValue(index, sample.second);
                parameters->SetValue(index++, sample.first);
            }
            Geom2dAPI_Interpolate interpolator(points, parameters, Standard_False,
                Precision::PConfusion());
            interpolator.Perform();
            if (!interpolator.IsDone()) {
                error = "sampled pcurve interpolation failed";
                return false;
            }
            const Handle(Geom2d_Curve) pcurve = interpolator.Curve();
            builder.UpdateEdge(probe, pcurve, face, 0.0);
            builder.Range(probe, face, first, last);
            BRepLib_CheckCurveOnSurface checker(probe, face);
            checker.Perform();
            if (!checker.IsDone() || !std::isfinite(checker.MaxDistance())) {
                error = "sampled pcurve rebuild failed deviation check";
                return false;
            }
            maxDeviation = checker.MaxDistance();
            // 在既有验收预算内，再限制插值引入的额外偏差；不直接花完宽松容差。
            const double fitTolerance = std::min(tolerance,
                maxProjectionDistance + std::max(Precision::Confusion(), tolerance * 0.05));
            if (maxDeviation <= fitTolerance) {
                if (!GeometryUtils::HasValidPCurveOnFace(probe, face, tolerance, &error)) {
                    return false;
                }
                const double measuredTolerance = std::min(tolerance,
                    std::max(Precision::Confusion(), maxDeviation * (1.0 + 1.0e-6)));
                // 只提交当前 face 的 pcurve，保留 canonical edge、3D 曲线和其他表示。
                builder.UpdateEdge(edge, pcurve, face, measuredTolerance);
                builder.Range(edge, face, first, last);
                ShapeFix_Edge().FixVertexTolerance(edge, face);
                return true;
            }

            std::vector<double> additional;
            for (auto it = samples.begin(), next = std::next(it); next != samples.end(); ++it, ++next) {
                for (int part = 1; part < 4; ++part) {
                    const double t = it->first + (next->first - it->first) * part / 4.0;
                    const gp_Pnt2d uv = pcurve->Value(t);
                    const double distance = curve->Value(t).Distance(surface->Value(uv.X(), uv.Y()));
                    if (!std::isfinite(distance) || distance > fitTolerance) {
                        additional.push_back(t);
                    }
                }
            }
            const double worstParameter = checker.MaxParameter();
            if (std::isfinite(worstParameter) && worstParameter > first && worstParameter < last) {
                additional.push_back(worstParameter);
            }
            const std::size_t previousSize = samples.size();
            std::sort(additional.begin(), additional.end());
            for (double t : additional) {
                if (!addSample(t)) {
                    return false;
                }
            }
            if (samples.size() == previousSize) {
                break;
            }
        }
        std::ostringstream message;
        message << "sampled pcurve rebuild did not converge: maxDeviation=" << maxDeviation
                << ", tolerance=" << tolerance << ", samples=" << samples.size();
        error = message.str();
    }
    catch (const Standard_Failure& failure) {
        error = "OCCT failure during sampled pcurve rebuild: ";
        error += failure.GetMessageString() ? failure.GetMessageString() : "OCCT failure";
    }
    return false;
}
} // namespace

double GeometryUtils::ComputeEdgeLength(const TopoDS_Edge& edge)
{
    if (edge.IsNull()) {
        return 0.0;
    }

    try {
        BRepAdaptor_Curve curveAdaptor(edge);
        return GCPnts_AbscissaPoint::Length(curveAdaptor, curveAdaptor.FirstParameter(), curveAdaptor.LastParameter());
    }
    catch (Standard_Failure&) {
        // 退化边或无几何边，回退到两点距离
        TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
        TopoDS_Vertex v2 = TopExp::LastVertex(edge);
        if (v1.IsNull() || v2.IsNull()) 
            return 0.0;
        return BRep_Tool::Pnt(v1).Distance(BRep_Tool::Pnt(v2));
    }
}

double GeometryUtils::ComputeIntersectionAngle(const TopoDS_Edge& e1, const TopoDS_Edge& e2, const TopoDS_Vertex commonV)
{
    if (e1.IsNull() || e2.IsNull() || commonV.IsNull()) 
        return -1.0;

    try {
        // 使用 BRepAdaptor_Curve 会自动处理 TopoDS_Edge的Orientation和Location
        BRepAdaptor_Curve bac1(e1);
        BRepAdaptor_Curve bac2(e2);

        gp_Pnt pnt1, pnt2;
        gp_Vec d1, d2;

        // 判断公共顶点在 e1 上的位置，并获取指向外部的切向量
        if (commonV.IsSame(TopExp::FirstVertex(e1))) {
            bac1.D1(bac1.FirstParameter(), pnt1, d1);
            d1.Reverse(); // 在起点，导数方向指向内部，反转使其向外
        }
        else if (commonV.IsSame(TopExp::LastVertex(e1))) {
            bac1.D1(bac1.LastParameter(), pnt1, d1);
            // 在终点，导数方向顺着曲线向外，无需反转
        }
        else {
            return -1.0;
        }

        // 判断公共顶点在 e2 上的位置，并获取指向外部的切向量
        if (commonV.IsSame(TopExp::FirstVertex(e2))) {
            bac2.D1(bac2.FirstParameter(), pnt2, d2);
            d2.Reverse();
        }
        else if (commonV.IsSame(TopExp::LastVertex(e2))) {
            bac2.D1(bac2.LastParameter(), pnt2, d2);
        }
        else {
            return -1.0;
        }

        // 检查向量有效性
        if (d1.Magnitude() <= Precision::Confusion() || d2.Magnitude() <= Precision::Confusion()) {
            return -1.0;
        }

        return d1.Angle(d2);
    }
    catch (Standard_Failure&) {
        return -1.0;
    }
}

double GeometryUtils::ComputeBoundingBoxDiagonal(const TopoDS_Shape& shape)
{
    if (shape.IsNull()) {
        return 0.0;
    }

    Bnd_Box bbox;
    BRepBndLib::Add(shape, bbox);

    if (bbox.IsVoid()) {
        return 0.0;
    }

    gp_Pnt minCorner = bbox.CornerMin();
    gp_Pnt maxCorner = bbox.CornerMax();

    return minCorner.Distance(maxCorner);
}

double GeometryUtils::ComputeBoundingBoxDiagonal(const std::vector<TopoDS_Shape>& shapes)
{
    if (shapes.empty()) {
        return 0.0;
    }

    Bnd_Box bbox;
    for (const auto& shape : shapes) {
        if (!shape.IsNull()) {
            BRepBndLib::Add(shape, bbox);
        }
    }

    if (bbox.IsVoid()) {
        return 0.0;
    }

    gp_Pnt minCorner = bbox.CornerMin();
    gp_Pnt maxCorner = bbox.CornerMax();

    return minCorner.Distance(maxCorner);
}

bool GeometryUtils::HasValidPCurveOnFace(const TopoDS_Edge& edge, const TopoDS_Face& face)
{
    double tolerance = Precision::Confusion();
    try {
        if (!edge.IsNull()) {
            tolerance = std::max(tolerance, BRep_Tool::Tolerance(edge));
        }
        if (!face.IsNull()) {
            tolerance = std::max(tolerance, BRep_Tool::Tolerance(face));
        }
    }
    catch (Standard_Failure&) {
        tolerance = Precision::Confusion();
    }

    return HasValidPCurveOnFace(edge, face, tolerance, nullptr);
}

bool GeometryUtils::HasValidPCurveOnFace(const TopoDS_Edge& edge,
                                         const TopoDS_Face& face,
                                         double tolerance,
                                         std::string* errorMessage)
{
    if (edge.IsNull() || face.IsNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = "edge or face is null";
        }
        return false;
    }

    const double safeTolerance = std::max(tolerance, Precision::Confusion());
    const double rangeTolerance = std::max(Precision::PConfusion(), Precision::Confusion());

    try {
        TopLoc_Location curveLocation;
        Standard_Real first3d = 0.0;
        Standard_Real last3d = 0.0;
        Handle(Geom_Curve) curve3d = BRep_Tool::Curve(edge, curveLocation, first3d, last3d);
        if (curve3d.IsNull()) {
            if (errorMessage != nullptr) {
                *errorMessage = "edge has no 3d curve";
            }
            return false;
        }

        Standard_Real first2d = 0.0;
        Standard_Real last2d = 0.0;
        Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(edge, face, first2d, last2d);
        if (pcurve.IsNull()) {
            if (errorMessage != nullptr) {
                *errorMessage = "edge has no pcurve on face";
            }
            return false;
        }

        Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
        if (surface.IsNull()) {
            if (errorMessage != nullptr) {
                *errorMessage = "face has no support surface";
            }
            return false;
        }

        if (!std::isfinite(first3d) || !std::isfinite(last3d) ||
            std::abs(last3d - first3d) <= rangeTolerance) {
            if (errorMessage != nullptr) {
                *errorMessage = "invalid 3d curve range";
            }
            return false;
        }

        if (!std::isfinite(first2d) || !std::isfinite(last2d) ||
            std::abs(last2d - first2d) <= rangeTolerance) {
            if (errorMessage != nullptr) {
                *errorMessage = "invalid pcurve range";
            }
            return false;
        }

        if (!BRep_Tool::SameParameter(edge) || !BRep_Tool::SameRange(edge)) {
            if (errorMessage != nullptr) {
                *errorMessage = "invalid same-parameter state";
            }
            return false;
        }

        BRepLib_CheckCurveOnSurface checker(edge, face);
        checker.Perform();
        if (!checker.IsDone()) {
            if (errorMessage != nullptr) {
                *errorMessage = "pcurve validation failed";
            }
            return false;
        }

        if (checker.MaxDistance() > safeTolerance) {
            if (errorMessage != nullptr) {
                *errorMessage = "pcurve projection mismatch";
            }
            return false;
        }
    }
    catch (Standard_Failure&) {
        if (errorMessage != nullptr) {
            *errorMessage = "OCCT failure while validating pcurve";
        }
        return false;
    }

    return true;
}

bool GeometryUtils::RepairPCurveForEdgeOnFace(TopoDS_Edge& edge,
                                              const TopoDS_Face& face,
                                              double tolerance,
                                              std::string* errorMessage)
{
    if (edge.IsNull() || face.IsNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = "Edge or face is null";
        }
        return false;
    }

    if (!std::isfinite(tolerance)) {
        if (errorMessage != nullptr) {
            *errorMessage = "pcurve repair tolerance is not finite";
        }
        return false;
    }
    const double safeTolerance = std::max(tolerance, Precision::Confusion());

    std::string validationError;
    if (HasValidPCurveOnFace(edge, face, safeTolerance, &validationError)) {
        return true;
    }

    std::string repairError;
    try {
        ShapeFix_Edge edgeFix;
        const Standard_Boolean isSeam = BRep_Tool::IsClosed(edge, face);
        edgeFix.FixAddPCurve(edge, face, isSeam, safeTolerance);
        edgeFix.FixReversed2d(edge, face);
        edgeFix.FixSameParameter(edge, face, safeTolerance);
        edgeFix.FixVertexTolerance(edge, face);
        validationError.clear();
        if (HasValidPCurveOnFace(edge, face, safeTolerance, &validationError)) {
            return true;
        }
        repairError = "pcurve repair failed validation: " + validationError;
    }
    catch (Standard_Failure& failure) {
        const Standard_CString message = failure.GetMessageString();
        repairError = "OCCT failure while repairing pcurve: ";
        repairError += message != nullptr ? message : "OCCT failure";
    }

    std::string rebuildError;
    if (rebuildPCurveByProjection(edge, face, safeTolerance, rebuildError)) {
        if (errorMessage != nullptr) {
            errorMessage->clear();
        }
        return true;
    }
    if (errorMessage != nullptr) {
        *errorMessage = repairError + "; " + rebuildError;
    }
    return false;
}

bool GeometryUtils::RepairPCurvesForEdge(TopoDS_Edge& edge,
                                         const TopTools_ListOfShape& faces,
                                         double tolerance,
                                         std::string* errorMessage)
{
    if (edge.IsNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = "edge is null";
        }
        return false;
    }

    bool allFacesValid = true;
    for (TopTools_ListIteratorOfListOfShape it(faces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        std::string localError;
        if (!RepairPCurveForEdgeOnFace(edge, face, tolerance, &localError)) {
            allFacesValid = false;
            if (errorMessage != nullptr && errorMessage->empty()) {
                *errorMessage = localError;
            }
        }
    }

    return allFacesValid;
}

double GeometryUtils::Clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

bool GeometryUtils::IsFinitePoint(const gp_Pnt& point)
{
    return std::isfinite(point.X())
        && std::isfinite(point.Y())
        && std::isfinite(point.Z());
}
