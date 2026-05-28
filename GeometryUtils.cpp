#include "GeometryUtils.h"

#include <BRep_Tool.hxx>
#include <BRepBndLib.hxx>
#include <BRepLib.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <Extrema_ExtPC.hxx>
#include <Geom2d_Curve.hxx>
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

#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

        TopLoc_Location surfaceLocation;
        Handle(Geom_Surface) surface = BRep_Tool::Surface(face, surfaceLocation);
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

        constexpr int sampleCount = 5;
        for (int sample = 0; sample < sampleCount; ++sample) {
            const double t = static_cast<double>(sample) / static_cast<double>(sampleCount - 1);
            const double param3d = first3d + (last3d - first3d) * t;
            const double param2d = first2d + (last2d - first2d) * t;

            gp_Pnt p3d = curve3d->Value(param3d);
            p3d.Transform(curveLocation.Transformation());

            const gp_Pnt2d uv = pcurve->Value(param2d);
            gp_Pnt pSurf = surface->Value(uv.X(), uv.Y());
            pSurf.Transform(surfaceLocation.Transformation());

            const double distance = p3d.Distance(pSurf);
            if (distance > safeTolerance) {
                if (errorMessage != nullptr) {
                    std::ostringstream os;
                    os << "pcurve projection mismatch: distance=" << distance
                       << ", tolerance=" << safeTolerance
                       << ", sample=" << sample;
                    *errorMessage = os.str();
                }
                return false;
            }
        }
    }
    catch (Standard_Failure& failure) {
        if (errorMessage != nullptr) {
            const Standard_CString message = failure.GetMessageString();
            *errorMessage = "OCCT failure while validating pcurve: ";
            *errorMessage += message != nullptr ? message : "OCCT failure";
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
            *errorMessage = "edge or face is null";
        }
        return false;
    }

    const double safeTolerance = std::max(tolerance, Precision::Confusion());

    std::string validationError;
    if (HasValidPCurveOnFace(edge, face, safeTolerance, &validationError)) {
        return true;
    }

    try {
        ShapeFix_Edge edgeFix;
        const Standard_Boolean isSeam = BRep_Tool::IsClosed(edge, face);
        edgeFix.FixAddPCurve(edge, face, isSeam, safeTolerance);
        edgeFix.FixReversed2d(edge, face);
        edgeFix.FixSameParameter(edge, face, safeTolerance);
        edgeFix.FixVertexTolerance(edge, face);

        // SameParameter acts on the edge and may affect this edge's pcurves on other faces.
        BRepLib::SameParameter(edge, safeTolerance);
    }
    catch (Standard_Failure& failure) {
        if (errorMessage != nullptr) {
            const Standard_CString message = failure.GetMessageString();
            *errorMessage = "OCCT failure while repairing pcurve: ";
            *errorMessage += message != nullptr ? message : "OCCT failure";
        }
        return false;
    }

    validationError.clear();
    if (!HasValidPCurveOnFace(edge, face, safeTolerance, &validationError)) {
        if (errorMessage != nullptr) {
            *errorMessage = "pcurve repair failed validation: " + validationError;
        }
        return false;
    }

    return true;
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
