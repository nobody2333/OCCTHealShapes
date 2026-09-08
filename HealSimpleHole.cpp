#include "HealSimpleHole.h"

#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <Bnd_Box.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFill_CurveConstraint.hxx>
#include <BRepFill_Filling.hxx>
#include <BRepGProp.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Extrema_ExtPC.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <GeomAbs_Shape.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomPlate_BuildPlateSurface.hxx>
#include <GeomPlate_MakeApprox.hxx>
#include <GeomPlate_PlateG0Criterion.hxx>
#include <GProp_GProps.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Pln.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shell.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TColgp_SequenceOfXY.hxx>
#include <TColgp_SequenceOfXYZ.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColGeom2d_HArray1OfCurve.hxx>

namespace
{
    struct FillingComplexityLevel
    {
        int maxDegree = 3;
        int maxSegments = 1;
        bool extendedFallback = false;
    };

    struct SurfaceApproximationMetrics
    {
        double approximationError = 0.0;
        double criterionError = 0.0;
        std::string method;
    };

    struct BoundaryEdgeUse
    {
        TopoDS_Edge edge;
        TopoDS_Face ownerFace;
        TopoDS_Vertex firstVertex;
        TopoDS_Vertex lastVertex;
    };

    struct BoundaryLoop
    {
        TopoDS_Wire wire;
        std::vector<TopoDS_Edge> edges;
        std::vector<TopoDS_Face> ownerFaces;
        std::vector<TopoDS_Vertex> vertices;
    };

    std::vector<int> steppedValues(int minimum, int maximum, int step)
    {
        std::vector<int> values;
        for (int value = minimum; value < maximum;) {
            values.push_back(value);
            if (value > maximum - step) {
                break;
            }
            value += step;
        }
        if (values.empty() || values.back() != maximum) {
            values.push_back(maximum);
        }
        return values;
    }

    std::vector<FillingComplexityLevel> buildComplexityLevels(
        int minDegree,
        int minSegments,
        int degreeStep,
        int segmentStep,
        int normalMaxDegree,
        int normalMaxSegments,
        bool includeExtendedFallback,
        int fallbackMaxDegree,
        int fallbackMaxSegments,
        int lowestDegreeMaxSegments)
    {
        const int finalMaxDegree = includeExtendedFallback
            ? fallbackMaxDegree
            : normalMaxDegree;
        const int finalMaxSegments = includeExtendedFallback
            ? fallbackMaxSegments
            : normalMaxSegments;
        const std::vector<int> degrees = steppedValues(
            minDegree,
            finalMaxDegree,
            degreeStep);

        std::vector<FillingComplexityLevel> levels;
        for (const int degree : degrees) {
            const int degreeMaxSegments = degree == minDegree
                ? std::max(finalMaxSegments, lowestDegreeMaxSegments)
                : finalMaxSegments;
            const std::vector<int> segments = steppedValues(
                minSegments,
                degreeMaxSegments,
                segmentStep);
            for (const int segmentCount : segments) {
                FillingComplexityLevel level;
                level.maxDegree = degree;
                level.maxSegments = segmentCount;
                level.extendedFallback = degree > normalMaxDegree
                    || segmentCount > normalMaxSegments;
                levels.push_back(level);
            }
        }
        return levels;
    }

    std::string failureText(const Standard_Failure& failure)
    {
        const char* message = failure.GetMessageString();
        return message != nullptr ? message : "OCCT failure";
    }

    bool checkSurfaceComplexity(
        const Handle(Geom_BSplineSurface)& surface,
        int requestedMaxDegree,
        int maxPolesPerDirection,
        int maxControlPoints,
        std::string& errorMessage)
    {
        if (surface.IsNull()) {
            errorMessage = "B-spline approximation returned a null surface";
            return false;
        }

        const int uDegree = surface->UDegree();
        const int vDegree = surface->VDegree();
        const int uPoles = surface->NbUPoles();
        const int vPoles = surface->NbVPoles();
        const long long controlPointCount = static_cast<long long>(uPoles)
            * static_cast<long long>(vPoles);
        if (uDegree > requestedMaxDegree || vDegree > requestedMaxDegree) {
            std::ostringstream message;
            message << "B-spline degree exceeds the requested level: uDegree="
                    << uDegree
                    << ", vDegree=" << vDegree
                    << ", requestedMaximum=" << requestedMaxDegree;
            errorMessage = message.str();
            return false;
        }
        if (uPoles > maxPolesPerDirection || vPoles > maxPolesPerDirection) {
            std::ostringstream message;
            message << "B-spline pole count exceeds the per-direction limit: uPoles="
                    << uPoles
                    << ", vPoles=" << vPoles
                    << ", maximum=" << maxPolesPerDirection;
            errorMessage = message.str();
            return false;
        }
        if (controlPointCount > maxControlPoints) {
            std::ostringstream message;
            message << "B-spline control-point count exceeds the limit: actual="
                    << controlPointCount
                    << ", maximum=" << maxControlPoints;
            errorMessage = message.str();
            return false;
        }
        return true;
    }

    struct StrictC0PlateData
    {
        Handle(GeomPlate_Surface) plateSurface;
        std::vector<Handle(Geom2d_Curve)> boundaryPCurves;
        TColgp_SequenceOfXY contour2d;
        TColgp_SequenceOfXYZ contour3d;
        double g0Error = 0.0;
        double boundaryTolerance = 0.0;
    };

    bool prepareStrictC0Plate(
        const std::vector<TopoDS_Edge>& boundaryEdges,
        int fillingDegree,
        int fillingPointsOnCurve,
        int fillingIterations,
        double fillingTolerance,
        double fillingAngularTolerance,
        double fillingCurvatureTolerance,
        double maxFillingG0Error,
        StrictC0PlateData& plateData,
        std::string& errorMessage)
    {
        plateData = StrictC0PlateData();
        if (boundaryEdges.empty()) {
            errorMessage = "C0 Plate boundary is empty";
            return false;
        }

        try {
            GeomPlate_BuildPlateSurface plateBuilder(
                std::max(3, fillingDegree),
                fillingPointsOnCurve,
                fillingIterations,
                fillingTolerance,
                fillingTolerance,
                fillingAngularTolerance,
                fillingCurvatureTolerance,
                Standard_False);

            for (const TopoDS_Edge& boundaryEdge : boundaryEdges) {
                Handle(BRepAdaptor_Curve) curveAdaptor = new BRepAdaptor_Curve();
                curveAdaptor->Initialize(boundaryEdge);
                const Handle(Adaptor3d_Curve)& curve3d = curveAdaptor;
                Handle(BRepFill_CurveConstraint) constraint =
                    new BRepFill_CurveConstraint(
                        curve3d,
                        0,
                        fillingPointsOnCurve,
                        fillingTolerance);
                plateBuilder.Add(constraint);
            }
            plateBuilder.SetNbBounds(
                static_cast<Standard_Integer>(boundaryEdges.size()));
            plateBuilder.Perform();
            if (!plateBuilder.IsDone() || plateBuilder.Surface().IsNull()) {
                errorMessage = "strict C0 Plate construction was not done";
                return false;
            }

            plateData.g0Error = plateBuilder.G0Error();
            if (!std::isfinite(plateData.g0Error)
                || plateData.g0Error > maxFillingG0Error) {
                std::ostringstream message;
                message << "strict C0 Plate G0 error is too large: error="
                        << plateData.g0Error
                        << ", maximum=" << maxFillingG0Error;
                errorMessage = message.str();
                return false;
            }

            plateBuilder.Disc2dContour(8, plateData.contour2d);
            plateBuilder.Disc3dContour(8, 0, plateData.contour3d);
            plateData.boundaryTolerance = std::max({
                fillingTolerance,
                maxFillingG0Error,
                plateData.g0Error + fillingTolerance });
            plateData.plateSurface = plateBuilder.Surface();

            const Handle(TColGeom2d_HArray1OfCurve) platePCurves =
                plateBuilder.Curves2d();
            if (platePCurves.IsNull()
                || platePCurves->Length()
                    != static_cast<Standard_Integer>(boundaryEdges.size())) {
                errorMessage = "strict C0 Plate returned an unexpected pcurve count";
                return false;
            }
            for (Standard_Integer index = platePCurves->Lower();
                 index <= platePCurves->Upper();
                 ++index) {
                const Handle(Geom2d_Curve) platePCurve =
                    platePCurves->Value(index);
                if (platePCurve.IsNull()) {
                    errorMessage = "strict C0 Plate returned a null boundary pcurve";
                    return false;
                }
                plateData.boundaryPCurves.push_back(platePCurve);
            }
        }
        catch (const Standard_Failure& failure) {
            errorMessage = "OCCT failure while preparing strict C0 Plate: "
                + failureText(failure);
            return false;
        }
        return true;
    }

    bool approximateStrictC0Plate(
        const StrictC0PlateData& plateData,
        double approximationTolerance,
        double maxApproximationDeviation,
        int maxDegree,
        int maxSegments,
        int maxPolesPerDirection,
        int maxControlPoints,
        TopoDS_Face& temporaryPatchFace,
        std::vector<Handle(Geom2d_Curve)>& boundaryPCurves,
        double& g0Error,
        double& g1Error,
        SurfaceApproximationMetrics& metrics,
        std::string& errorMessage)
    {
        temporaryPatchFace.Nullify();
        boundaryPCurves.clear();
        g0Error = plateData.g0Error;
        g1Error = 0.0;
        metrics = SurfaceApproximationMetrics();
        if (plateData.plateSurface.IsNull()) {
            errorMessage = "prepared strict C0 Plate surface is null";
            return false;
        }

        try {
            GeomPlate_PlateG0Criterion boundaryCriterion(
                plateData.contour2d,
                plateData.contour3d,
                plateData.boundaryTolerance);
            GeomPlate_MakeApprox approximation(
                plateData.plateSurface,
                boundaryCriterion,
                approximationTolerance,
                maxSegments,
                maxDegree,
                GeomAbs_C0,
                1.0);
            metrics.approximationError = approximation.ApproxError();
            metrics.criterionError = approximation.CriterionError();
            metrics.method = "GeomPlate_MakeApprox";
            if (!std::isfinite(metrics.approximationError)
                || metrics.approximationError > maxApproximationDeviation) {
                std::ostringstream message;
                message << "B-spline to Plate approximation deviation is too large: actual="
                        << metrics.approximationError
                        << ", maximum=" << maxApproximationDeviation
                        << ", target=" << approximationTolerance;
                errorMessage = message.str();
                return false;
            }
            const Handle(Geom_BSplineSurface) surface = approximation.Surface();
            if (!checkSurfaceComplexity(
                    surface,
                    maxDegree,
                    maxPolesPerDirection,
                    maxControlPoints,
                    errorMessage)) {
                return false;
            }

            boundaryPCurves = plateData.boundaryPCurves;
            BRep_Builder builder;
            builder.MakeFace(
                temporaryPatchFace,
                surface,
                TopLoc_Location(),
                approximationTolerance);
            if (temporaryPatchFace.IsNull()) {
                errorMessage = "strict C0 Plate returned a null temporary face";
                return false;
            }
        }
        catch (const Standard_Failure& failure) {
            errorMessage = "OCCT failure while approximating strict C0 Plate: "
                + failureText(failure);
            return false;
        }
        return true;
    }

    int listExtent(const TopTools_ListOfShape& shapes)
    {
        int count = 0;
        for (TopTools_ListIteratorOfListOfShape it(shapes); it.More(); it.Next()) {
            ++count;
        }
        return count;
    }

    bool containsSameFace(
        const std::vector<TopoDS_Face>& faces,
        const TopoDS_Face& face)
    {
        for (const TopoDS_Face& existingFace : faces) {
            if (existingFace.IsSame(face)) {
                return true;
            }
        }
        return false;
    }

    bool containsSameEdge(
        const std::vector<TopoDS_Edge>& edges,
        const TopoDS_Edge& edge)
    {
        for (const TopoDS_Edge& existingEdge : edges) {
            if (existingEdge.IsSame(edge)) {
                return true;
            }
        }
        return false;
    }

    TopoDS_Face firstFaceInList(const TopTools_ListOfShape& faces)
    {
        for (TopTools_ListIteratorOfListOfShape it(faces); it.More(); it.Next()) {
            const TopoDS_Face face = TopoDS::Face(it.Value());
            if (!face.IsNull()) {
                return face;
            }
        }
        return TopoDS_Face();
    }

    bool getEdgeVertices(
        const TopoDS_Edge& edge,
        TopoDS_Vertex& firstVertex,
        TopoDS_Vertex& lastVertex)
    {
        firstVertex.Nullify();
        lastVertex.Nullify();
        if (edge.IsNull()) {
            return false;
        }

        TopExp::Vertices(edge, firstVertex, lastVertex, Standard_True);
        return !firstVertex.IsNull() && !lastVertex.IsNull();
    }

    TopoDS_Edge orientEdgeFromTo(
        const TopoDS_Edge& edge,
        const TopoDS_Vertex& fromVertex,
        const TopoDS_Vertex& toVertex)
    {
        TopoDS_Edge orientedEdge = edge;
        TopoDS_Vertex firstVertex;
        TopoDS_Vertex lastVertex;
        if (!getEdgeVertices(orientedEdge, firstVertex, lastVertex)) {
            return orientedEdge;
        }

        if (firstVertex.IsSame(toVertex) && lastVertex.IsSame(fromVertex)) {
            orientedEdge.Reverse();
        }
        return orientedEdge;
    }

    bool edgeHasVertex(
        const BoundaryEdgeUse& edgeUse,
        const TopoDS_Vertex& vertex)
    {
        return !vertex.IsNull()
            && (edgeUse.firstVertex.IsSame(vertex)
                || edgeUse.lastVertex.IsSame(vertex));
    }

    TopoDS_Vertex otherVertex(
        const BoundaryEdgeUse& edgeUse,
        const TopoDS_Vertex& vertex)
    {
        if (edgeUse.firstVertex.IsSame(vertex)) {
            return edgeUse.lastVertex;
        }
        if (edgeUse.lastVertex.IsSame(vertex)) {
            return edgeUse.firstVertex;
        }
        return TopoDS_Vertex();
    }

    bool isUsableBoundaryEdge(
        const TopoDS_Edge& edge,
        const TopoDS_Face& ownerFace)
    {
        if (edge.IsNull() || ownerFace.IsNull()) {
            return false;
        }

        try {
            if (BRep_Tool::Degenerated(edge)
                || BRep_Tool::IsClosed(edge)
                || BRep_Tool::IsClosed(edge, ownerFace)) {
                return false;
            }
        }
        catch (const Standard_Failure&) {
            return false;
        }

        return GeometryUtils::ComputeEdgeLength(edge) > Precision::Confusion();
    }

    std::vector<BoundaryLoop> buildBoundaryLoops(
        const TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap)
    {
        std::vector<BoundaryEdgeUse> boundaryEdges;
        for (int edgeIndex = 1; edgeIndex <= edgeFaceMap.Extent(); ++edgeIndex) {
            const TopTools_ListOfShape& adjacentFaces =
                edgeFaceMap.FindFromIndex(edgeIndex);
            if (listExtent(adjacentFaces) != 1) {
                continue;
            }

            const TopoDS_Edge edge = TopoDS::Edge(edgeFaceMap.FindKey(edgeIndex));
            const TopoDS_Face ownerFace = firstFaceInList(adjacentFaces);
            if (!isUsableBoundaryEdge(edge, ownerFace)) {
                continue;
            }

            BoundaryEdgeUse edgeUse;
            edgeUse.edge = edge;
            edgeUse.ownerFace = ownerFace;
            if (!getEdgeVertices(edge, edgeUse.firstVertex, edgeUse.lastVertex)) {
                continue;
            }
            boundaryEdges.push_back(edgeUse);
        }

        std::vector<BoundaryLoop> loops;
        std::vector<bool> consumed(boundaryEdges.size(), false);
        for (std::size_t startIndex = 0;
             startIndex < boundaryEdges.size();
             ++startIndex) {
            if (consumed[startIndex]) {
                continue;
            }

            BoundaryLoop loop;
            const TopoDS_Vertex loopStart = boundaryEdges[startIndex].firstVertex;
            TopoDS_Vertex currentVertex = boundaryEdges[startIndex].lastVertex;
            loop.edges.push_back(orientEdgeFromTo(
                boundaryEdges[startIndex].edge,
                loopStart,
                currentVertex));
            loop.ownerFaces.push_back(boundaryEdges[startIndex].ownerFace);
            loop.vertices.push_back(loopStart);
            loop.vertices.push_back(currentVertex);
            consumed[startIndex] = true;

            bool closed = false;
            bool failed = false;
            while (!failed && loop.edges.size() <= boundaryEdges.size()) {
                if (currentVertex.IsSame(loopStart)) {
                    closed = true;
                    break;
                }

                int nextIndex = -1;
                for (std::size_t edgeIndex = 0;
                     edgeIndex < boundaryEdges.size();
                     ++edgeIndex) {
                    if (consumed[edgeIndex]
                        || !edgeHasVertex(boundaryEdges[edgeIndex], currentVertex)) {
                        continue;
                    }

                    // 分支边界无法唯一确定孔洞轮廓，整条链保守跳过。
                    if (nextIndex >= 0) {
                        failed = true;
                        break;
                    }
                    nextIndex = static_cast<int>(edgeIndex);
                }

                if (failed || nextIndex < 0) {
                    failed = true;
                    break;
                }

                const BoundaryEdgeUse& nextUse =
                    boundaryEdges[static_cast<std::size_t>(nextIndex)];
                const TopoDS_Vertex nextVertex = otherVertex(nextUse, currentVertex);
                if (nextVertex.IsNull()) {
                    failed = true;
                    break;
                }

                loop.edges.push_back(orientEdgeFromTo(
                    nextUse.edge,
                    currentVertex,
                    nextVertex));
                loop.ownerFaces.push_back(nextUse.ownerFace);
                loop.vertices.push_back(nextVertex);
                consumed[static_cast<std::size_t>(nextIndex)] = true;
                currentVertex = nextVertex;
            }

            if (!closed
                || loop.edges.size() < 3
                || loop.vertices.size() != loop.edges.size() + 1) {
                continue;
            }

            BRepBuilderAPI_MakeWire wireMaker;
            for (const TopoDS_Edge& edge : loop.edges) {
                wireMaker.Add(edge);
                if (!wireMaker.IsDone()) {
                    failed = true;
                    break;
                }
            }
            if (failed || !wireMaker.IsDone()) {
                continue;
            }

            loop.wire = wireMaker.Wire();
            if (loop.wire.IsNull() || !BRep_Tool::IsClosed(loop.wire)) {
                continue;
            }
            loops.push_back(loop);
        }
        return loops;
    }

    bool nearestDistanceToCurve(
        const gp_Pnt& point,
        const BRepAdaptor_Curve& targetCurve,
        double& nearestDistance)
    {
        nearestDistance = std::numeric_limits<double>::max();
        try {
            Extrema_ExtPC extrema(
                point,
                targetCurve,
                targetCurve.FirstParameter(),
                targetCurve.LastParameter(),
                Precision::PConfusion());
            if (extrema.IsDone()) {
                for (int solution = 1; solution <= extrema.NbExt(); ++solution) {
                    nearestDistance = std::min(
                        nearestDistance,
                        std::sqrt(std::max(0.0, extrema.SquareDistance(solution))));
                }
            }

            nearestDistance = std::min(
                nearestDistance,
                point.Distance(targetCurve.Value(targetCurve.FirstParameter())));
            nearestDistance = std::min(
                nearestDistance,
                point.Distance(targetCurve.Value(targetCurve.LastParameter())));
        }
        catch (const Standard_Failure&) {
            return false;
        }

        return std::isfinite(nearestDistance)
            && nearestDistance < std::numeric_limits<double>::max();
    }

    bool buildContinuousProjectedPCurve(
        const TopoDS_Edge& edge,
        const Handle(Geom_Surface)& surface,
        int sampleCount,
        double projectionTolerance,
        Handle(Geom2d_Curve)& pcurve,
        std::string& errorMessage)
    {
        pcurve.Nullify();
        if (edge.IsNull() || surface.IsNull()) {
            errorMessage = "edge or projection surface is null";
            return false;
        }

        sampleCount = std::max(5, sampleCount);
        const BRepAdaptor_Curve edgeCurve(edge);
        const double first = edgeCurve.FirstParameter();
        const double last = edgeCurve.LastParameter();
        if (!std::isfinite(first)
            || !std::isfinite(last)
            || last - first <= Precision::PConfusion()) {
            errorMessage = "invalid edge range for sampled pcurve projection";
            return false;
        }

        Handle(TColgp_HArray1OfPnt2d) uvPoints =
            new TColgp_HArray1OfPnt2d(1, sampleCount);
        Handle(TColStd_HArray1OfReal) parameters =
            new TColStd_HArray1OfReal(1, sampleCount);
        gp_Pnt2d previousUV;
        bool hasPreviousUV = false;
        Standard_Real uMin = 0.0;
        Standard_Real uMax = 0.0;
        Standard_Real vMin = 0.0;
        Standard_Real vMax = 0.0;
        surface->Bounds(uMin, uMax, vMin, vMax);
        if (!std::isfinite(uMin)
            || !std::isfinite(uMax)
            || !std::isfinite(vMin)
            || !std::isfinite(vMax)
            || uMax - uMin <= Precision::PConfusion()
            || vMax - vMin <= Precision::PConfusion()) {
            errorMessage = "invalid patch surface parameter bounds";
            return false;
        }
        const double localHalfU = 0.08 * (uMax - uMin);
        const double localHalfV = 0.08 * (vMax - vMin);

        for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
            const double ratio = static_cast<double>(sampleIndex)
                / static_cast<double>(sampleCount - 1);
            const double parameter = first + (last - first) * ratio;
            const gp_Pnt point = edgeCurve.Value(parameter);
            gp_Pnt2d selectedUV;
            bool projectionDone = false;
            if (hasPreviousUV) {
                double centerU = previousUV.X();
                double centerV = previousUV.Y();
                if (surface->IsUPeriodic()) {
                    const double period = surface->UPeriod();
                    while (centerU < uMin) {
                        centerU += period;
                    }
                    while (centerU > uMax) {
                        centerU -= period;
                    }
                }
                if (surface->IsVPeriodic()) {
                    const double period = surface->VPeriod();
                    while (centerV < vMin) {
                        centerV += period;
                    }
                    while (centerV > vMax) {
                        centerV -= period;
                    }
                }

                const double localUMin = std::max(uMin, centerU - localHalfU);
                const double localUMax = std::min(uMax, centerU + localHalfU);
                const double localVMin = std::max(vMin, centerV - localHalfV);
                const double localVMax = std::min(vMax, centerV + localHalfV);
                GeomAPI_ProjectPointOnSurf localProjector(
                    point,
                    surface,
                    localUMin,
                    localUMax,
                    localVMin,
                    localVMax,
                    projectionTolerance);
                if (localProjector.IsDone() && localProjector.NbPoints() > 0) {
                    double u = 0.0;
                    double v = 0.0;
                    localProjector.LowerDistanceParameters(u, v);
                    selectedUV.SetCoord(u, v);
                    projectionDone = true;
                }
            }

            if (!projectionDone) {
                GeomAPI_ProjectPointOnSurf globalProjector(
                    point,
                    surface,
                    projectionTolerance);
                if (!globalProjector.IsDone() || globalProjector.NbPoints() < 1) {
                    errorMessage = "sampled point projection onto patch surface failed";
                    return false;
                }
                double u = 0.0;
                double v = 0.0;
                globalProjector.LowerDistanceParameters(u, v);
                selectedUV.SetCoord(u, v);
            }

            if (hasPreviousUV && surface->IsUPeriodic()) {
                const double period = surface->UPeriod();
                while (selectedUV.X() - previousUV.X() > 0.5 * period) {
                    selectedUV.SetX(selectedUV.X() - period);
                }
                while (selectedUV.X() - previousUV.X() < -0.5 * period) {
                    selectedUV.SetX(selectedUV.X() + period);
                }
            }
            if (hasPreviousUV && surface->IsVPeriodic()) {
                const double period = surface->VPeriod();
                while (selectedUV.Y() - previousUV.Y() > 0.5 * period) {
                    selectedUV.SetY(selectedUV.Y() - period);
                }
                while (selectedUV.Y() - previousUV.Y() < -0.5 * period) {
                    selectedUV.SetY(selectedUV.Y() + period);
                }
            }
            if (!std::isfinite(selectedUV.X()) || !std::isfinite(selectedUV.Y())) {
                errorMessage = "continuous UV projection returned a non-finite point";
                return false;
            }

            uvPoints->SetValue(sampleIndex + 1, selectedUV);
            parameters->SetValue(sampleIndex + 1, parameter);
            previousUV = selectedUV;
            hasPreviousUV = true;
        }

        Geom2dAPI_Interpolate interpolator(
            uvPoints,
            parameters,
            Standard_False,
            std::max(Precision::PConfusion(), projectionTolerance));
        interpolator.Perform();
        if (!interpolator.IsDone() || interpolator.Curve().IsNull()) {
            errorMessage = "sampled UV interpolation failed";
            return false;
        }

        pcurve = interpolator.Curve();
        return true;
    }

    bool measurePCurveDeviation(
        const TopoDS_Edge& edge,
        const Handle(Geom_Surface)& surface,
        const Handle(Geom2d_Curve)& pcurve,
        int sampleCount,
        double& maxDeviation)
    {
        maxDeviation = std::numeric_limits<double>::infinity();
        if (edge.IsNull() || surface.IsNull() || pcurve.IsNull()) {
            return false;
        }

        try {
            sampleCount = std::max(5, sampleCount);
            const BRepAdaptor_Curve edgeCurve(edge);
            const double first = edgeCurve.FirstParameter();
            const double last = edgeCurve.LastParameter();
            if (!std::isfinite(first)
                || !std::isfinite(last)
                || last - first <= Precision::PConfusion()) {
                return false;
            }

            double measuredMaximum = 0.0;
            for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
                const double ratio = static_cast<double>(sampleIndex)
                    / static_cast<double>(sampleCount - 1);
                const double parameter = first + (last - first) * ratio;
                const gp_Pnt point3d = edgeCurve.Value(parameter);
                const gp_Pnt2d uv = pcurve->Value(parameter);
                if (!std::isfinite(uv.X()) || !std::isfinite(uv.Y())) {
                    return false;
                }
                const gp_Pnt surfacePoint = surface->Value(uv.X(), uv.Y());
                const double deviation = point3d.Distance(surfacePoint);
                if (!std::isfinite(deviation)) {
                    return false;
                }
                measuredMaximum = std::max(measuredMaximum, deviation);
            }
            maxDeviation = measuredMaximum;
            return true;
        }
        catch (const Standard_Failure&) {
            return false;
        }
    }

    bool findSingleOwnerShell(
        const std::vector<TopoDS_Face>& ownerFaces,
        const TopTools_IndexedDataMapOfShapeListOfShape& faceShellMap,
        TopoDS_Shell& ownerShell)
    {
        ownerShell.Nullify();
        for (const TopoDS_Face& face : ownerFaces) {
            if (!faceShellMap.Contains(face)) {
                return false;
            }

            const TopTools_ListOfShape& shells = faceShellMap.FindFromKey(face);
            if (listExtent(shells) != 1) {
                return false;
            }

            TopoDS_Shell faceShell;
            for (TopTools_ListIteratorOfListOfShape it(shells); it.More(); it.Next()) {
                faceShell = TopoDS::Shell(it.Value());
                break;
            }
            if (faceShell.IsNull()) {
                return false;
            }

            if (ownerShell.IsNull()) {
                ownerShell = faceShell;
            }
            else if (!ownerShell.IsSame(faceShell)) {
                return false;
            }
        }
        return !ownerShell.IsNull();
    }

    bool computeWireDiameter(const TopoDS_Wire& wire, double& diameter)
    {
        diameter = 0.0;
        Bnd_Box box;
        BRepBndLib::AddOptimal(wire, box, Standard_False, Standard_True);
        if (box.IsVoid() || box.IsWhole()) {
            return false;
        }

        Standard_Real xMin = 0.0;
        Standard_Real yMin = 0.0;
        Standard_Real zMin = 0.0;
        Standard_Real xMax = 0.0;
        Standard_Real yMax = 0.0;
        Standard_Real zMax = 0.0;
        box.Get(xMin, yMin, zMin, xMax, yMax, zMax);
        const double dx = xMax - xMin;
        const double dy = yMax - yMin;
        const double dz = zMax - zMin;
        diameter = std::sqrt(dx * dx + dy * dy + dz * dz);
        return std::isfinite(diameter);
    }

    int countWireEdgesByConnectivity(const TopoDS_Wire& wire)
    {
        int count = 0;
        for (BRepTools_WireExplorer explorer(wire); explorer.More(); explorer.Next()) {
            ++count;
        }
        return count;
    }

    bool shapeContainsFace(const TopoDS_Shape& shape, const TopoDS_Face& face)
    {
        for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
            if (explorer.Current().IsSame(face)) {
                return true;
            }
        }
        return false;
    }

    int countUniqueFaces(const TopoDS_Shape& shape)
    {
        TopTools_MapOfShape faces;
        for (TopExp_Explorer explorer(shape, TopAbs_FACE);
             explorer.More();
             explorer.Next()) {
            faces.Add(explorer.Current());
        }
        return faces.Extent();
    }

    bool extractSingleShell(
        const TopoDS_Shape& shape,
        TopoDS_Shell& shell)
    {
        shell.Nullify();
        if (shape.IsNull()) {
            return false;
        }
        if (shape.ShapeType() == TopAbs_SHELL) {
            shell = TopoDS::Shell(shape);
            return true;
        }

        TopTools_MapOfShape shells;
        for (TopExp_Explorer explorer(shape, TopAbs_SHELL);
             explorer.More();
             explorer.Next()) {
            if (shells.Add(explorer.Current())) {
                shell = TopoDS::Shell(explorer.Current());
            }
        }
        if (shells.Extent() != 1 || shell.IsNull()) {
            shell.Nullify();
            return false;
        }

        // 如果 shell 外还存在裸 face，说明 patch 并未真正 sewing 进 shell。
        if (countUniqueFaces(shape) != countUniqueFaces(shell)) {
            shell.Nullify();
            return false;
        }
        return true;
    }

    TopoDS_Face mappedFaceAfterSewing(
        BRepBuilderAPI_Sewing& sewing,
        const TopoDS_Face& sourceFace)
    {
        TopoDS_Shape mappedShape = sewing.ModifiedSubShape(sourceFace);
        if (!mappedShape.IsNull() && mappedShape.ShapeType() == TopAbs_FACE) {
            return TopoDS::Face(mappedShape);
        }

        const TopoDS_Shape& directlyModified = sewing.Modified(sourceFace);
        if (!directlyModified.IsNull()
            && directlyModified.ShapeType() == TopAbs_FACE) {
            return TopoDS::Face(directlyModified);
        }
        return TopoDS_Face();
    }
}

HealSimpleHole::HealSimpleHole(const Params& params)
    : params_(params)
{
    params_.minBoundaryEdgeCount = std::max(4, params_.minBoundaryEdgeCount);
    params_.maxBoundaryEdgeCount = std::max(
        params_.minBoundaryEdgeCount,
        params_.maxBoundaryEdgeCount);
    params_.minOwnerFaceCount = std::max(1, params_.minOwnerFaceCount);
    params_.minBoundaryDiameter = std::max(0.0, params_.minBoundaryDiameter);
    if (params_.maxBoundaryDiameter > 0.0) {
        params_.maxBoundaryDiameter = std::max(
            params_.minBoundaryDiameter,
            params_.maxBoundaryDiameter);
    }
    params_.minBoundaryPerimeter = std::max(0.0, params_.minBoundaryPerimeter);
    if (params_.maxBoundaryPerimeter > 0.0) {
        params_.maxBoundaryPerimeter = std::max(
            params_.minBoundaryPerimeter,
            params_.maxBoundaryPerimeter);
    }
    params_.minBoundaryCornerAngle = std::max(
        0.0,
        std::min(3.14159265358979323846, params_.minBoundaryCornerAngle));
    params_.boundarySampleCountPerEdge = std::max(
        1,
        params_.boundarySampleCountPerEdge);
    params_.localEdgeExclusionRatio = std::max(
        0.0,
        std::min(0.49, params_.localEdgeExclusionRatio));
    params_.boundarySeparationPercentile = std::max(
        0.0,
        std::min(1.0, params_.boundarySeparationPercentile));
    params_.minBoundarySeparation = std::max(
        0.0,
        params_.minBoundarySeparation);
    params_.minBoundarySeparationRatio = std::max(
        0.0,
        params_.minBoundarySeparationRatio);
    params_.minPatchArea = std::max(0.0, params_.minPatchArea);
    if (params_.maxPatchArea > 0.0) {
        params_.maxPatchArea = std::max(params_.minPatchArea, params_.maxPatchArea);
    }
    params_.pcurveTolerance = std::max(
        Precision::Confusion(),
        params_.pcurveTolerance);
    params_.maxPCurveTolerance = std::max(
        params_.pcurveTolerance,
        params_.maxPCurveTolerance);
    params_.planarityTolerance = std::min(
        std::max(Precision::Confusion(), params_.planarityTolerance),
        params_.maxPCurveTolerance);
    params_.fillingTolerance = std::max(
        Precision::Confusion(),
        params_.fillingTolerance);
    params_.maxFillingG0Error = std::max(
        params_.fillingTolerance,
        params_.maxFillingG0Error);
    params_.maxRebuiltBoundaryDeviation = std::max(
        params_.fillingTolerance,
        params_.maxRebuiltBoundaryDeviation);
    params_.surfaceApproximationTolerance = std::min(
        params_.maxRebuiltBoundaryDeviation,
        std::max(params_.fillingTolerance,
                 params_.surfaceApproximationTolerance));
    params_.maxSurfaceApproximationDeviation = std::max(
        params_.surfaceApproximationTolerance,
        params_.maxSurfaceApproximationDeviation);
    params_.maxFillingG1Error = std::max(0.0, params_.maxFillingG1Error);
    params_.fillingAngularTolerance = std::max(
        Precision::Angular(),
        params_.fillingAngularTolerance);
    params_.fillingCurvatureTolerance = std::max(
        Precision::Confusion(),
        params_.fillingCurvatureTolerance);
    params_.fillingDegree = std::max(2, params_.fillingDegree);
    params_.fillingPointsOnCurve = std::max(3, params_.fillingPointsOnCurve);
    params_.fillingIterations = std::max(1, params_.fillingIterations);
    params_.approximationMinDegree = std::max(
        2,
        params_.approximationMinDegree);
    params_.approximationMinSegments = std::max(
        1,
        params_.approximationMinSegments);
    params_.approximationDegreeStep = std::max(
        1,
        params_.approximationDegreeStep);
    params_.approximationSegmentStep = std::max(
        1,
        params_.approximationSegmentStep);
    params_.fillingMaxDegree = std::max(
        std::max(params_.fillingDegree, params_.approximationMinDegree),
        params_.fillingMaxDegree);
    params_.fillingMaxSegments = std::max(
        params_.approximationMinSegments,
        params_.fillingMaxSegments);
    params_.fallbackFillingMaxDegree = std::max(
        params_.fillingMaxDegree,
        params_.fallbackFillingMaxDegree);
    params_.fallbackFillingMaxSegments = std::max(
        params_.fillingMaxSegments,
        params_.fallbackFillingMaxSegments);
    params_.lowestDegreeMaxSegments = std::max(
        params_.fallbackFillingMaxSegments,
        params_.lowestDegreeMaxSegments);
    params_.maxFillingPolesPerDirection = std::max(
        1,
        params_.maxFillingPolesPerDirection);
    params_.maxFillingControlPoints = std::max(
        1,
        params_.maxFillingControlPoints);
    params_.finalBoundaryProjectionSamples = std::max(
        5,
        params_.finalBoundaryProjectionSamples);
    params_.localSewingMinTolerance = std::max(
        Precision::Confusion(),
        params_.localSewingMinTolerance);
    params_.localSewingTolerance = std::max(
        params_.localSewingMinTolerance,
        params_.localSewingTolerance);
    params_.localSewingMaxTolerance = std::max(
        params_.localSewingTolerance,
        params_.localSewingMaxTolerance);
    params_.shellFixPrecision = std::max(
        Precision::Confusion(),
        params_.shellFixPrecision);
    params_.shellFixMaxTolerance = std::max(
        params_.shellFixPrecision,
        params_.shellFixMaxTolerance);
    params_.maxRepairCount = std::max(1, params_.maxRepairCount);
}

void HealSimpleHole::setInputShape(const TopoDS_Shape& shape)
{
    inputShape_ = shape;
    resultShape_ = shape;
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();
}

int HealSimpleHole::perform()
{
    resultShape_ = inputShape_;
    modified_ = false;
    repairedCount_ = 0;
    messages_.clear();

    if (inputShape_.IsNull()) {
        appendMessage("Input shape is null");
        return 0;
    }

    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Input shape", inputShape_);

    {
        std::ostringstream parameters;
        parameters << "Parameters: minBoundaryCornerAngle="
                   << params_.minBoundaryCornerAngle
                   << ", boundarySampleCountPerEdge="
                   << params_.boundarySampleCountPerEdge
                   << ", localEdgeExclusionRatio="
                   << params_.localEdgeExclusionRatio
                   << ", boundarySeparationPercentile="
                   << params_.boundarySeparationPercentile
                   << ", minBoundarySeparation="
                   << params_.minBoundarySeparation
                   << ", minBoundarySeparationRatio="
                   << params_.minBoundarySeparationRatio
                   << ", minBoundaryDiameter="
                   << params_.minBoundaryDiameter
                   << ", maxBoundaryDiameter="
                   << params_.maxBoundaryDiameter
                   << ", tryPlanarPatch="
                   << (params_.tryPlanarPatch ? "true" : "false")
                   << ", tryBSplineG1Patch="
                   << (params_.tryBSplineG1Patch ? "true" : "false")
                   << ", tryBSplineC0Patch="
                   << (params_.tryBSplineC0Patch ? "true" : "false")
                   << ", useLocalSewingFallback="
                   << (params_.useLocalSewingFallback ? "true" : "false")
                   << ", fillingTolerance=" << params_.fillingTolerance
                   << ", surfaceApproximationTolerance="
                   << params_.surfaceApproximationTolerance
                   << ", maxSurfaceApproximationDeviation="
                   << params_.maxSurfaceApproximationDeviation
                   << ", maxFillingG0Error=" << params_.maxFillingG0Error
                   << ", maxRebuiltBoundaryDeviation="
                   << params_.maxRebuiltBoundaryDeviation
                   << ", maxFillingG1Error=" << params_.maxFillingG1Error
                   << ", fillingPointsOnCurve=" << params_.fillingPointsOnCurve
                   << ", fillingIterations=" << params_.fillingIterations
                   << ", reprojectBoundaryOnFinalSurface="
                   << (params_.reprojectBoundaryOnFinalSurface
                           ? "true"
                           : "false")
                   << ", finalBoundaryProjectionSamples="
                   << params_.finalBoundaryProjectionSamples
                   << ", approximationMinDegree="
                   << params_.approximationMinDegree
                   << ", approximationMinSegments="
                   << params_.approximationMinSegments
                   << ", approximationDegreeStep="
                   << params_.approximationDegreeStep
                   << ", approximationSegmentStep="
                   << params_.approximationSegmentStep
                   << ", fillingMaxDegree=" << params_.fillingMaxDegree
                   << ", fillingMaxSegments=" << params_.fillingMaxSegments
                   << ", tryHighComplexityC0Fallback="
                   << (params_.tryHighComplexityC0Fallback ? "true" : "false")
                   << ", fallbackFillingMaxDegree="
                   << params_.fallbackFillingMaxDegree
                   << ", fallbackFillingMaxSegments="
                   << params_.fallbackFillingMaxSegments
                   << ", lowestDegreeMaxSegments="
                   << params_.lowestDegreeMaxSegments
                   << ", maxFillingPolesPerDirection="
                   << params_.maxFillingPolesPerDirection
                   << ", maxFillingControlPoints="
                   << params_.maxFillingControlPoints
                   << ", localSewingTolerance="
                   << params_.localSewingTolerance;
        appendMessage(parameters.str());
    }

    TopoDS_Shape currentShape = inputShape_;
    for (int iteration = 1; iteration <= params_.maxRepairCount; ++iteration) {
        RepairDiagnostics::AppendIterationHeader(messages_, iteration, iteration > 1);
        const std::vector<Candidate> candidates = extractCandidates(currentShape);
        if (candidates.empty()) {
            appendMessage("No simple-hole boundary candidates found");
            break;
        }

        // 候选最终结果统一与本轮修复前的 currentShape 比较，而非已修改 pcurve 的 trial。
        // 关闭错误增量拒绝时，保持原行为，不额外收集 BRepCheck 基线。
        std::string baselineError;
        const RepairDiagnostics::FailureList existingFailures = params_.rejectInvalidResult
            ? RepairDiagnostics::CollectBRepCheckFailures(currentShape, &baselineError)
            : RepairDiagnostics::FailureList();
        if (!baselineError.empty()) {
            appendMessage("Repair stopped: failed to collect existing BRepCheck failures: "
                + baselineError);
            break;
        }

        {
            std::ostringstream message;
            message << "Candidate count=" << candidates.size();
            appendMessage(message.str());
        }

        bool committed = false;
        for (std::size_t candidateIndex = 0;
             candidateIndex < candidates.size();
             ++candidateIndex) {
            const Candidate& candidate = candidates[candidateIndex];
            std::ostringstream candidateInfo;
            candidateInfo << "Candidate " << candidateIndex
                          << ": edgeCount=" << candidate.boundaryEdges.size()
                          << ", ownerFaceCount=" << candidate.ownerFaces.size()
                          << ", perimeter=" << candidate.perimeter
                          << ", diameter=" << candidate.diameter
                          << ", minimumCornerAngle="
                          << candidate.minimumCornerAngle
                          << ", minimumSampledSeparation="
                          << candidate.minimumSampledSeparation
                          << ", meanSampledSeparation="
                          << candidate.meanSampledSeparation
                          << ", boundarySeparation="
                          << candidate.boundarySeparation
                          << ", boundarySeparationRatio="
                          << candidate.boundarySeparationRatio;
            appendMessage(candidateInfo.str());

            TopoDS_Shape repairedShape;
            std::string errorMessage;
            if (!tryRepairCandidate(
                    currentShape,
                    existingFailures,
                    candidate,
                    repairedShape,
                    errorMessage)) {
                appendMessage(
                    RepairDiagnostics::StrategyFailureBlock(
                        "HoleRepairStrategies",
                        errorMessage));
                continue;
            }

            currentShape = repairedShape;
            ++repairedCount_;
            modified_ = true;
            committed = true;

            std::ostringstream successMessage;
            successMessage << "Simple hole repaired: edgeCount="
                           << candidate.boundaryEdges.size()
                           << ", perimeter=" << candidate.perimeter
                           << ", diameter=" << candidate.diameter
                           << ", boundarySeparation="
                           << candidate.boundarySeparation
                           << ", boundarySeparationRatio="
                           << candidate.boundarySeparationRatio;
            appendMessage(successMessage.str());
            break;
        }

        if (!committed) {
            appendMessage("No candidate could be committed in this iteration");
            break;
        }
    }

    resultShape_ = currentShape;
    RepairDiagnostics::AppendShapeCheckSummary(messages_, "Result shape", resultShape_);
    return repairedCount_;
}

int HealSimpleHole::perform(const TopoDS_Shape& shape)
{
    setInputShape(shape);
    return perform();
}

std::vector<HealSimpleHole::Candidate> HealSimpleHole::extractCandidates(
    const TopoDS_Shape& shape)
{
    std::vector<Candidate> candidates;
    if (shape.IsNull()) {
        return candidates;
    }

    try {
        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        TopTools_IndexedDataMapOfShapeListOfShape faceShellMap;
        TopologyUtils::BuildEdgeFaceMap(shape, edgeFaceMap);
        TopologyUtils::BuildFaceShellMap(shape, faceShellMap);

        const std::vector<BoundaryLoop> boundaryLoops = buildBoundaryLoops(edgeFaceMap);
        int rejectedBySize = 0;
        int rejectedByOwnerShell = 0;
        int rejectedBySharpCorner = 0;
        int rejectedBySeparation = 0;
        int rejectedByOtherClassification = 0;

        for (const BoundaryLoop& loop : boundaryLoops) {
            Candidate candidate;
            candidate.boundaryWire = loop.wire;
            candidate.boundaryEdges = loop.edges;
            candidate.boundaryVertices = loop.vertices;
            candidate.boundaryOwnerFaces = loop.ownerFaces;
            if (candidate.boundaryWire.IsNull()
                || !BRep_Tool::IsClosed(candidate.boundaryWire)) {
                continue;
            }

            bool candidateValid = true;
            candidate.perimeter = 0.0;
            for (std::size_t edgeIndex = 0;
                 edgeIndex < candidate.boundaryEdges.size();
                 ++edgeIndex) {
                const TopoDS_Edge& edge = candidate.boundaryEdges[edgeIndex];
                const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
                if (!std::isfinite(edgeLength)
                    || edgeLength <= Precision::Confusion()) {
                    candidateValid = false;
                    break;
                }

                candidate.perimeter += edgeLength;

                const TopoDS_Face& ownerFace = loop.ownerFaces[edgeIndex];
                if (!ownerFace.IsNull()
                    && !containsSameFace(candidate.ownerFaces, ownerFace)) {
                    candidate.ownerFaces.push_back(ownerFace);
                }
            }

            if (!candidateValid) {
                continue;
            }

            const int directEdgeCount = static_cast<int>(candidate.boundaryEdges.size());
            if (directEdgeCount < params_.minBoundaryEdgeCount
                || directEdgeCount > params_.maxBoundaryEdgeCount
                || countWireEdgesByConnectivity(candidate.boundaryWire) != directEdgeCount
                || static_cast<int>(candidate.ownerFaces.size()) < params_.minOwnerFaceCount) {
                ++rejectedBySize;
                continue;
            }

            if (!findSingleOwnerShell(
                    candidate.ownerFaces,
                    faceShellMap,
                    candidate.ownerShell)) {
                ++rejectedByOwnerShell;
                continue;
            }

            if (!computeWireDiameter(candidate.boundaryWire, candidate.diameter)) {
                continue;
            }

            if (candidate.diameter < params_.minBoundaryDiameter
                || (params_.maxBoundaryDiameter > 0.0
                    && candidate.diameter > params_.maxBoundaryDiameter)
                || candidate.perimeter < params_.minBoundaryPerimeter
                || (params_.maxBoundaryPerimeter > 0.0
                    && candidate.perimeter > params_.maxBoundaryPerimeter)) {
                ++rejectedBySize;
                continue;
            }

            std::string classificationError;
            if (!evaluateHoleBoundary(candidate, classificationError)) {
                if (classificationError.find("sharp corner") != std::string::npos) {
                    ++rejectedBySharpCorner;
                }
                else if (classificationError.find("boundary separation") != std::string::npos) {
                    ++rejectedBySeparation;
                }
                else {
                    ++rejectedByOtherClassification;
                }

                if (params_.verbose) {
                    appendMessage("Rejected boundary loop: " + classificationError);
                }
                continue;
            }

            candidate.confidence = candidate.boundarySeparationRatio
                + 0.01 / std::max(candidate.diameter, Precision::Confusion());
            candidates.push_back(candidate);
        }

        std::ostringstream summary;
        summary << "Custom boundary classification: closedLoops="
                << boundaryLoops.size()
                << ", accepted=" << candidates.size()
                << ", rejectedBySize=" << rejectedBySize
                << ", rejectedByOwnerShell=" << rejectedByOwnerShell
                << ", rejectedBySharpCorner=" << rejectedBySharpCorner
                << ", rejectedBySeparation=" << rejectedBySeparation
                << ", rejectedByOtherClassification="
                << rejectedByOtherClassification;
        appendMessage(summary.str());
    }
    catch (const Standard_Failure& failure) {
        appendMessage(
            "OCCT failure while extracting simple-hole candidates: "
            + failureText(failure));
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const Candidate& first, const Candidate& second) {
            return first.confidence > second.confidence;
        });
    return candidates;
}

bool HealSimpleHole::evaluateHoleBoundary(
    Candidate& candidate,
    std::string& errorMessage) const
{
    const int edgeCount = static_cast<int>(candidate.boundaryEdges.size());
    if (edgeCount < 4
        || candidate.boundaryVertices.size()
            != candidate.boundaryEdges.size() + 1) {
        errorMessage = "insufficient boundary topology for hole classification";
        return false;
    }

    candidate.minimumCornerAngle = 3.14159265358979323846;
    for (int edgeIndex = 0; edgeIndex < edgeCount; ++edgeIndex) {
        const int previousEdgeIndex = (edgeIndex + edgeCount - 1) % edgeCount;
        const TopoDS_Vertex& commonVertex =
            candidate.boundaryVertices[static_cast<std::size_t>(edgeIndex)];
        const double angle = GeometryUtils::ComputeIntersectionAngle(
            candidate.boundaryEdges[static_cast<std::size_t>(previousEdgeIndex)],
            candidate.boundaryEdges[static_cast<std::size_t>(edgeIndex)],
            commonVertex);
        if (!std::isfinite(angle) || angle < 0.0) {
            errorMessage = "failed to compute boundary corner angle";
            return false;
        }

        candidate.minimumCornerAngle = std::min(
            candidate.minimumCornerAngle,
            angle);
        if (angle + Precision::Angular() < params_.minBoundaryCornerAngle) {
            std::ostringstream message;
            message << "sharp corner detected: angle=" << angle
                    << ", minimumAllowed=" << params_.minBoundaryCornerAngle
                    << ", vertexIndex=" << edgeIndex;
            errorMessage = message.str();
            return false;
        }
    }

    int exclusionCount = std::max(
        1,
        static_cast<int>(std::ceil(
            static_cast<double>(edgeCount)
            * params_.localEdgeExclusionRatio)));
    exclusionCount = std::min(
        exclusionCount,
        std::max(1, (edgeCount - 2) / 2));

    std::vector<double> sampledSeparations;
    for (int sourceIndex = 0; sourceIndex < edgeCount; ++sourceIndex) {
        const BRepAdaptor_Curve sourceCurve(
            candidate.boundaryEdges[static_cast<std::size_t>(sourceIndex)]);
        for (int sampleIndex = 0;
             sampleIndex < params_.boundarySampleCountPerEdge;
             ++sampleIndex) {
            // 不采首尾端点，避免公共顶点的零距离影响边界对侧距离。
            const double normalized =
                static_cast<double>(sampleIndex + 1)
                / static_cast<double>(params_.boundarySampleCountPerEdge + 1);
            const double parameter = sourceCurve.FirstParameter()
                + (sourceCurve.LastParameter() - sourceCurve.FirstParameter())
                    * normalized;
            const gp_Pnt sourcePoint = sourceCurve.Value(parameter);

            double nearestNonLocalDistance = std::numeric_limits<double>::max();
            bool foundProjectionTarget = false;
            for (int targetIndex = 0; targetIndex < edgeCount; ++targetIndex) {
                const int indexDifference = std::abs(targetIndex - sourceIndex);
                const int cyclicDistance = std::min(
                    indexDifference,
                    edgeCount - indexDifference);
                if (cyclicDistance <= exclusionCount) {
                    continue;
                }

                const BRepAdaptor_Curve targetCurve(
                    candidate.boundaryEdges[static_cast<std::size_t>(targetIndex)]);
                double projectionDistance = 0.0;
                if (!nearestDistanceToCurve(
                        sourcePoint,
                        targetCurve,
                        projectionDistance)) {
                    continue;
                }

                nearestNonLocalDistance = std::min(
                    nearestNonLocalDistance,
                    projectionDistance);
                foundProjectionTarget = true;
            }

            if (!foundProjectionTarget
                || !std::isfinite(nearestNonLocalDistance)) {
                errorMessage = "failed to project samples to non-local boundary edges";
                return false;
            }
            sampledSeparations.push_back(nearestNonLocalDistance);
        }
    }

    if (sampledSeparations.empty()) {
        errorMessage = "boundary separation sampling produced no values";
        return false;
    }

    std::sort(sampledSeparations.begin(), sampledSeparations.end());
    candidate.minimumSampledSeparation = sampledSeparations.front();
    double separationSum = 0.0;
    for (const double separation : sampledSeparations) {
        separationSum += separation;
    }
    candidate.meanSampledSeparation = separationSum
        / static_cast<double>(sampledSeparations.size());

    const std::size_t percentileIndex = static_cast<std::size_t>(std::floor(
        params_.boundarySeparationPercentile
        * static_cast<double>(sampledSeparations.size() - 1)));
    candidate.boundarySeparation = sampledSeparations[percentileIndex];
    candidate.boundarySeparationRatio = candidate.boundarySeparation
        / std::max(candidate.diameter, Precision::Confusion());

    if (candidate.boundarySeparation + Precision::Confusion()
            < params_.minBoundarySeparation
        || candidate.boundarySeparationRatio + Precision::Confusion()
            < params_.minBoundarySeparationRatio) {
        std::ostringstream message;
        message << "boundary separation too small: percentileDistance="
                << candidate.boundarySeparation
                << ", separationRatio=" << candidate.boundarySeparationRatio
                << ", minimumDistance=" << params_.minBoundarySeparation
                << ", minimumRatio=" << params_.minBoundarySeparationRatio;
        errorMessage = message.str();
        return false;
    }

    return true;
}

bool HealSimpleHole::buildPlanarPatchFace(
    const Candidate& candidate,
    TopoDS_Face& patchFace,
    double& patchArea,
    std::string& errorMessage) const
{
    patchFace.Nullify();
    patchArea = 0.0;
    if (candidate.boundaryWire.IsNull()) {
        errorMessage = "boundary wire is null";
        return false;
    }

    try {
        BRepBuilderAPI_FindPlane planeFinder(
            candidate.boundaryWire,
            params_.planarityTolerance);
        if (!planeFinder.Found() || planeFinder.Plane().IsNull()) {
            errorMessage = "boundary is not planar within planarityTolerance";
            return false;
        }

        BRepBuilderAPI_MakeFace faceMaker(
            planeFinder.Plane()->Pln(),
            candidate.boundaryWire,
            Standard_True);
        if (!faceMaker.IsDone()) {
            errorMessage = "BRepBuilderAPI_MakeFace failed for planar boundary";
            return false;
        }

        patchFace = faceMaker.Face();
        if (patchFace.IsNull()) {
            errorMessage = "planar face builder returned a null face";
            return false;
        }

        if (!finalizePatchFace(
                candidate,
                patchFace,
                params_.planarityTolerance,
                patchArea,
                errorMessage)) {
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while constructing planar patch face: "
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealSimpleHole::buildTemporaryBSplinePatchFace(
    const Candidate& candidate,
    bool useG1Continuity,
    BRepFill_Filling& filling,
    TopoDS_Face& temporaryPatchFace,
    std::vector<Handle(Geom2d_Curve)>& boundaryPCurves,
    double& g0Error,
    double& g1Error,
    std::string& errorMessage,
    int maxDegree,
    int maxSegments) const
{
    temporaryPatchFace.Nullify();
    boundaryPCurves.clear();
    g0Error = 0.0;
    g1Error = 0.0;
    if (candidate.boundaryEdges.size() != candidate.boundaryOwnerFaces.size()) {
        errorMessage = "boundary edge and owner-face arrays have different sizes";
        return false;
    }

    if (!useG1Continuity) {
        return buildStrictC0BSplinePatchFace(
            candidate,
            temporaryPatchFace,
            boundaryPCurves,
            g0Error,
            g1Error,
            errorMessage,
            maxDegree,
            maxSegments);
    }

    try {
        filling.SetResolParam(
            params_.fillingDegree,
            params_.fillingPointsOnCurve,
            params_.fillingIterations,
            Standard_False);
        filling.SetConstrParam(
            params_.fillingTolerance,
            params_.fillingTolerance,
            params_.fillingAngularTolerance,
            params_.fillingCurvatureTolerance);
        filling.SetApproxParam(
            maxDegree,
            maxSegments);

        for (std::size_t edgeIndex = 0;
             edgeIndex < candidate.boundaryEdges.size();
             ++edgeIndex) {
            const TopoDS_Edge& edge = candidate.boundaryEdges[edgeIndex];
            const TopoDS_Face& ownerFace =
                candidate.boundaryOwnerFaces[edgeIndex];
            if (ownerFace.IsNull()) {
                errorMessage = "G1 filling boundary has a null owner face";
                return false;
            }
            filling.Add(edge, ownerFace, GeomAbs_G1, Standard_True);
        }

        filling.Build();
        if (!filling.IsDone()) {
            errorMessage = "BRepFill_Filling G1 build was not done";
            return false;
        }

        temporaryPatchFace = filling.Face();
        if (temporaryPatchFace.IsNull()) {
            errorMessage = "BRepFill_Filling returned a null face";
            return false;
        }

        g0Error = filling.G0Error();
        if (!std::isfinite(g0Error) || g0Error > params_.maxFillingG0Error) {
            std::ostringstream message;
            message << "B-spline filling G0 error is too large: error="
                    << g0Error
                    << ", maximum=" << params_.maxFillingG0Error;
            errorMessage = message.str();
            return false;
        }

        g1Error = filling.G1Error();
        if (!std::isfinite(g1Error) || g1Error > params_.maxFillingG1Error) {
            std::ostringstream message;
            message << "B-spline filling G1 error is too large: error="
                    << g1Error
                    << ", maximum=" << params_.maxFillingG1Error;
            errorMessage = message.str();
            return false;
        }

        TopLoc_Location surfaceLocation;
        const Handle(Geom_Surface) surface =
            BRep_Tool::Surface(temporaryPatchFace, surfaceLocation);
        if (surface.IsNull()) {
            errorMessage = "B-spline filling face has no support surface";
            return false;
        }
        const Handle(Geom_BSplineSurface) bsplineSurface =
            Handle(Geom_BSplineSurface)::DownCast(surface);
        if (bsplineSurface.IsNull()) {
            errorMessage = "BRepFill_Filling result surface is not a Geom_BSplineSurface";
            return false;
        }
        if (!checkSurfaceComplexity(
                bsplineSurface,
                maxDegree,
                params_.maxFillingPolesPerDirection,
                params_.maxFillingControlPoints,
                errorMessage)) {
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = std::string("OCCT failure while building G1 B-spline patch: ")
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealSimpleHole::buildStrictC0BSplinePatchFace(
    const Candidate& candidate,
    TopoDS_Face& temporaryPatchFace,
    std::vector<Handle(Geom2d_Curve)>& boundaryPCurves,
    double& g0Error,
    double& g1Error,
    std::string& errorMessage,
    int maxDegree,
    int maxSegments) const
{
    StrictC0PlateData plateData;
    SurfaceApproximationMetrics metrics;
    if (!prepareStrictC0Plate(
            candidate.boundaryEdges,
            params_.fillingDegree,
            params_.fillingPointsOnCurve,
            params_.fillingIterations,
            params_.fillingTolerance,
            params_.fillingAngularTolerance,
            params_.fillingCurvatureTolerance,
            params_.maxFillingG0Error,
            plateData,
            errorMessage)) {
        return false;
    }
    return approximateStrictC0Plate(
        plateData,
        params_.surfaceApproximationTolerance,
        params_.maxSurfaceApproximationDeviation,
        maxDegree,
        maxSegments,
        params_.maxFillingPolesPerDirection,
        params_.maxFillingControlPoints,
        temporaryPatchFace,
        boundaryPCurves,
        g0Error,
        g1Error,
        metrics,
        errorMessage);
}

bool HealSimpleHole::rebuildBSplinePatchWithOriginalBoundary(
    const Candidate& candidate,
    const TopoDS_Face& temporaryPatchFace,
    const std::vector<Handle(Geom2d_Curve)>& boundaryPCurves,
    double g0Error,
    TopoDS_Face& patchFace,
    double& patchArea,
    int& reprojectedEdgeCount,
    int& reusedPlatePCurveCount,
    double& maxSelectedPCurveDeviation,
    std::string& errorMessage) const
{
    patchFace.Nullify();
    patchArea = 0.0;
    reprojectedEdgeCount = 0;
    reusedPlatePCurveCount = 0;
    maxSelectedPCurveDeviation = 0.0;
    if (temporaryPatchFace.IsNull()) {
        errorMessage = "temporary B-spline patch face is null";
        return false;
    }

    try {
        const Handle(Geom_Surface) locatedFillingSurface =
            BRep_Tool::Surface(temporaryPatchFace);
        if (locatedFillingSurface.IsNull()) {
            errorMessage = "temporary B-spline patch has no support surface";
            return false;
        }

        double maxEdgeTolerance = Precision::Confusion();
        for (const TopoDS_Edge& edge : candidate.boundaryEdges) {
            maxEdgeTolerance = std::max(
                maxEdgeTolerance,
                BRep_Tool::Tolerance(edge));
        }
        const double requiredTolerance = std::max({
            params_.pcurveTolerance,
            params_.fillingTolerance,
            g0Error,
            maxEdgeTolerance });
        // 复用输入边时不能把该边已有的拓扑容差当作本次修复新增的容差。
        // 配置上限只限制算法额外放大的部分，已有边容差始终作为最低可接受上限。
        const double allowedTolerance = std::max(
            params_.maxFillingG0Error,
            std::max(params_.maxPCurveTolerance, maxEdgeTolerance));
        if (requiredTolerance > allowedTolerance) {
            std::ostringstream message;
            message << "B-spline pcurve projection requires excessive tolerance: required="
                    << requiredTolerance
                    << ", maximum=" << allowedTolerance
                    << ", configuredMaximum=" << params_.maxPCurveTolerance
                    << ", existingEdgeTolerance=" << maxEdgeTolerance;
            errorMessage = message.str();
            return false;
        }

        BRep_Builder builder;
        for (std::size_t edgeIndex = 0;
             edgeIndex < candidate.boundaryEdges.size();
             ++edgeIndex) {
            const TopoDS_Edge& originalEdgeUse =
                candidate.boundaryEdges[edgeIndex];
            Standard_Real first3d = 0.0;
            Standard_Real last3d = 0.0;
            Handle(Geom_Curve) curve3d =
                BRep_Tool::Curve(originalEdgeUse, first3d, last3d);
            if (curve3d.IsNull()
                || !std::isfinite(first3d)
                || !std::isfinite(last3d)
                || last3d - first3d <= Precision::PConfusion()) {
                errorMessage = "boundary edge has no valid 3d curve range";
                return false;
            }

            const double projectionPrecision = std::max(
                Precision::Confusion(),
                params_.pcurveTolerance);
            Handle(Geom2d_Curve) platePCurve;
            if (boundaryPCurves.size() == candidate.boundaryEdges.size()) {
                platePCurve = boundaryPCurves[edgeIndex];
            }

            Handle(Geom2d_Curve) finalSurfacePCurve;
            std::string projectionError;
            if (params_.reprojectBoundaryOnFinalSurface) {
                buildContinuousProjectedPCurve(
                    originalEdgeUse,
                    locatedFillingSurface,
                    params_.finalBoundaryProjectionSamples,
                    projectionPrecision,
                    finalSurfacePCurve,
                    projectionError);
            }

            double plateDeviation = std::numeric_limits<double>::infinity();
            const bool plateDeviationMeasured = measurePCurveDeviation(
                originalEdgeUse,
                locatedFillingSurface,
                platePCurve,
                params_.finalBoundaryProjectionSamples,
                plateDeviation);
            double projectedDeviation = std::numeric_limits<double>::infinity();
            const bool projectedDeviationMeasured = measurePCurveDeviation(
                originalEdgeUse,
                locatedFillingSurface,
                finalSurfacePCurve,
                params_.finalBoundaryProjectionSamples,
                projectedDeviation);

            Handle(Geom2d_Curve) selectedPCurve;
            double selectedDeviation = std::numeric_limits<double>::infinity();
            if (projectedDeviationMeasured
                && (!plateDeviationMeasured
                    || projectedDeviation <= plateDeviation)) {
                selectedPCurve = finalSurfacePCurve;
                selectedDeviation = projectedDeviation;
                ++reprojectedEdgeCount;
            }
            else if (plateDeviationMeasured) {
                selectedPCurve = platePCurve;
                selectedDeviation = plateDeviation;
                ++reusedPlatePCurveCount;
            }

            if (selectedPCurve.IsNull()) {
                if (!params_.reprojectBoundaryOnFinalSurface) {
                    projectionError = "final-surface reprojection is disabled";
                }
                else if (projectionError.empty()) {
                    projectionError = "projected and Plate pcurves could not be measured";
                }
                errorMessage = "no measurable boundary pcurve is available: "
                    + projectionError;
                return false;
            }
            maxSelectedPCurveDeviation = std::max(
                maxSelectedPCurveDeviation,
                selectedDeviation);

            builder.UpdateEdge(
                originalEdgeUse,
                selectedPCurve,
                locatedFillingSurface,
                TopLoc_Location(),
                requiredTolerance);
            builder.Range(
                originalEdgeUse,
                locatedFillingSurface,
                TopLoc_Location(),
                first3d,
                last3d);
        }

        // MakeFace 会根据曲面参数域把闭环定向为有限外边界；底层 Builder::Add 不会做该检查。
        BRepBuilderAPI_MakeFace faceMaker(
            locatedFillingSurface,
            candidate.boundaryWire,
            Standard_True);
        if (!faceMaker.IsDone()) {
            errorMessage = "BRepBuilderAPI_MakeFace failed for B-spline boundary";
            return false;
        }
        patchFace = faceMaker.Face();
        if (patchFace.IsNull()) {
            errorMessage = "BRepBuilderAPI_MakeFace returned a null B-spline face";
            return false;
        }
        builder.UpdateFace(patchFace, requiredTolerance);

        // filling 报告的 G0Error 来自内部离散检查；重建为原边后，再用最终 edge/face
        // 组合测量一次真实偏差。只允许在“已有边容差 + filling 工作容差”内补偿。
        const double allowedMeasuredTolerance = std::max({
            params_.maxPCurveTolerance,
            params_.maxRebuiltBoundaryDeviation,
            maxEdgeTolerance + params_.fillingTolerance });
        double measuredTolerance = requiredTolerance;
        for (const TopoDS_Edge& boundaryEdge : candidate.boundaryEdges) {
            BRepLib_CheckCurveOnSurface checker(boundaryEdge, patchFace);
            checker.Perform();
            if (!checker.IsDone() || !std::isfinite(checker.MaxDistance())) {
                errorMessage = "failed to measure rebuilt B-spline boundary deviation";
                return false;
            }

            const double edgeRequiredTolerance = std::max(
                BRep_Tool::Tolerance(boundaryEdge),
                checker.MaxDistance() + Precision::Confusion());
            if (edgeRequiredTolerance > allowedMeasuredTolerance) {
                std::ostringstream message;
                message << "rebuilt B-spline boundary deviation exceeds tolerance budget: required="
                        << edgeRequiredTolerance
                        << ", maximum=" << allowedMeasuredTolerance;
                errorMessage = message.str();
                return false;
            }

            builder.UpdateEdge(boundaryEdge, edgeRequiredTolerance);
            measuredTolerance = std::max(
                measuredTolerance,
                edgeRequiredTolerance);
        }
        builder.UpdateFace(patchFace, measuredTolerance);

        if (!finalizePatchFace(
                candidate,
                patchFace,
                measuredTolerance,
                patchArea,
                errorMessage)) {
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while transferring B-spline pcurves: "
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealSimpleHole::finalizePatchFace(
    const Candidate& candidate,
    TopoDS_Face& patchFace,
    double baseTolerance,
    double& patchArea,
    std::string& errorMessage) const
{
    patchArea = 0.0;
    if (patchFace.IsNull()) {
        errorMessage = "patch face is null";
        return false;
    }

    try {
        double maxEdgeTolerance = Precision::Confusion();
        for (const TopoDS_Edge& edge : candidate.boundaryEdges) {
            maxEdgeTolerance = std::max(
                maxEdgeTolerance,
                BRep_Tool::Tolerance(edge));
        }

        const double requiredTolerance = std::max({
            params_.pcurveTolerance,
            baseTolerance,
            maxEdgeTolerance });
        // 已有边容差不是补面过程新引入的误差，不应被固定配置上限反向拒绝。
        const double allowedTolerance = std::max(
            params_.maxPCurveTolerance,
            maxEdgeTolerance);
        if (requiredTolerance > allowedTolerance) {
            std::ostringstream message;
            message << "required pcurve tolerance exceeds maxPCurveTolerance: required="
                    << requiredTolerance
                    << ", maximum=" << allowedTolerance
                    << ", configuredMaximum=" << params_.maxPCurveTolerance
                    << ", existingEdgeTolerance=" << maxEdgeTolerance;
            errorMessage = message.str();
            return false;
        }
        const double pcurveTolerance = std::max(
            requiredTolerance,
            Precision::Confusion());

        BRep_Builder builder;
        builder.UpdateFace(patchFace, pcurveTolerance);

        std::vector<TopoDS_Edge> patchEdges;
        for (TopExp_Explorer edgeExplorer(patchFace, TopAbs_EDGE);
             edgeExplorer.More();
             edgeExplorer.Next()) {
            TopoDS_Edge edgeUse = TopoDS::Edge(edgeExplorer.Current());
            if (containsSameEdge(patchEdges, edgeUse)) {
                continue;
            }
            patchEdges.push_back(edgeUse);

            std::string pcurveError;
            if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                    edgeUse,
                    patchFace,
                    pcurveTolerance,
                    &pcurveError)) {
                errorMessage = "failed to prepare boundary pcurve on patch face: "
                    + pcurveError;
                return false;
            }
        }

        if (patchEdges.size() != candidate.boundaryEdges.size()) {
            errorMessage = "patch face boundary edge count differs from candidate wire";
            return false;
        }
        for (const TopoDS_Edge& candidateEdge : candidate.boundaryEdges) {
            if (!containsSameEdge(patchEdges, candidateEdge)) {
                errorMessage = "patch face did not reuse every candidate boundary edge";
                return false;
            }
        }

        BRepCheck_Analyzer faceAnalyzer(patchFace, Standard_True);
        if (!faceAnalyzer.IsValid()) {
            errorMessage = "constructed patch face is invalid: "
                + RepairDiagnostics::DescribeBRepCheckFailures(
                    faceAnalyzer,
                    patchFace);
            return false;
        }

        GProp_GProps properties;
        BRepGProp::SurfaceProperties(patchFace, properties);
        patchArea = properties.Mass();
        if (!std::isfinite(patchArea)
            || patchArea <= Precision::SquareConfusion()) {
            errorMessage = "constructed patch face has invalid area";
            return false;
        }
        if (patchArea < params_.minPatchArea
            || (params_.maxPatchArea > 0.0
                && patchArea > params_.maxPatchArea)) {
            std::ostringstream message;
            message << "patch area outside configured range: area=" << patchArea;
            errorMessage = message.str();
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while finalizing patch face: "
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealSimpleHole::addPatchFaceToOwnerShell(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    const TopoDS_Face& patchFace,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    repairedShape.Nullify();
    if (currentShape.IsNull()
        || candidate.ownerShell.IsNull()
        || patchFace.IsNull()) {
        errorMessage = "current shape, owner shell, or patch face is null";
        return false;
    }

    try {
        BRep_Builder builder;
        TopoDS_Shell rebuiltShell;
        builder.MakeShell(rebuiltShell);

        int existingFaceCount = 0;
        for (TopoDS_Iterator childIterator(
                 candidate.ownerShell,
                 Standard_False,
                 Standard_True);
             childIterator.More();
             childIterator.Next()) {
            if (childIterator.Value().ShapeType() != TopAbs_FACE) {
                errorMessage = "owner shell contains a direct non-face child";
                return false;
            }
            builder.Add(rebuiltShell, childIterator.Value());
            ++existingFaceCount;
        }
        if (existingFaceCount == 0) {
            errorMessage = "owner shell has no direct faces";
            return false;
        }

        builder.Add(rebuiltShell, patchFace);

        TopoDS_Shell acceptedShell = rebuiltShell;
        if (params_.fixShellOrientation) {
            ShapeFix_Shell shellFixer;
            shellFixer.SetPrecision(params_.shellFixPrecision);
            shellFixer.SetMinTolerance(Precision::Confusion());
            shellFixer.SetMaxTolerance(params_.shellFixMaxTolerance);
            shellFixer.FixFaceOrientation(
                rebuiltShell,
                Standard_False,
                Standard_False);

            const TopoDS_Shape fixedShape = shellFixer.Shape();
            if (fixedShape.IsNull() || fixedShape.ShapeType() != TopAbs_SHELL) {
                errorMessage = "ShapeFix_Shell did not return one shell";
                return false;
            }
            acceptedShell = TopoDS::Shell(fixedShape);
        }

        acceptedShell.Closed(BRep_Tool::IsClosed(acceptedShell));

        TopoDS_Shell sourceKey = candidate.ownerShell;
        sourceKey.Orientation(TopAbs_FORWARD);
        TopoDS_Shell replacementValue = acceptedShell;
        replacementValue.Orientation(TopAbs_FORWARD);

        BRepTools_ReShape reshaper;
        reshaper.Replace(sourceKey, replacementValue);
        repairedShape = reshaper.Apply(currentShape);
        if (repairedShape.IsNull()) {
            errorMessage = "owner shell replacement produced a null shape";
            return false;
        }
        if (repairedShape.ShapeType() != currentShape.ShapeType()) {
            errorMessage = "owner shell replacement changed top-level shape type";
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while adding patch face to owner shell: "
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealSimpleHole::tryCommitPatchFace(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const Candidate& candidate,
    const TopoDS_Face& patchFace,
    const std::vector<TopoDS_Edge>& resultBoundaryEdges,
    int freeEdgesBefore,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage) const
{
    repairedShape.Nullify();
    if (!addPatchFaceToOwnerShell(
            currentShape,
            candidate,
            patchFace,
            repairedShape,
            errorMessage)) {
        return false;
    }

    if (!validateRepairResult(
            existingFailures,
            repairedShape,
            candidate,
            patchFace,
            resultBoundaryEdges,
            freeEdgesBefore,
            errorMessage)) {
        repairedShape.Nullify();
        return false;
    }
    return true;
}

bool HealSimpleHole::tryRepairByLocalSewing(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);

    try {
        BRepFill_Filling filling;
        TopoDS_Face temporaryPatchFace;
        std::vector<Handle(Geom2d_Curve)> boundaryPCurves;
        double g0Error = 0.0;
        double g1Error = 0.0;
        if (!buildTemporaryBSplinePatchFace(
                candidate,
                Standard_False,
                filling,
                temporaryPatchFace,
                boundaryPCurves,
                g0Error,
                g1Error,
                errorMessage,
                params_.tryHighComplexityC0Fallback
                    ? params_.fallbackFillingMaxDegree
                    : params_.fillingMaxDegree,
                params_.tryHighComplexityC0Fallback
                    ? params_.fallbackFillingMaxSegments
                    : params_.fillingMaxSegments)) {
            return false;
        }

        BRepBuilderAPI_Sewing sewing(
            params_.localSewingTolerance,
            Standard_True,
            Standard_True,
            Standard_False,
            Standard_False);
        sewing.SetTolerance(params_.localSewingTolerance);
        sewing.SetMinTolerance(params_.localSewingMinTolerance);
        sewing.SetMaxTolerance(params_.localSewingMaxTolerance);
        sewing.SetFaceMode(Standard_True);
        sewing.SetSameParameterMode(Standard_True);
        sewing.SetFloatingEdgesMode(Standard_False);
        sewing.SetLocalTolerancesMode(Standard_True);
        sewing.SetNonManifoldMode(Standard_False);
        sewing.Add(candidate.ownerShell);
        sewing.Add(temporaryPatchFace);
        sewing.Perform();

        const TopoDS_Shape rawSewedShape = sewing.SewedShape();
        if (rawSewedShape.IsNull()) {
            errorMessage = "local sewing returned a null shape";
            return false;
        }

        TopoDS_Shell sewedShell;
        if (!extractSingleShell(rawSewedShape, sewedShell)) {
            std::ostringstream message;
            message << "local sewing did not produce exactly one shell without loose faces"
                    << ": resultType="
                    << TopologyUtils::ShapeTypeName(rawSewedShape)
                    << ", freeEdges=" << sewing.NbFreeEdges()
                    << ", contiguousEdges=" << sewing.NbContigousEdges()
                    << ", multipleEdges=" << sewing.NbMultipleEdges();
            errorMessage = message.str();
            return false;
        }

        const int expectedFaceCount = countUniqueFaces(candidate.ownerShell) + 1;
        if (countUniqueFaces(sewedShell) != expectedFaceCount) {
            std::ostringstream message;
            message << "local sewing face count mismatch: expected="
                    << expectedFaceCount
                    << ", actual=" << countUniqueFaces(sewedShell);
            errorMessage = message.str();
            return false;
        }

        TopoDS_Face mappedPatchFace = mappedFaceAfterSewing(
            sewing,
            temporaryPatchFace);
        if (mappedPatchFace.IsNull()) {
            errorMessage = "cannot map temporary patch face through local sewing";
            return false;
        }

        std::vector<TopoDS_Edge> mappedBoundaryEdges;
        for (const TopoDS_Edge& originalEdge : candidate.boundaryEdges) {
            TopoDS_Shape mappedShape = sewing.ModifiedSubShape(originalEdge);
            if (mappedShape.IsNull()) {
                mappedShape = originalEdge;
            }
            if (mappedShape.ShapeType() != TopAbs_EDGE) {
                errorMessage = "local sewing split a boundary edge into a non-edge result";
                return false;
            }

            const TopoDS_Edge mappedEdge = TopoDS::Edge(mappedShape);
            if (!containsSameEdge(mappedBoundaryEdges, mappedEdge)) {
                mappedBoundaryEdges.push_back(mappedEdge);
            }
        }
        if (mappedBoundaryEdges.size() != candidate.boundaryEdges.size()) {
            errorMessage = "local sewing merged multiple boundary edges unexpectedly";
            return false;
        }

        TopoDS_Shell acceptedShell = sewedShell;
        if (params_.fixShellOrientation) {
            ShapeFix_Shell shellFixer;
            shellFixer.SetPrecision(params_.shellFixPrecision);
            shellFixer.SetMinTolerance(Precision::Confusion());
            shellFixer.SetMaxTolerance(params_.shellFixMaxTolerance);
            shellFixer.FixFaceOrientation(
                sewedShell,
                Standard_False,
                Standard_False);
            const TopoDS_Shape fixedShape = shellFixer.Shape();
            if (fixedShape.IsNull() || fixedShape.ShapeType() != TopAbs_SHELL) {
                errorMessage = "ShapeFix_Shell failed after local sewing";
                return false;
            }
            acceptedShell = TopoDS::Shell(fixedShape);
        }
        acceptedShell.Closed(BRep_Tool::IsClosed(acceptedShell));

        TopoDS_Shell sourceKey = candidate.ownerShell;
        sourceKey.Orientation(TopAbs_FORWARD);
        TopoDS_Shell replacementValue = acceptedShell;
        replacementValue.Orientation(TopAbs_FORWARD);

        BRepTools_ReShape reshaper;
        reshaper.Replace(sourceKey, replacementValue);
        repairedShape = reshaper.Apply(currentShape);
        if (repairedShape.IsNull()
            || repairedShape.ShapeType() != currentShape.ShapeType()) {
            errorMessage = "local sewing shell replacement produced an invalid root shape";
            repairedShape.Nullify();
            return false;
        }

        if (!validateRepairResult(
                existingFailures,
                repairedShape,
                candidate,
                mappedPatchFace,
                mappedBoundaryEdges,
                freeEdgesBefore,
                errorMessage)) {
            repairedShape.Nullify();
            return false;
        }

        std::ostringstream success;
        success << "Local sewing patch accepted: g0Error=" << g0Error
                << ", freeEdges=" << sewing.NbFreeEdges()
                << ", contiguousEdges=" << sewing.NbContigousEdges();
        appendMessage(success.str());
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure during local patch sewing: "
            + failureText(failure);
        return false;
    }

    return true;
}

bool HealSimpleHole::createTrialCandidate(
    const TopoDS_Shape& currentShape,
    const Candidate& candidate,
    TopoDS_Shape& trialShape,
    Candidate& trialCandidate,
    std::string& errorMessage) const
{
    trialShape.Nullify();
    trialCandidate = candidate;
    trialCandidate.boundaryWire.Nullify();
    trialCandidate.boundaryEdges.clear();
    trialCandidate.boundaryVertices.clear();
    trialCandidate.boundaryOwnerFaces.clear();
    trialCandidate.ownerFaces.clear();
    trialCandidate.ownerShell.Nullify();
    if (currentShape.IsNull() || candidate.ownerShell.IsNull()) {
        errorMessage = "cannot create a trial candidate from a null shape or shell";
        return false;
    }

    try {
        // 完整复制 owner shell 的拓扑和几何；试算中的 SameParameter、PCurve 和
        // tolerance 修改全部留在 trial 内，失败时可以真正丢弃。
        BRepBuilderAPI_Copy copier(
            candidate.ownerShell,
            Standard_True,
            Standard_False);
        const TopoDS_Shape copiedOwnerShape = copier.Shape();
        if (copiedOwnerShape.IsNull()
            || copiedOwnerShape.ShapeType() != TopAbs_SHELL) {
            errorMessage = "owner shell trial copy did not produce a shell";
            return false;
        }
        trialCandidate.ownerShell = TopoDS::Shell(copiedOwnerShape);

        const auto mapSubShape = [&copier](
            const TopoDS_Shape& source,
            TopAbs_ShapeEnum expectedType,
            TopoDS_Shape& mapped) {
            mapped = copier.ModifiedShape(source);
            if (mapped.IsNull() || mapped.ShapeType() != expectedType) {
                return false;
            }
            mapped.Orientation(source.Orientation());
            return true;
        };

        BRepBuilderAPI_MakeWire wireBuilder;
        for (const TopoDS_Edge& sourceEdge : candidate.boundaryEdges) {
            TopoDS_Shape mappedShape;
            if (!mapSubShape(sourceEdge, TopAbs_EDGE, mappedShape)) {
                errorMessage = "cannot map a boundary edge to the trial shell";
                return false;
            }
            const TopoDS_Edge mappedEdge = TopoDS::Edge(mappedShape);
            trialCandidate.boundaryEdges.push_back(mappedEdge);
            wireBuilder.Add(mappedEdge);
        }
        if (!wireBuilder.IsDone()) {
            errorMessage = "cannot rebuild the boundary wire on the trial shell";
            return false;
        }
        trialCandidate.boundaryWire = wireBuilder.Wire();

        for (const TopoDS_Vertex& sourceVertex : candidate.boundaryVertices) {
            TopoDS_Shape mappedShape;
            if (!mapSubShape(sourceVertex, TopAbs_VERTEX, mappedShape)) {
                errorMessage = "cannot map a boundary vertex to the trial shell";
                return false;
            }
            trialCandidate.boundaryVertices.push_back(
                TopoDS::Vertex(mappedShape));
        }
        for (const TopoDS_Face& sourceFace : candidate.boundaryOwnerFaces) {
            TopoDS_Shape mappedShape;
            if (!mapSubShape(sourceFace, TopAbs_FACE, mappedShape)) {
                errorMessage = "cannot map a boundary owner face to the trial shell";
                return false;
            }
            trialCandidate.boundaryOwnerFaces.push_back(
                TopoDS::Face(mappedShape));
        }
        for (const TopoDS_Face& sourceFace : candidate.ownerFaces) {
            TopoDS_Shape mappedShape;
            if (!mapSubShape(sourceFace, TopAbs_FACE, mappedShape)) {
                errorMessage = "cannot map an owner face to the trial shell";
                return false;
            }
            trialCandidate.ownerFaces.push_back(TopoDS::Face(mappedShape));
        }

        TopoDS_Shell sourceKey = candidate.ownerShell;
        sourceKey.Orientation(TopAbs_FORWARD);
        TopoDS_Shell replacementValue = trialCandidate.ownerShell;
        replacementValue.Orientation(TopAbs_FORWARD);
        BRepTools_ReShape reshaper;
        reshaper.Replace(sourceKey, replacementValue);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()
            || trialShape.ShapeType() != currentShape.ShapeType()) {
            errorMessage = "trial shell replacement produced an invalid root shape";
            trialShape.Nullify();
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "OCCT failure while creating a trial candidate: "
            + failureText(failure);
        trialShape.Nullify();
        return false;
    }
    return true;
}

bool HealSimpleHole::tryRepairCandidate(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    std::vector<std::string> strategyFailures;

    if (params_.tryPlanarPatch) {
        TopoDS_Shape trialShape;
        Candidate trialCandidate;
        TopoDS_Face patchFace;
        double patchArea = 0.0;
        std::string strategyError;
        if (createTrialCandidate(
                currentShape,
                candidate,
                trialShape,
                trialCandidate,
                strategyError)
            && buildPlanarPatchFace(
                trialCandidate,
                patchFace,
                patchArea,
                strategyError)
            && tryCommitPatchFace(
                trialShape,
                existingFailures,
                trialCandidate,
                patchFace,
                trialCandidate.boundaryEdges,
                TopologyUtils::CountFreeEdges(trialShape),
                repairedShape,
                strategyError)) {
            std::ostringstream success;
            success << "Planar patch accepted: area=" << patchArea
                    << ", freeEdgesBefore="
                    << TopologyUtils::CountFreeEdges(trialShape)
                    << ", freeEdgesAfter="
                    << TopologyUtils::CountFreeEdges(repairedShape);
            appendMessage(success.str());
            return true;
        }
        strategyFailures.push_back(
            RepairDiagnostics::StrategyFailureBlock(
                "PlanarPatch",
                strategyError));
    }

    const auto tryBSplineStrategy = [this,
                                     &currentShape,
                                     &existingFailures,
                                     &candidate,
                                     &repairedShape](
                                        bool useG1,
                                        const FillingComplexityLevel& level,
                                        const StrictC0PlateData* preparedC0Plate,
                                        std::string& strategyError) {
        TopoDS_Shape trialShape;
        Candidate trialCandidate;
        BRepFill_Filling filling;
        TopoDS_Face temporaryPatchFace;
        std::vector<Handle(Geom2d_Curve)> boundaryPCurves;
        double g0Error = 0.0;
        double g1Error = 0.0;
        SurfaceApproximationMetrics approximationMetrics;
        if (useG1) {
            if (!createTrialCandidate(
                    currentShape,
                    candidate,
                    trialShape,
                    trialCandidate,
                    strategyError)
                || !buildTemporaryBSplinePatchFace(
                    trialCandidate,
                    Standard_True,
                    filling,
                    temporaryPatchFace,
                    boundaryPCurves,
                    g0Error,
                    g1Error,
                    strategyError,
                    level.maxDegree,
                    level.maxSegments)) {
                return false;
            }
        }
        else {
            if (preparedC0Plate == nullptr
                || !approximateStrictC0Plate(
                    *preparedC0Plate,
                    params_.surfaceApproximationTolerance,
                    params_.maxSurfaceApproximationDeviation,
                    level.maxDegree,
                    level.maxSegments,
                    params_.maxFillingPolesPerDirection,
                    params_.maxFillingControlPoints,
                    temporaryPatchFace,
                    boundaryPCurves,
                    g0Error,
                    g1Error,
                    approximationMetrics,
                    strategyError)
                || !createTrialCandidate(
                    currentShape,
                    candidate,
                    trialShape,
                    trialCandidate,
                    strategyError)) {
                return false;
            }
        }

        TopoDS_Face patchFace;
        double patchArea = 0.0;
        int reprojectedEdgeCount = 0;
        int reusedPlatePCurveCount = 0;
        double maxSelectedPCurveDeviation = 0.0;
        if (!rebuildBSplinePatchWithOriginalBoundary(
                trialCandidate,
                temporaryPatchFace,
                boundaryPCurves,
                g0Error,
                patchFace,
                patchArea,
                reprojectedEdgeCount,
                reusedPlatePCurveCount,
                maxSelectedPCurveDeviation,
                strategyError)) {
            return false;
        }

        if (!tryCommitPatchFace(
                trialShape,
                existingFailures,
                trialCandidate,
                patchFace,
                trialCandidate.boundaryEdges,
                TopologyUtils::CountFreeEdges(trialShape),
                repairedShape,
                strategyError)) {
            return false;
        }

        std::ostringstream success;
        success << "B-spline "
                << (useG1 ? "G1" : "progressive C0")
                << " patch accepted: area=" << patchArea
                << ", requestedMaxDegree=" << level.maxDegree
                << ", requestedMaxSegments=" << level.maxSegments
                << ", extendedFallback="
                << (level.extendedFallback ? "true" : "false")
                << ", g0Error=" << g0Error;
        if (useG1) {
            success << ", g1Error=" << g1Error;
        }
        else {
            success << ", approximationError="
                    << approximationMetrics.approximationError
                    << ", criterionError="
                    << approximationMetrics.criterionError
                    << ", approximationMethod="
                    << approximationMetrics.method;
        }
        success << ", reprojectedEdges=" << reprojectedEdgeCount
                << ", reusedPlatePCurves=" << reusedPlatePCurveCount
                << ", sampledPCurveDeviation="
                << maxSelectedPCurveDeviation;
        const Handle(Geom_BSplineSurface) patchSurface =
            Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(patchFace));
        if (!patchSurface.IsNull()) {
            success << ", uDegree=" << patchSurface->UDegree()
                    << ", vDegree=" << patchSurface->VDegree()
                    << ", uPoles=" << patchSurface->NbUPoles()
                    << ", vPoles=" << patchSurface->NbVPoles()
                    << ", uC1="
                    << (patchSurface->IsCNu(1) ? "true" : "false")
                    << ", vC1="
                    << (patchSurface->IsCNv(1) ? "true" : "false");
        }
        success << ", freeEdgesAfter="
                << TopologyUtils::CountFreeEdges(repairedShape);
        appendMessage(success.str());
        return true;
    };

    if (params_.tryBSplineG1Patch) {
        const std::vector<FillingComplexityLevel> levels =
            buildComplexityLevels(
                params_.approximationMinDegree,
                params_.approximationMinSegments,
                params_.approximationDegreeStep,
                params_.approximationSegmentStep,
                params_.fillingMaxDegree,
                params_.fillingMaxSegments,
                Standard_False,
                params_.fallbackFillingMaxDegree,
                params_.fallbackFillingMaxSegments,
                params_.fillingMaxSegments);
        for (const FillingComplexityLevel& level : levels) {
            std::string strategyError;
            if (tryBSplineStrategy(
                    Standard_True,
                    level,
                    nullptr,
                    strategyError)) {
                return true;
            }
            std::ostringstream strategyName;
            strategyName << "BSplineG1[d=" << level.maxDegree
                         << ",s=" << level.maxSegments << "]";
            strategyFailures.push_back(
                RepairDiagnostics::StrategyFailureBlock(
                    strategyName.str(),
                    strategyError));
            if (params_.verbose) {
                appendMessage(strategyName.str() + " failed: " + strategyError);
            }
        }
    }

    if (params_.tryBSplineC0Patch) {
        StrictC0PlateData preparedC0Plate;
        std::string preparationError;
        const bool platePrepared = prepareStrictC0Plate(
            candidate.boundaryEdges,
            params_.fillingDegree,
            params_.fillingPointsOnCurve,
            params_.fillingIterations,
            params_.fillingTolerance,
            params_.fillingAngularTolerance,
            params_.fillingCurvatureTolerance,
            params_.maxFillingG0Error,
            preparedC0Plate,
            preparationError);
        if (!platePrepared) {
            strategyFailures.push_back(
                RepairDiagnostics::StrategyFailureBlock(
                    "BSplineC0Plate",
                    preparationError));
        }
        const std::vector<FillingComplexityLevel> levels =
            buildComplexityLevels(
                params_.approximationMinDegree,
                params_.approximationMinSegments,
                params_.approximationDegreeStep,
                params_.approximationSegmentStep,
                params_.fillingMaxDegree,
                params_.fillingMaxSegments,
                params_.tryHighComplexityC0Fallback,
                params_.fallbackFillingMaxDegree,
                params_.fallbackFillingMaxSegments,
                params_.lowestDegreeMaxSegments);
        for (const FillingComplexityLevel& level : levels) {
            if (!platePrepared) {
                break;
            }
            std::string strategyError;
            if (tryBSplineStrategy(
                    Standard_False,
                    level,
                    &preparedC0Plate,
                    strategyError)) {
                return true;
            }
            std::ostringstream strategyName;
            strategyName << "BSplineC0[d=" << level.maxDegree
                         << ",s=" << level.maxSegments << "]";
            strategyFailures.push_back(
                RepairDiagnostics::StrategyFailureBlock(
                    strategyName.str(),
                    strategyError));
            if (params_.verbose) {
                appendMessage(strategyName.str() + " failed: " + strategyError);
            }
        }
    }

    if (params_.useLocalSewingFallback) {
        TopoDS_Shape trialShape;
        Candidate trialCandidate;
        std::string strategyError;
        if (createTrialCandidate(
                currentShape,
                candidate,
                trialShape,
                trialCandidate,
                strategyError)
            && tryRepairByLocalSewing(
                trialShape,
                existingFailures,
                trialCandidate,
                repairedShape,
                strategyError)) {
            return true;
        }
        strategyFailures.push_back(
            RepairDiagnostics::StrategyFailureBlock(
                "LocalSewingFallback",
                strategyError));
    }

    std::ostringstream combinedError;
    for (std::size_t failureIndex = 0;
         failureIndex < strategyFailures.size();
         ++failureIndex) {
        if (failureIndex > 0) {
            combinedError << "\n";
        }
        combinedError << strategyFailures[failureIndex];
    }
    errorMessage = combinedError.str();
    if (errorMessage.empty()) {
        errorMessage = "all hole repair strategies are disabled";
    }
    return false;
}

bool HealSimpleHole::validateRepairResult(
    const TopoDS_Shape& beforeShape,
    const TopoDS_Shape& afterShape,
    const Candidate& candidate,
    const TopoDS_Face& patchFace,
    const std::vector<TopoDS_Edge>& resultBoundaryEdges,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (beforeShape.IsNull() || afterShape.IsNull() || patchFace.IsNull()) {
        errorMessage = "before shape, after shape, or patch face is null";
        return false;
    }

    if (!validateRepairBoundary(
            afterShape, candidate, patchFace, resultBoundaryEdges, freeEdgesBefore, errorMessage)) {
        return false;
    }

    if (params_.rejectInvalidResult) {
        std::string existingError;
        const RepairDiagnostics::FailureList existingFailures =
            RepairDiagnostics::CollectBRepCheckFailures(beforeShape, &existingError);
        if (!existingError.empty()) {
            errorMessage = "failed to collect existing BRepCheck failures: "
                + existingError;
            return false;
        }

        return validateBRepCheckResult(existingFailures, afterShape, errorMessage);
    }
    return true;
}

// 复用本轮错误基线；各候选及各策略的结果仍分别执行完整检查。
bool HealSimpleHole::validateRepairResult(
    const RepairDiagnostics::FailureList& existingFailures,
    const TopoDS_Shape& afterShape,
    const Candidate& candidate,
    const TopoDS_Face& patchFace,
    const std::vector<TopoDS_Edge>& resultBoundaryEdges,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (!validateRepairBoundary(
            afterShape, candidate, patchFace, resultBoundaryEdges, freeEdgesBefore, errorMessage)) {
        return false;
    }
    return !params_.rejectInvalidResult
        || validateBRepCheckResult(existingFailures, afterShape, errorMessage);
}

// 保留原有边界检查顺序，供独立 trial 和本轮缓存两个入口共用。
bool HealSimpleHole::validateRepairBoundary(
    const TopoDS_Shape& afterShape,
    const Candidate& candidate,
    const TopoDS_Face& patchFace,
    const std::vector<TopoDS_Edge>& resultBoundaryEdges,
    int freeEdgesBefore,
    std::string& errorMessage) const
{
    if (afterShape.IsNull() || patchFace.IsNull()) {
        errorMessage = "before shape, after shape, or patch face is null";
        return false;
    }

    if (!shapeContainsFace(afterShape, patchFace)) {
        errorMessage = "result shape does not contain the constructed patch face";
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(afterShape, edgeFaceMap);
    if (resultBoundaryEdges.size() != candidate.boundaryEdges.size()) {
        errorMessage = "result boundary edge count differs from candidate boundary";
        return false;
    }

    for (const TopoDS_Edge& boundaryEdge : resultBoundaryEdges) {
        if (!edgeFaceMap.Contains(boundaryEdge)) {
            errorMessage = "result shape no longer contains a boundary edge";
            return false;
        }

        const TopTools_ListOfShape& adjacentFaces =
            edgeFaceMap.FindFromKey(boundaryEdge);
        const int adjacentFaceCount = listExtent(adjacentFaces);
        if (adjacentFaceCount != 2) {
            std::ostringstream message;
            message << "boundary edge is not shared by exactly two faces after patch: count="
                    << adjacentFaceCount;
            errorMessage = message.str();
            return false;
        }

        bool patchUsesBoundaryEdge = false;
        for (TopTools_ListIteratorOfListOfShape it(adjacentFaces);
             it.More();
             it.Next()) {
            if (it.Value().IsSame(patchFace)) {
                patchUsesBoundaryEdge = true;
                break;
            }
        }
        if (!patchUsesBoundaryEdge) {
            errorMessage = "constructed patch face is not adjacent to a result boundary edge";
            return false;
        }

        // 最终 shell 可能经过 ShapeFix 或 sewing，必须在提交前重新检查补面侧 pcurve。
        const double pcurveCheckTolerance = std::max({
            Precision::Confusion(),
            params_.pcurveTolerance,
            BRep_Tool::Tolerance(boundaryEdge),
            BRep_Tool::Tolerance(patchFace) });
        std::string pcurveError;
        if (!GeometryUtils::HasValidPCurveOnFace(
                boundaryEdge,
                patchFace,
                pcurveCheckTolerance,
                &pcurveError)) {
            errorMessage = "invalid boundary pcurve on committed patch face: "
                + pcurveError;
            return false;
        }
    }

    const int freeEdgesAfter = TopologyUtils::CountFreeEdges(afterShape);
    if (params_.rejectFreeEdgeIncrease) {
        const int expectedMaximum = std::max(
            0,
            freeEdgesBefore - static_cast<int>(resultBoundaryEdges.size()));
        if (freeEdgesAfter > expectedMaximum) {
            std::ostringstream message;
            message << "free edge count did not decrease by the patched boundary size: before="
                    << freeEdgesBefore
                    << ", after=" << freeEdgesAfter
                    << ", boundaryEdges=" << resultBoundaryEdges.size();
            errorMessage = message.str();
            return false;
        }
    }

    return true;
}

// 这里只收集 afterShape，错误增量仍按子形状类型与状态数量比较。
bool HealSimpleHole::validateBRepCheckResult(
    const RepairDiagnostics::FailureList& existingFailures,
    const TopoDS_Shape& afterShape,
    std::string& errorMessage) const
{
    std::string resultError;
    const RepairDiagnostics::FailureList resultFailures =
        RepairDiagnostics::CollectBRepCheckFailures(afterShape, &resultError);
    if (!resultError.empty()) {
        errorMessage = "failed to collect result BRepCheck failures: "
            + resultError;
        return false;
    }

    const RepairDiagnostics::FailureList increasedFailures =
        RepairDiagnostics::FindBRepCheckFailureCountIncrease(
            existingFailures,
            resultFailures);
    if (!increasedFailures.empty()) {
        errorMessage = "new BRepCheck failures after simple-hole patch: "
            + RepairDiagnostics::DescribeBRepCheckFailures(increasedFailures);
        return false;
    }
    return true;
}

void HealSimpleHole::appendMessage(const std::string& message)
{
    if (!message.empty()) {
        messages_.push_back(message);
    }
}
