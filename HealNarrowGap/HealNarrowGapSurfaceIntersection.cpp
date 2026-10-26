#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <BRepAlgoAPI_Section.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools_ReShape.hxx>
#include <Extrema_ExtPC.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Geometry.hxx>
#include <gp_Pnt.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>

#include <TopTools_ListOfShape.hxx>

namespace {
Handle(Geom_Curve) curveForEdgeInWorld(const TopoDS_Edge& edge, double& first, double& last)
{
    first = 0.0;
    last = 0.0;
    if (edge.IsNull()) {
        return Handle(Geom_Curve)();
    }

    TopLoc_Location location;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, location, first, last);
    if (curve.IsNull()) {
        return curve;
    }

    if (!location.IsIdentity()) {
        Handle(Geom_Geometry) transformed = curve->Transformed(location.Transformation());
        curve = Handle(Geom_Curve)::DownCast(transformed);
    }

    return curve;
}

bool projectPointOnCurve(
    const gp_Pnt& point,
    const Handle(Geom_Curve)& curve,
    double first,
    double last,
    double& parameter,
    double& distance)
{
    parameter = first;
    distance = std::numeric_limits<double>::max();
    if (curve.IsNull() || !std::isfinite(first) || !std::isfinite(last) ||
        std::abs(last - first) <= Precision::PConfusion()) {
        return false;
    }

    GeomAPI_ProjectPointOnCurve projector(point, curve, first, last);
    if (projector.NbPoints() <= 0) {
        return false;
    }

    parameter = projector.LowerDistanceParameter();
    distance = projector.LowerDistance();
    return std::isfinite(parameter) && std::isfinite(distance);
}

bool nearestDistanceToEdge(
    const gp_Pnt& point,
    BRepAdaptor_Curve& targetCurve,
    double& nearestDistance)
{
    nearestDistance = std::numeric_limits<double>::max();

    Extrema_ExtPC extrema(
        point,
        targetCurve,
        targetCurve.FirstParameter(),
        targetCurve.LastParameter(),
        Precision::PConfusion());

    if (extrema.IsDone()) {
        for (int extIndex = 1; extIndex <= extrema.NbExt(); ++extIndex) {
            const double squareDistance = extrema.SquareDistance(extIndex);
            if (squareDistance < nearestDistance * nearestDistance) {
                nearestDistance = std::sqrt(std::max(0.0, squareDistance));
            }
        }

        Standard_Real firstSquareDistance = 0.0;
        Standard_Real lastSquareDistance = 0.0;
        gp_Pnt firstPoint;
        gp_Pnt lastPoint;
        extrema.TrimmedSquareDistances(firstSquareDistance, lastSquareDistance, firstPoint, lastPoint);
        nearestDistance = std::min(nearestDistance, std::sqrt(std::max(0.0, firstSquareDistance)));
        nearestDistance = std::min(nearestDistance, std::sqrt(std::max(0.0, lastSquareDistance)));
    }

    if (!std::isfinite(nearestDistance) || nearestDistance == std::numeric_limits<double>::max()) {
        constexpr int fallbackSamples = 32;
        for (int sample = 0; sample < fallbackSamples; ++sample) {
            const double t = static_cast<double>(sample) / static_cast<double>(fallbackSamples - 1);
            const double param = targetCurve.FirstParameter()
                + (targetCurve.LastParameter() - targetCurve.FirstParameter()) * t;
            const double distance = point.Distance(targetCurve.Value(param));
            nearestDistance = std::min(nearestDistance, distance);
        }
    }

    return std::isfinite(nearestDistance) && nearestDistance < std::numeric_limits<double>::max();
}

} // namespace

bool HealNarrowGap::buildSurfaceIntersectionEdgeForCandidate(
    const Candidate& candidate,
    double endpointDistanceTolerance,
    TopoDS_Edge& intersectionEdge,
    std::string& errorMessage) const
{
    intersectionEdge.Nullify();
    if (candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
        candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "invalid candidate for surface intersection";
        return false;
    }

    try {
        TopoDS_Vertex aStart;
        TopoDS_Vertex aEnd;
        TopoDS_Vertex bStart;
        TopoDS_Vertex bEnd;
        if (!getEdgeVerticesOriented(candidate.edgeA, aStart, aEnd) ||
            !getEdgeVerticesOriented(candidate.edgeB, bStart, bEnd)) {
            errorMessage = "failed to get oriented edge vertices for surface intersection";
            return false;
        }

        if (!((aStart.IsSame(bStart) && aEnd.IsSame(bEnd)) ||
              (aStart.IsSame(bEnd) && aEnd.IsSame(bStart)))) {
            errorMessage = "candidate edges do not share endpoints for surface intersection";
            return false;
        }

        BRepAlgoAPI_Section section(candidate.faceA, candidate.faceB, Standard_False);
        if (myParams.surfaceIntersectionFuzzyTol > 0.0) {
            section.SetFuzzyValue(myParams.surfaceIntersectionFuzzyTol);
        }
        section.Approximation(Standard_True);
        section.ComputePCurveOn1(Standard_True);
        section.ComputePCurveOn2(Standard_True);
        section.Build();
        if (!section.IsDone() || section.HasErrors()) {
            errorMessage = "BRepAlgoAPI_Section failed";
            return false;
        }

        const TopoDS_Shape sectionShape = section.Shape();
        if (sectionShape.IsNull()) {
            errorMessage = "surface intersection returned null shape";
            return false;
        }

        BRepAdaptor_Curve edgeCurveA(candidate.edgeA);
        BRepAdaptor_Curve edgeCurveB(candidate.edgeB);
        const gp_Pnt startPoint = BRep_Tool::Pnt(aStart);
        const gp_Pnt endPoint = BRep_Tool::Pnt(aEnd);
        const double endpointTol = std::max({
            Precision::Confusion(),
            endpointDistanceTolerance,
            BRep_Tool::Tolerance(aStart),
            BRep_Tool::Tolerance(aEnd)
            });
        const int sampleCount = std::max(8, std::min(32, computeSampleCountForCandidate(candidate)));

        bool found = false;
        double bestScore = std::numeric_limits<double>::max();
        TopoDS_Edge bestEdge;
        int sectionEdgeCount = 0;
        int rejectedByEndpointDistance = 0;
        int rejectedByGapDistance = 0;
        int rejectedByMakeEdge = 0;
        double bestRejectedStartDistance = std::numeric_limits<double>::max();
        double bestRejectedEndDistance = std::numeric_limits<double>::max();
        double bestRejectedEndpointDistance = std::numeric_limits<double>::max();
        double bestRejectedGapDistance = std::numeric_limits<double>::max();

        for (TopExp_Explorer edgeExp(sectionShape, TopAbs_EDGE); edgeExp.More(); edgeExp.Next()) {
            const TopoDS_Edge sectionEdge = TopoDS::Edge(edgeExp.Current());
            if (sectionEdge.IsNull() || isDegeneratedEdge(sectionEdge) || !has3DCurve(sectionEdge)) {
                continue;
            }

            double first = 0.0;
            double last = 0.0;
            Handle(Geom_Curve) curve = curveForEdgeInWorld(sectionEdge, first, last);
            if (curve.IsNull()) {
                continue;
            }
            ++sectionEdgeCount;

            double startParam = first;
            double endParam = last;
            double startDistance = 0.0;
            double endDistance = 0.0;
            if (!projectPointOnCurve(startPoint, curve, first, last, startParam, startDistance) ||
                !projectPointOnCurve(endPoint, curve, first, last, endParam, endDistance)) {
                continue;
            }

            if (std::abs(endParam - startParam) <= Precision::PConfusion()) {
                continue;
            }

            // 构造 edge 前先校验原端点是否真正落在 section curve 上。
            // 否则 MakeEdge 可能成功，但后续 pcurve 校验会在端点或近端采样处失败。
            const double endpointDistance = std::max(startDistance, endDistance);
            if (endpointDistance > endpointTol) {
                ++rejectedByEndpointDistance;
                if (endpointDistance < bestRejectedEndpointDistance) {
                    bestRejectedEndpointDistance = endpointDistance;
                    bestRejectedStartDistance = startDistance;
                    bestRejectedEndDistance = endDistance;
                }
                continue;
            }

            double maxDistanceToGapEdges = 0.0;
            double meanDistanceToGapEdges = 0.0;
            bool distanceFailed = false;
            for (int sample = 0; sample < sampleCount; ++sample) {
                const double t = static_cast<double>(sample) / static_cast<double>(sampleCount - 1);
                const double param = startParam + (endParam - startParam) * t;
                const gp_Pnt point = curve->Value(param);

                double distanceToA = 0.0;
                double distanceToB = 0.0;
                if (!nearestDistanceToEdge(point, edgeCurveA, distanceToA) ||
                    !nearestDistanceToEdge(point, edgeCurveB, distanceToB)) {
                    distanceFailed = true;
                    break;
                }

                const double distanceToGap = std::min(distanceToA, distanceToB);
                maxDistanceToGapEdges = std::max(maxDistanceToGapEdges, distanceToGap);
                meanDistanceToGapEdges += distanceToGap;
            }

            if (distanceFailed) {
                continue;
            }
            meanDistanceToGapEdges /= static_cast<double>(sampleCount);

            // section edge 必须靠近原 gap 区域；否则复杂 face-face 相交时可能选到其他交线。
            if (maxDistanceToGapEdges > myParams.surfaceIntersectionMaxGap) {
                ++rejectedByGapDistance;
                bestRejectedGapDistance = std::min(bestRejectedGapDistance, maxDistanceToGapEdges);
                continue;
            }

            Handle(Geom_Curve) edgeCurve = curve;
            double edgeStartParam = startParam;
            double edgeEndParam = endParam;
            if (edgeEndParam < edgeStartParam) {
                // 保证新 edge 的底层曲线方向和候选起点/终点一致，后续 face use 再沿用原 wire 方向。
                edgeCurve = curve->Reversed();
                edgeStartParam = curve->ReversedParameter(startParam);
                edgeEndParam = curve->ReversedParameter(endParam);
            }

            BRepBuilderAPI_MakeEdge edgeMaker(edgeCurve, aStart, aEnd, edgeStartParam, edgeEndParam);
            if (!edgeMaker.IsDone()) {
                // 如果原顶点不在交线上，说明交线不能在不重建相邻边的情况下直接用于当前 wire。
                ++rejectedByMakeEdge;
                continue;
            }

            const double endpointScore = std::max(startDistance, endDistance);
            const double score =
                endpointScore * 4.0 +
                meanDistanceToGapEdges +
                maxDistanceToGapEdges * 0.25;
            if (score < bestScore) {
                bestScore = score;
                bestEdge = edgeMaker.Edge();
                found = true;
            }
        }

        if (!found || bestEdge.IsNull()) {
            std::ostringstream os;
            os << "no usable surface intersection edge close to candidate gap"
               << ": sectionEdges=" << sectionEdgeCount
               << ", endpointDistanceTolerance=" << endpointTol
               << ", rejectedByEndpointDistance=" << rejectedByEndpointDistance
               << ", rejectedByGapDistance=" << rejectedByGapDistance
               << ", rejectedByMakeEdge=" << rejectedByMakeEdge;
            if (bestRejectedEndpointDistance < std::numeric_limits<double>::max()) {
                os << ", bestRejectedStartDistance=" << bestRejectedStartDistance
                   << ", bestRejectedEndDistance=" << bestRejectedEndDistance;
            }
            if (bestRejectedGapDistance < std::numeric_limits<double>::max()) {
                os << ", bestRejectedGapDistance=" << bestRejectedGapDistance;
            }
            errorMessage = os.str();
            return false;
        }

        intersectionEdge = bestEdge;
    }
    catch (Standard_Failure& failure) {
        errorMessage = "surface intersection edge construction failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealNarrowGap::tryRepairBySurfaceIntersection(
    const TopoDS_Shape& currentShape,
    const RepairDiagnostics::FailureList& existingFailures,
    const Candidate& candidate,
    TopoDS_Shape& repairedShape,
    std::string& errorMessage)
{
    repairedShape.Nullify();
    if (currentShape.IsNull()) {
        errorMessage = "Current shape is null";
        return false;
    }
    if (candidate.faceA.IsNull() || candidate.faceB.IsNull() ||
        candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "Invalid candidate for surface intersection";
        return false;
    }
    if (candidate.faceA.IsSame(candidate.faceB)) {
        errorMessage = "surface intersection does not support two edges on the same face";
        return false;
    }
    if (candidate.maxGap > myParams.surfaceIntersectionMaxGap) {
        std::ostringstream os;
        os << "surface intersection gap too large: maxGap=" << candidate.maxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
           << ", candidateMaxGap=" << myParams.candidateMaxGap;
        errorMessage = os.str();
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(currentShape);
    // 新 edge 的 pcurve 精度至少为配置值，并保留原两条 edge 的最大 tolerance。
    const double pcurveTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(candidate.edgeA),
        BRep_Tool::Tolerance(candidate.edgeB),
        myParams.pcurveTol
        });

    const double existingFaceTol = std::max({
        Precision::Confusion(),
        BRep_Tool::Tolerance(candidate.faceA),
        BRep_Tool::Tolerance(candidate.faceB)
        });
    const double faceTol = std::min(
        std::max(existingFaceTol, myParams.faceTol),
        myParams.maxFaceTol);

    TopoDS_Edge intersectionEdge;
    if (!buildSurfaceIntersectionEdgeForCandidate(
            candidate,
            myParams.maxEndpointTolerance,
            intersectionEdge,
            errorMessage)) {
        return false;
    }

    TopoDS_Edge edgeAInWire;
    if (!findEdgeUseInFaceWire(candidate.faceA, candidate.edgeA, edgeAInWire, errorMessage)) {
        errorMessage = "edgeA not found in faceA wire for surface intersection: " + errorMessage;
        return false;
    }

    TopoDS_Edge edgeBInWire;
    if (!findEdgeUseInFaceWire(candidate.faceB, candidate.edgeB, edgeBInWire, errorMessage)) {
        errorMessage = "edgeB not found in faceB wire for surface intersection: " + errorMessage;
        return false;
    }

    // 新 edge 的使用方向跟随被替换 edge 在各自 wire 中的方向，保持 face 内侧判断依赖的 wire 方向关系。
    TopoDS_Edge intersectionUseA = intersectionEdge;
    intersectionUseA.Orientation(edgeAInWire.Orientation());
    TopoDS_Edge intersectionUseB = intersectionEdge;
    intersectionUseB.Orientation(edgeBInWire.Orientation());

    if (!ensurePCurveOnFace(intersectionUseA, candidate.faceA, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add surface intersection pcurve on faceA: " + errorMessage;
        return false;
    }
    if (!ensurePCurveOnFace(intersectionUseB, candidate.faceB, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add surface intersection pcurve on faceB: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceA;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceA,
            candidate.edgeA,
            intersectionUseA,
            faceTol,
            rebuiltFaceA,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceA with surface intersection edge: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceB;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceB,
            candidate.edgeB,
            intersectionUseB,
            faceTol,
            rebuiltFaceB,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceB with surface intersection edge: " + errorMessage;
        return false;
    }

    BRepCheck_Analyzer analyzerFaceA(rebuiltFaceA);
    if (!analyzerFaceA.IsValid()) {
        appendMessage(
            messages_,
            std::string("Warning: intermediate rebuilt faceA is invalid before final validation:")
                + RepairDiagnostics::DescribeBRepCheckFailures(analyzerFaceA, rebuiltFaceA));
    }
    BRepCheck_Analyzer analyzerFaceB(rebuiltFaceB);
    if (!analyzerFaceB.IsValid()) {
        appendMessage(
            messages_,
            std::string("Warning: intermediate rebuilt faceB is invalid before final validation:")
                + RepairDiagnostics::DescribeBRepCheckFailures(analyzerFaceB, rebuiltFaceB));
    }

    TopoDS_Shape trialShape;
    try {
        BRepTools_ReShape reshaper;
        reshaper.Replace(candidate.faceA, rebuiltFaceA);
        reshaper.Replace(candidate.faceB, rebuiltFaceB);
        trialShape = reshaper.Apply(currentShape);
        if (trialShape.IsNull()) {
            errorMessage = "Replaced shape is null after surface intersection";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Surface intersection shape replacement failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape repairedEdgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(trialShape, repairedEdgeFaceMap);
    TopTools_ListOfShape intersectionFaces;
    if (!findFacesForSameEdge(repairedEdgeFaceMap, intersectionEdge, intersectionFaces)) {
        errorMessage = "surface intersection edge not found after shape replacement";
        return false;
    }

    std::string pcurveRepairError;
    for (TopTools_ListIteratorOfListOfShape it(intersectionFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        TopoDS_Edge intersectionUseInFace;
        if (!findEdgeUseInFaceWire(face, intersectionEdge, intersectionUseInFace, pcurveRepairError)) {
            errorMessage = "Failed to find surface intersection edge use after replacement: " + pcurveRepairError;
            return false;
        }

        if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                intersectionUseInFace,
                face,
                pcurveTol,
                &pcurveRepairError)) {
            errorMessage = "Failed to repair surface intersection pcurve after replacement: "
                + pcurveRepairError;
            return false;
        }
    }

    if (myParams.runShapeFixAfterRepair) {
        try {
            const double fixPrecision = std::max(Precision::Confusion(), myParams.shapeFixPrecision);
            const double fixMaxTol = std::max(fixPrecision, myParams.shapeFixMaxTolerance);
            ShapeFix_Shape fixer(trialShape);
            fixer.SetPrecision(fixPrecision);
            fixer.SetMaxTolerance(fixMaxTol);
            fixer.FixWireTool()->FixEdgeCurvesMode() = 1;
            fixer.FixFaceTool()->FixMissingSeamMode() = 1;
            fixer.Perform();
            TopoDS_Shape fixedShape = fixer.Shape();
            if (!fixedShape.IsNull()) {
                trialShape = fixedShape;
            }
        }
        catch (Standard_Failure& failure) {
            errorMessage = "ShapeFix failed after surface intersection: " + RepairDiagnostics::FailureText(failure);
            return false;
        }
    }

    if (!validateEdgeReplacement(
            existingFailures,
            trialShape,
            intersectionEdge,
            freeEdgesBefore,
            pcurveTol,
            errorMessage)) {
        errorMessage = "surface intersection validation failed: " + errorMessage;
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by surface intersection: candidateMaxGap=" << myParams.candidateMaxGap
           << ", replaceOneSideMaxGap=" << myParams.replaceOneSideMaxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
           << ", surfaceIntersectionFuzzyTol=" << myParams.surfaceIntersectionFuzzyTol
           << ", midlineMaxGap=" << myParams.midlineMaxGap
           << ", candidate.maxGap=" << candidate.maxGap
           << ", candidate.meanGap=" << candidate.meanGap
           << ", sampleCount=" << candidate.usedSampleCount
           << ", effectivePCurveTol=" << pcurveTol
           << ", faceTol=" << faceTol
           << ", shapeFixPrecision=" << myParams.shapeFixPrecision
           << ", shapeFixMaxTolerance=" << myParams.shapeFixMaxTolerance;
        appendMessage(messages_, os.str());
    }

    repairedShape = trialShape;
    return true;
}

