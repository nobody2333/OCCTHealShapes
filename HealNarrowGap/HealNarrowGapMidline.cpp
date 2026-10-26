#include "HealNarrowGap.h"
#include "GeometryUtils.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools_ReShape.hxx>
#include <Extrema_ExtPC.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <Geom_BSplineCurve.hxx>
#include <gp_Pnt.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TColgp_Array1OfPnt.hxx>

#include <TopTools_ListOfShape.hxx>

bool HealNarrowGap::buildMidlineEdgeForCandidate(
    const Candidate& candidate,
    TopoDS_Edge& midlineEdge,
    std::string& errorMessage) const
{
    midlineEdge.Nullify();
    if (candidate.edgeA.IsNull() || candidate.edgeB.IsNull()) {
        errorMessage = "candidate edge is null";
        return false;
    }

    try {
        // 从较可信的一侧按弧长取点，再在另一条裁剪后的 3D edge 上找最近点。
        // 不用两条曲线各自的相同比例参数推断几何对应关系。
        const EdgeChoice sourceChoice = chooseEdgeToKeep(candidate);
        const TopoDS_Edge& sourceEdge = sourceChoice.keptEdge;
        const TopoDS_Edge& targetEdge = sourceChoice.edgeToReplace;
        TopoDS_Vertex sourceStart;
        TopoDS_Vertex sourceEnd;
        TopoDS_Vertex targetStart;
        TopoDS_Vertex targetEnd;
        if (!getEdgeVerticesOriented(sourceEdge, sourceStart, sourceEnd) ||
            !getEdgeVerticesOriented(targetEdge, targetStart, targetEnd)) {
            errorMessage = "failed to get oriented edge vertices for midline";
            return false;
        }

        TopoDS_Vertex matchingTargetStart;
        TopoDS_Vertex matchingTargetEnd;
        if (sourceStart.IsSame(targetStart) && sourceEnd.IsSame(targetEnd)) {
            matchingTargetStart = targetStart;
            matchingTargetEnd = targetEnd;
        }
        else if (sourceStart.IsSame(targetEnd) && sourceEnd.IsSame(targetStart)) {
            matchingTargetStart = targetEnd;
            matchingTargetEnd = targetStart;
        }
        else {
            errorMessage = "candidate edges do not share endpoints for midline";
            return false;
        }

        const double sourceFirst = BRep_Tool::Parameter(sourceStart, sourceEdge);
        const double sourceLast = BRep_Tool::Parameter(sourceEnd, sourceEdge);
        const double targetFirst = BRep_Tool::Parameter(matchingTargetStart, targetEdge);
        const double targetLast = BRep_Tool::Parameter(matchingTargetEnd, targetEdge);
        if (!std::isfinite(sourceFirst) || !std::isfinite(sourceLast) ||
            !std::isfinite(targetFirst) || !std::isfinite(targetLast) ||
            std::abs(sourceLast - sourceFirst) <= Precision::PConfusion() ||
            std::abs(targetLast - targetFirst) <= Precision::PConfusion()) {
            errorMessage = "candidate edge has an invalid curve parameter range for midline";
            return false;
        }

        BRepAdaptor_Curve sourceCurve(sourceEdge);
        BRepAdaptor_Curve targetCurve(targetEdge);
        const int sampleCount = std::max(3, computeSampleCountForCandidate(candidate));
        GCPnts_UniformAbscissa sourceSampler(
            sourceCurve,
            sampleCount,
            std::min(sourceFirst, sourceLast),
            std::max(sourceFirst, sourceLast),
            Precision::Confusion());
        if (!sourceSampler.IsDone() || sourceSampler.NbPoints() != sampleCount) {
            errorMessage = "failed to sample source edge uniformly by arc length";
            return false;
        }

        const double targetDirection = targetLast > targetFirst ? 1.0 : -1.0;
        const double targetParameterSpan = std::abs(targetLast - targetFirst);
        const double progressTolerance = std::min(
            Precision::PConfusion(),
            targetParameterSpan / (10.0 * sampleCount));
        double previousTargetProgress = 0.0;

        std::vector<gp_Pnt> midPoints;
        midPoints.reserve(static_cast<std::size_t>(sampleCount));
        for (int i = 0; i < sampleCount; ++i) {
            gp_Pnt point;
            if (i == 0) {
                // 端点沿用原拓扑顶点，保证新 edge 能与原 wire 中相邻 edge 精确连接。
                point = BRep_Tool::Pnt(sourceStart);
            }
            else if (i + 1 == sampleCount) {
                point = BRep_Tool::Pnt(sourceEnd);
            }
            else {
                const int sampleIndex = sourceFirst <= sourceLast
                    ? i + 1
                    : sampleCount - i;
                const gp_Pnt sourcePoint = sourceCurve.Value(
                    sourceSampler.Parameter(sampleIndex));
                Extrema_ExtPC projection(
                    sourcePoint,
                    targetCurve,
                    std::min(targetFirst, targetLast),
                    std::max(targetFirst, targetLast),
                    Precision::PConfusion());
                if (!projection.IsDone()) {
                    errorMessage = "midline projection onto the target edge did not finish";
                    return false;
                }

                double nearestDistance = std::numeric_limits<double>::infinity();
                double forwardDistance = std::numeric_limits<double>::infinity();
                double forwardProgress = 0.0;
                gp_Pnt forwardPoint;
                const auto considerProjection = [&](double parameter,
                                                    const gp_Pnt& projectedPoint) {
                    const double distance = sourcePoint.Distance(projectedPoint);
                    if (!std::isfinite(parameter) || !std::isfinite(distance)) {
                        return;
                    }
                    nearestDistance = std::min(nearestDistance, distance);
                    const double progress = targetDirection * (parameter - targetFirst);
                    if (progress <= previousTargetProgress + progressTolerance ||
                        progress >= targetParameterSpan - progressTolerance ||
                        distance >= forwardDistance) {
                        return;
                    }
                    forwardDistance = distance;
                    forwardProgress = progress;
                    forwardPoint = projectedPoint;
                };

                for (int solution = 1; solution <= projection.NbExt(); ++solution) {
                    const double parameter = projection.Point(solution).Parameter();
                    considerProjection(parameter, targetCurve.Value(parameter));
                }
                considerProjection(targetFirst, targetCurve.Value(targetFirst));
                considerProjection(targetLast, targetCurve.Value(targetLast));

                if (!std::isfinite(forwardDistance) ||
                    forwardDistance > nearestDistance + Precision::Confusion()) {
                    errorMessage = "nearest midline projection reverses or stalls on the target edge";
                    return false;
                }
                if (forwardDistance > myParams.midlineMaxGap) {
                    errorMessage = "projected midline point exceeds midlineMaxGap";
                    return false;
                }
                previousTargetProgress = forwardProgress;
                point.SetCoord(
                    0.5 * (sourcePoint.X() + forwardPoint.X()),
                    0.5 * (sourcePoint.Y() + forwardPoint.Y()),
                    0.5 * (sourcePoint.Z() + forwardPoint.Z()));
            }

            if (!midPoints.empty() &&
                point.Distance(midPoints.back()) <= Precision::Confusion()) {
                if (i + 1 == sampleCount) {
                    errorMessage = "midline edge has zero-length end segment";
                    return false;
                }
                continue;
            }
            midPoints.push_back(point);
        }

        if (midPoints.size() < 2) {
            errorMessage = "midline edge has insufficient points";
            return false;
        }

        if (midPoints.size() == 2) {
            BRepBuilderAPI_MakeEdge edgeMaker(sourceStart, sourceEnd);
            if (!edgeMaker.IsDone()) {
                errorMessage = "failed to create straight midline edge";
                return false;
            }
            midlineEdge = edgeMaker.Edge();
        }
        else {
            TColgp_Array1OfPnt points(1, static_cast<Standard_Integer>(midPoints.size()));
            for (std::size_t i = 0; i < midPoints.size(); ++i) {
                points.SetValue(static_cast<Standard_Integer>(i + 1), midPoints[i]);
            }

            GeomAPI_PointsToBSpline splineBuilder(points);
            Handle(Geom_BSplineCurve) midlineCurve = splineBuilder.Curve();
            if (midlineCurve.IsNull()) {
                errorMessage = "failed to interpolate midline curve";
                return false;
            }

            BRepBuilderAPI_MakeEdge edgeMaker(midlineCurve, sourceStart, sourceEnd);
            if (!edgeMaker.IsDone()) {
                errorMessage = "failed to create midline edge";
                return false;
            }
            midlineEdge = edgeMaker.Edge();
        }

        if (midlineEdge.IsNull()) {
            errorMessage = "midline edge is null";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "midline edge construction failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealNarrowGap::tryRepairByMidlineCollapse(
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
        errorMessage = "Invalid candidate for midline collapse";
        return false;
    }
    if (candidate.faceA.IsSame(candidate.faceB)) {
        errorMessage = "midline collapse does not support two edges on the same face";
        return false;
    }
    if (candidate.maxGap > myParams.midlineMaxGap) {
        std::ostringstream os;
        os << "midline gap too large: maxGap=" << candidate.maxGap
           << ", midlineMaxGap=" << myParams.midlineMaxGap
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

    TopoDS_Edge midlineEdge;
    if (!buildMidlineEdgeForCandidate(candidate, midlineEdge, errorMessage)) {
        return false;
    }

    TopoDS_Edge edgeAInWire;
    if (!findEdgeUseInFaceWire(candidate.faceA, candidate.edgeA, edgeAInWire, errorMessage)) {
        errorMessage = "edgeA not found in faceA wire for midline collapse: " + errorMessage;
        return false;
    }

    TopoDS_Edge edgeBInWire;
    if (!findEdgeUseInFaceWire(candidate.faceB, candidate.edgeB, edgeBInWire, errorMessage)) {
        errorMessage = "edgeB not found in faceB wire for midline collapse: " + errorMessage;
        return false;
    }

    // 新 edge 的使用方向跟随被替换 edge 在各自 wire 中的方向，而不是重新用几何端点推断。
    TopoDS_Edge midlineUseA = midlineEdge;
    midlineUseA.Orientation(edgeAInWire.Orientation());
    TopoDS_Edge midlineUseB = midlineEdge;
    midlineUseB.Orientation(edgeBInWire.Orientation());

    if (!ensurePCurveOnFace(midlineUseA, candidate.faceA, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add midline pcurve on faceA: " + errorMessage;
        return false;
    }
    if (!ensurePCurveOnFace(midlineUseB, candidate.faceB, pcurveTol, errorMessage)) {
        errorMessage = "Failed to add midline pcurve on faceB: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceA;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceA,
            candidate.edgeA,
            midlineUseA,
            faceTol,
            rebuiltFaceA,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceA with midline edge: " + errorMessage;
        return false;
    }

    TopoDS_Face rebuiltFaceB;
    if (!rebuildFaceByReplacingEdge(
            candidate.faceB,
            candidate.edgeB,
            midlineUseB,
            faceTol,
            rebuiltFaceB,
            errorMessage)) {
        errorMessage = "Failed to rebuild faceB with midline edge: " + errorMessage;
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
            errorMessage = "Replaced shape is null after midline collapse";
            return false;
        }
    }
    catch (Standard_Failure& failure) {
        errorMessage = "Midline shape replacement failed: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    TopTools_IndexedDataMapOfShapeListOfShape repairedEdgeFaceMap;
    TopologyUtils::BuildEdgeFaceMap(trialShape, repairedEdgeFaceMap);
    TopTools_ListOfShape midlineFaces;
    if (!findFacesForSameEdge(repairedEdgeFaceMap, midlineEdge, midlineFaces)) {
        errorMessage = "midline edge not found after shape replacement";
        return false;
    }

    std::string pcurveRepairError;
    for (TopTools_ListIteratorOfListOfShape it(midlineFaces); it.More(); it.Next()) {
        const TopoDS_Face face = TopoDS::Face(it.Value());
        if (face.IsNull()) {
            continue;
        }

        TopoDS_Edge midlineUseInFace;
        if (!findEdgeUseInFaceWire(face, midlineEdge, midlineUseInFace, pcurveRepairError)) {
            errorMessage = "Failed to find midline edge use after replacement: " + pcurveRepairError;
            return false;
        }

        if (!GeometryUtils::RepairPCurveForEdgeOnFace(
                midlineUseInFace,
                face,
                pcurveTol,
                &pcurveRepairError)) {
            errorMessage = "Failed to repair midline pcurve after replacement: " + pcurveRepairError;
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
            errorMessage = "ShapeFix failed after midline collapse: " + RepairDiagnostics::FailureText(failure);
            return false;
        }
    }

    if (!validateEdgeReplacement(
            existingFailures,
            trialShape,
            midlineEdge,
            freeEdgesBefore,
            pcurveTol,
            errorMessage)) {
        errorMessage = "midline collapse validation failed: " + errorMessage;
        return false;
    }

    {
        std::ostringstream os;
        os << "Candidate repair accepted by midline collapse: candidateMaxGap=" << myParams.candidateMaxGap
           << ", replaceOneSideMaxGap=" << myParams.replaceOneSideMaxGap
           << ", surfaceIntersectionMaxGap=" << myParams.surfaceIntersectionMaxGap
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


