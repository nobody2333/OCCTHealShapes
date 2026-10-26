#include "HealClosedNarrowGap.h"
#include "RepairDiagnostics.h"
#include "TopologyUtils.h"

#include "GeometryUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_Curve.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS.hxx>
#include <TColgp_HArray1OfPnt.hxx>

namespace {
constexpr double kStationMergeTolerance =
    128.0 * std::numeric_limits<double>::epsilon();

} // namespace

bool HealClosedNarrowGap::isCandidateEndConsistent(
    const GapEnd& end,
    const TopoDS_Vertex& vertexA,
    const TopoDS_Vertex& vertexB) const
{
    if (vertexA.IsNull() || vertexB.IsNull() ||
        end.vertexA.IsNull() || end.vertexB.IsNull() ||
        !end.vertexA.IsSame(vertexA) ||
        !end.vertexB.IsSame(vertexB)) {
        return false;
    }

    if (end.type == GapEndType::SharedTip) {
        return vertexA.IsSame(vertexB) && end.connectorEdge.IsNull();
    }
    if (end.type != GapEndType::ShortEdge ||
        vertexA.IsSame(vertexB) ||
        end.connectorEdge.IsNull() ||
        end.connectorOwnerFace.IsNull()) {
        return false;
    }

    TopoDS_Vertex connectorFirst;
    TopoDS_Vertex connectorLast;
    if (!TopologyUtils::GetEdgeVertices(
            end.connectorEdge,
            connectorFirst,
            connectorLast)) {
        return false;
    }
    return (connectorFirst.IsSame(vertexA) &&
            connectorLast.IsSame(vertexB)) ||
        (connectorFirst.IsSame(vertexB) &&
            connectorLast.IsSame(vertexA));
}

bool HealClosedNarrowGap::buildTargetEndVertex(
    const GapEnd& end,
    TopoDS_Vertex& targetVertex,
    std::string& errorMessage) const
{
    targetVertex.Nullify();
    if (end.vertexA.IsNull() || end.vertexB.IsNull()) {
        errorMessage = "target endpoint construction failure: source endpoint is null";
        return false;
    }

    if (end.type == GapEndType::SharedTip) {
        if (!end.vertexA.IsSame(end.vertexB) || !end.connectorEdge.IsNull()) {
            errorMessage = "target endpoint construction failure: shared tip topology is inconsistent";
            return false;
        }
        targetVertex = end.vertexA;
        return true;
    }
    if (end.type != GapEndType::ShortEdge ||
        end.vertexA.IsSame(end.vertexB) ||
        end.connectorEdge.IsNull() ||
        end.connectorOwnerFace.IsNull()) {
        errorMessage = "target endpoint construction failure: short end edge topology is inconsistent";
        return false;
    }

    try {
        const gp_Pnt pointA = BRep_Tool::Pnt(end.vertexA);
        const gp_Pnt pointB = BRep_Tool::Pnt(end.vertexB);
        if (!GeometryUtils::IsFinitePoint(pointA) || !GeometryUtils::IsFinitePoint(pointB)) {
            errorMessage = "target endpoint construction failure: short end edge has a non-finite endpoint";
            return false;
        }

        const gp_Pnt midpoint(
            0.5 * (pointA.X() + pointB.X()),
            0.5 * (pointA.Y() + pointB.Y()),
            0.5 * (pointA.Z() + pointB.Z()));
        if (!GeometryUtils::IsFinitePoint(midpoint)) {
            errorMessage = "target endpoint construction failure: short end midpoint is not finite";
            return false;
        }

        // 两个短边端点都会替换为同一个 vertex，容差只覆盖已有误差和实际移动距离。
        // 不能混入 face/pcurve 修复容差，否则后续组 wire 时可能把附近的其他顶点误认为连接点。
        const double vertexTolerance = std::max({
            Precision::Confusion(),
            BRep_Tool::Tolerance(end.vertexA),
            BRep_Tool::Tolerance(end.vertexB),
            BRep_Tool::Tolerance(end.connectorEdge),
            pointA.Distance(midpoint) + Precision::Confusion(),
            pointB.Distance(midpoint) + Precision::Confusion()
        });

        BRep_Builder builder;
        builder.MakeVertex(targetVertex, midpoint, vertexTolerance);
        if (targetVertex.IsNull()) {
            errorMessage = "target endpoint construction failure: failed to create short end midpoint vertex";
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "target endpoint construction failure: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::buildMidlineChain(
    const ClosedGapCandidate& candidate,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    midlineData = MidlineChainData();
    if (candidate.chainA.edges.empty() || candidate.chainB.edges.empty()) {
        errorMessage = "candidate chains are empty";
        return false;
    }
    if (candidate.chainA.vertices.size() != candidate.chainA.edges.size() + 1 ||
        candidate.chainB.vertices.size() != candidate.chainB.edges.size() + 1) {
        errorMessage = "candidate chain vertex count is invalid";
        return false;
    }
    if (candidate.chainA.ownerFaces.size() != candidate.chainA.edges.size() ||
        candidate.chainB.ownerFaces.size() != candidate.chainB.edges.size()) {
        errorMessage = "candidate chain owner-face count is invalid";
        return false;
    }
    if (!std::isfinite(candidate.chainA.length) ||
        !std::isfinite(candidate.chainB.length) ||
        candidate.chainA.length <= Precision::Confusion() ||
        candidate.chainB.length <= Precision::Confusion()) {
        errorMessage = "candidate chain total length is invalid";
        return false;
    }
    if (!isCandidateEndConsistent(
            candidate.startEnd,
            candidate.chainA.vertices.front(),
            candidate.chainB.vertices.front()) ||
        !isCandidateEndConsistent(
            candidate.finishEnd,
            candidate.chainA.vertices.back(),
            candidate.chainB.vertices.back())) {
        errorMessage = "candidate chain endpoints do not match their endpoint definitions";
        return false;
    }
    for (const TopoDS_Face& ownerFace : candidate.chainA.ownerFaces) {
        if (ownerFace.IsNull()) {
            errorMessage = "chainA contains an edge without owner face";
            return false;
        }
    }
    for (const TopoDS_Face& ownerFace : candidate.chainB.ownerFaces) {
        if (ownerFace.IsNull()) {
            errorMessage = "chainB contains an edge without owner face";
            return false;
        }
    }

    MidlineConstructionMode constructionMode =
        MidlineConstructionMode::Auto;
    if (!resolveMidlineConstructionMode(
            candidate,
            constructionMode,
            errorMessage)) {
        return false;
    }

    if (constructionMode == MidlineConstructionMode::LegacyOneToOne) {
        return buildMidlineChainOneToOne(candidate, midlineData, errorMessage);
    }
    if (constructionMode ==
        MidlineConstructionMode::ProjectionCorrespondence) {
        return buildMidlineChainByProjectionCorrespondence(
            candidate,
            midlineData,
            errorMessage);
    }

    errorMessage = "resolved midline construction mode is invalid";
    return false;
}

bool HealClosedNarrowGap::canUseLegacyOneToOnePath(
    const ClosedGapCandidate& candidate) const
{
    // Legacy 路径保留原有语义，只处理两侧共享两个尖点的对称拓扑。
    if (candidate.startEnd.type != GapEndType::SharedTip ||
        candidate.finishEnd.type != GapEndType::SharedTip ||
        candidate.chainA.edges.size() != candidate.chainB.edges.size() ||
        candidate.chainA.vertices.size() != candidate.chainB.vertices.size() ||
        candidate.chainA.vertices.size() != candidate.chainA.edges.size() + 1) {
        return false;
    }

    try {
        for (std::size_t vertexIndex = 1;
             vertexIndex + 1 < candidate.chainA.vertices.size();
             ++vertexIndex) {
            const TopoDS_Vertex& vertexA = candidate.chainA.vertices[vertexIndex];
            const TopoDS_Vertex& vertexB = candidate.chainB.vertices[vertexIndex];
            if (vertexA.IsNull() || vertexB.IsNull()) {
                return false;
            }

            const gp_Pnt pointA = BRep_Tool::Pnt(vertexA);
            const gp_Pnt pointB = BRep_Tool::Pnt(vertexB);
            if (!GeometryUtils::IsFinitePoint(pointA) || !GeometryUtils::IsFinitePoint(pointB) ||
                pointA.Distance(pointB) >
                    params_.maxGapWidth + params_.midlineRepairTolerance) {
                return false;
            }
        }
    }
    catch (Standard_Failure&) {
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::resolveMidlineConstructionMode(
    const ClosedGapCandidate& candidate,
    MidlineConstructionMode& mode,
    std::string& errorMessage) const
{
    switch (params_.midlineConstructionMode) {
    case MidlineConstructionMode::Auto:
        // 自动模式保留原有行为：可安全一一配对时继续使用旧路径。
        mode = canUseLegacyOneToOnePath(candidate)
            ? MidlineConstructionMode::LegacyOneToOne
            : MidlineConstructionMode::ProjectionCorrespondence;
        return true;

    case MidlineConstructionMode::LegacyOneToOne:
        if (!canUseLegacyOneToOnePath(candidate)) {
            errorMessage =
                "forced LegacyOneToOne mode requires equal source edge counts "
                "and safely paired internal vertices";
            return false;
        }
        mode = MidlineConstructionMode::LegacyOneToOne;
        return true;

    case MidlineConstructionMode::ProjectionCorrespondence:
        // 投影对应是通用路径，因此也允许显式用于对称候选。
        mode = MidlineConstructionMode::ProjectionCorrespondence;
        return true;

    default:
        errorMessage = "unknown midline construction mode";
        return false;
    }
}

bool HealClosedNarrowGap::buildMidlineChainOneToOne(
    const ClosedGapCandidate& candidate,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    midlineData = MidlineChainData();
    midlineData.constructionMode = MidlineConstructionMode::LegacyOneToOne;

    try {
        const std::size_t edgeCount = candidate.chainA.edges.size();
        std::vector<TopoDS_Vertex> midVertices;
        midVertices.reserve(edgeCount + 1);

        BRep_Builder builder;
        for (std::size_t vertexIndex = 0; vertexIndex <= edgeCount; ++vertexIndex) {
            if (vertexIndex == 0) {
                // 两侧起始尖点已经拓扑重合，直接沿用原 vertex。
                midVertices.push_back(candidate.chainA.vertices.front());
                continue;
            }
            if (vertexIndex == edgeCount) {
                // 两侧结束尖点同样沿用原 vertex。
                midVertices.push_back(candidate.chainA.vertices.back());
                continue;
            }

            const TopoDS_Vertex& vertexA = candidate.chainA.vertices[vertexIndex];
            const TopoDS_Vertex& vertexB = candidate.chainB.vertices[vertexIndex];
            if (vertexA.IsNull() || vertexB.IsNull()) {
                errorMessage = "target vertex construction failure: "
                    "candidate contains a null corresponding internal vertex";
                return false;
            }

            const gp_Pnt pointA = BRep_Tool::Pnt(vertexA);
            const gp_Pnt pointB = BRep_Tool::Pnt(vertexB);
            if (!GeometryUtils::IsFinitePoint(pointA) || !GeometryUtils::IsFinitePoint(pointB)) {
                errorMessage = "target vertex construction failure: "
                    "candidate contains a non-finite internal vertex";
                return false;
            }
            const double vertexGap = pointA.Distance(pointB);
            if (vertexGap > params_.maxGapWidth + params_.midlineRepairTolerance) {
                errorMessage = "target vertex construction failure: "
                    "corresponding internal chain vertices are too far apart";
                return false;
            }

            const gp_Pnt midpoint(
                0.5 * (pointA.X() + pointB.X()),
                0.5 * (pointA.Y() + pointB.Y()),
                0.5 * (pointA.Z() + pointB.Z()));

            // 新顶点还会替换两侧旧顶点，因此容差必须覆盖旧顶点到中点的移动距离。
            const double vertexTolerance = std::max({
                Precision::Confusion(),
                params_.midlineRepairTolerance,
                BRep_Tool::Tolerance(vertexA),
                BRep_Tool::Tolerance(vertexB),
                midpoint.Distance(pointA) + Precision::Confusion(),
                midpoint.Distance(pointB) + Precision::Confusion()
            });

            TopoDS_Vertex midVertex;
            builder.MakeVertex(midVertex, midpoint, vertexTolerance);
            if (midVertex.IsNull()) {
                errorMessage = "target vertex construction failure: "
                    "failed to create midline vertex";
                return false;
            }
            midVertices.push_back(midVertex);
        }

        midlineData.chain.vertices = midVertices;
        midlineData.replacementsForChainA.resize(edgeCount);
        midlineData.replacementsForChainB.resize(edgeCount);
        midlineData.targetVertexIndicesForChainA.resize(edgeCount + 1);
        midlineData.targetVertexIndicesForChainB.resize(edgeCount + 1);
        for (std::size_t vertexIndex = 0; vertexIndex <= edgeCount; ++vertexIndex) {
            midlineData.targetVertexIndicesForChainA[vertexIndex] = vertexIndex;
            midlineData.targetVertexIndicesForChainB[vertexIndex] = vertexIndex;
        }

        auto sampleOrientedEdge = [&](const EdgeChain& chain,
                                      std::size_t edgeIndex,
                                      int sampleCount,
                                      std::vector<gp_Pnt>& sampledPoints,
                                      const char* chainName) {
            sampledPoints.clear();
            sampledPoints.reserve(static_cast<std::size_t>(sampleCount));

            const TopoDS_Edge& edge = chain.edges[edgeIndex];
            const TopoDS_Vertex& startVertex = chain.vertices[edgeIndex];
            const TopoDS_Vertex& endVertex = chain.vertices[edgeIndex + 1];
            const double startParameter = BRep_Tool::Parameter(startVertex, edge);
            const double endParameter = BRep_Tool::Parameter(endVertex, edge);
            if (!std::isfinite(startParameter) || !std::isfinite(endParameter) ||
                std::abs(endParameter - startParameter) <= Precision::PConfusion()) {
                errorMessage = std::string(chainName) + " edge has an invalid parameter range";
                return false;
            }

            BRepAdaptor_Curve curve(edge);
            const double lowerParameter = std::min(startParameter, endParameter);
            const double upperParameter = std::max(startParameter, endParameter);
            GCPnts_UniformAbscissa uniformSampler(
                curve,
                sampleCount,
                lowerParameter,
                upperParameter,
                Precision::Confusion());
            const bool useUniformAbscissa =
                uniformSampler.IsDone() && uniformSampler.NbPoints() == sampleCount;

            for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
                double parameter = 0.0;
                if (useUniformAbscissa) {
                    const int uniformIndex = startParameter <= endParameter
                        ? sampleIndex + 1
                        : sampleCount - sampleIndex;
                    parameter = uniformSampler.Parameter(uniformIndex);
                }
                else {
                    const double normalized = static_cast<double>(sampleIndex) /
                        static_cast<double>(sampleCount - 1);
                    parameter = startParameter + (endParameter - startParameter) * normalized;
                }
                sampledPoints.push_back(curve.Value(parameter));
                if (!GeometryUtils::IsFinitePoint(sampledPoints.back())) {
                    errorMessage = std::string(chainName)
                        + " edge produced a non-finite sample point";
                    return false;
                }
            }
            return true;
        };

        for (std::size_t edgeIndex = 0; edgeIndex < edgeCount; ++edgeIndex) {
            const double lengthA = GeometryUtils::ComputeEdgeLength(candidate.chainA.edges[edgeIndex]);
            const double lengthB = GeometryUtils::ComputeEdgeLength(candidate.chainB.edges[edgeIndex]);
            if (!std::isfinite(lengthA) || !std::isfinite(lengthB) ||
                lengthA <= Precision::Confusion() ||
                lengthB <= Precision::Confusion()) {
                errorMessage = "target edge construction failure: "
                    "candidate contains an invalid source edge pair";
                return false;
            }

            const int sampleCount = computeMidlineSampleCount(lengthA, lengthB);
            std::vector<gp_Pnt> sampledA;
            std::vector<gp_Pnt> sampledB;
            if (!sampleOrientedEdge(
                    candidate.chainA, edgeIndex, sampleCount, sampledA, "chainA") ||
                !sampleOrientedEdge(
                    candidate.chainB, edgeIndex, sampleCount, sampledB, "chainB")) {
                return false;
            }

            std::vector<gp_Pnt> guidancePoints;
            guidancePoints.reserve(static_cast<std::size_t>(sampleCount));
            const gp_Pnt startPoint = BRep_Tool::Pnt(midVertices[edgeIndex]);
            const gp_Pnt endPoint = BRep_Tool::Pnt(midVertices[edgeIndex + 1]);
            guidancePoints.push_back(startPoint);
            for (int sampleIndex = 1; sampleIndex + 1 < sampleCount; ++sampleIndex) {
                const gp_Pnt midpoint(
                    0.5 * (sampledA[static_cast<std::size_t>(sampleIndex)].X()
                         + sampledB[static_cast<std::size_t>(sampleIndex)].X()),
                    0.5 * (sampledA[static_cast<std::size_t>(sampleIndex)].Y()
                         + sampledB[static_cast<std::size_t>(sampleIndex)].Y()),
                    0.5 * (sampledA[static_cast<std::size_t>(sampleIndex)].Z()
                         + sampledB[static_cast<std::size_t>(sampleIndex)].Z()));
                if (guidancePoints.back().Distance(midpoint) > Precision::Confusion() &&
                    midpoint.Distance(endPoint) > Precision::Confusion()) {
                    guidancePoints.push_back(midpoint);
                }
            }
            if (guidancePoints.back().Distance(endPoint) <= Precision::Confusion()) {
                errorMessage = "target edge construction failure: "
                    "midline curve guidance points collapse at the segment endpoint";
                return false;
            }
            guidancePoints.push_back(endPoint);

            Handle(TColgp_HArray1OfPnt) interpolationPoints =
                new TColgp_HArray1OfPnt(1, static_cast<Standard_Integer>(guidancePoints.size()));
            for (std::size_t pointIndex = 0; pointIndex < guidancePoints.size(); ++pointIndex) {
                interpolationPoints->SetValue(
                    static_cast<Standard_Integer>(pointIndex + 1),
                    guidancePoints[pointIndex]);
            }

            // 采样中点只约束曲线形状，不创建额外拓扑顶点或 edge。
            GeomAPI_Interpolate interpolator(
                interpolationPoints,
                Standard_False,
                Precision::Confusion());
            interpolator.Perform();
            if (!interpolator.IsDone() || interpolator.Curve().IsNull()) {
                errorMessage = "target edge construction failure: "
                    "failed to interpolate a midline curve from sampled midpoints";
                return false;
            }

            Handle(Geom_Curve) midlineCurve = interpolator.Curve();
            BRepBuilderAPI_MakeEdge edgeMaker(
                midlineCurve,
                midVertices[edgeIndex],
                midVertices[edgeIndex + 1]);
            if (!edgeMaker.IsDone()) {
                errorMessage = "target edge construction failure: "
                    "failed to create midline edge";
                return false;
            }

            const TopoDS_Edge edge = edgeMaker.Edge();
            if (edge.IsNull()) {
                errorMessage = "target edge construction failure: "
                    "created midline edge is null";
                return false;
            }
            const double midlineEdgeLength = GeometryUtils::ComputeEdgeLength(edge);
            if (!std::isfinite(midlineEdgeLength) ||
                midlineEdgeLength <= Precision::Confusion()) {
                errorMessage = "target edge construction failure: "
                    "created midline edge has invalid length";
                return false;
            }
            midlineData.chain.edges.push_back(edge);
            midlineData.chain.length += midlineEdgeLength;
            const std::size_t midlineEdgeIndex = midlineData.chain.edges.size() - 1;

            // source edge 与中线 edge 一一对应，采样密度不再改变拓扑分段数量。
            midlineData.replacementsForChainA[edgeIndex].push_back(midlineEdgeIndex);
            midlineData.replacementsForChainB[edgeIndex].push_back(midlineEdgeIndex);
        }

        if (midlineData.chain.edges.empty() ||
            !std::isfinite(midlineData.chain.length) ||
            midlineData.chain.length <= Precision::Confusion()) {
            errorMessage = "midline chain is empty";
            return false;
        }

        if (!validateReplacementMapping(
                midlineData.replacementsForChainA,
                candidate.chainA.edges.size(),
                midlineData.chain.edges.size(),
                "chainA",
                errorMessage) ||
            !validateReplacementMapping(
                midlineData.replacementsForChainB,
                candidate.chainB.edges.size(),
                midlineData.chain.edges.size(),
                "chainB",
                errorMessage)) {
            return false;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "legacy one-to-one midline construction failed: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::buildNormalizedVertexStations(
    const EdgeChain& chain,
    const char* chainName,
    ChainStationData& stationData,
    std::string& errorMessage) const
{
    stationData = ChainStationData();
    const std::string name = chainName != nullptr ? chainName : "chain";
    if (chain.edges.empty() ||
        chain.vertices.size() != chain.edges.size() + 1 ||
        chain.ownerFaces.size() != chain.edges.size()) {
        errorMessage = "station construction failure: " + name
            + " topology counts are inconsistent";
        return false;
    }
    if (!std::isfinite(chain.length) ||
        chain.length <= Precision::Confusion()) {
        errorMessage = "station construction failure: " + name
            + " total length is invalid";
        return false;
    }

    stationData.edgeLengths.reserve(chain.edges.size());
    for (const TopoDS_Edge& edge : chain.edges) {
        const double edgeLength = GeometryUtils::ComputeEdgeLength(edge);
        if (!std::isfinite(edgeLength) ||
            edgeLength <= Precision::Confusion()) {
            errorMessage = "station construction failure: " + name
                + " contains an invalid source edge length";
            return false;
        }

        stationData.edgeLengths.push_back(edgeLength);
        stationData.totalLength += edgeLength;
        if (!std::isfinite(stationData.totalLength)) {
            errorMessage = "station construction failure: " + name
                + " accumulated length is not finite";
            return false;
        }
    }

    if (stationData.totalLength <= Precision::Confusion()) {
        errorMessage = "station construction failure: " + name
            + " computed total length is invalid";
        return false;
    }

    stationData.vertexStations.reserve(chain.vertices.size());
    stationData.vertexStations.push_back(0.0);
    double accumulatedLength = 0.0;
    for (std::size_t edgeIndex = 0;
         edgeIndex + 1 < stationData.edgeLengths.size();
         ++edgeIndex) {
        accumulatedLength += stationData.edgeLengths[edgeIndex];
        const double station = accumulatedLength / stationData.totalLength;
        if (!std::isfinite(station) ||
            station <= stationData.vertexStations.back() ||
            station >= 1.0) {
            errorMessage = "station construction failure: " + name
                + " vertex stations are not strictly increasing";
            return false;
        }
        stationData.vertexStations.push_back(station);
    }
    stationData.vertexStations.push_back(1.0);

    if (stationData.vertexStations.size() != chain.vertices.size()) {
        errorMessage = "station construction failure: " + name
            + " vertex station count is inconsistent";
        return false;
    }
    return true;
}

bool HealClosedNarrowGap::projectPointToChain(
    const gp_Pnt& sourcePoint,
    const EdgeChain& targetChain,
    const ChainStationData& targetStationData,
    ChainProjectionResult& projection,
    std::string& errorMessage) const
{
    projection = ChainProjectionResult();
    const std::size_t invalidIndex =
        std::numeric_limits<std::size_t>::max();
    if (!GeometryUtils::IsFinitePoint(sourcePoint) ||
        targetChain.edges.empty() ||
        targetChain.vertices.size() != targetChain.edges.size() + 1 ||
        targetStationData.edgeLengths.size() != targetChain.edges.size() ||
        targetStationData.vertexStations.size() != targetChain.vertices.size() ||
        !std::isfinite(targetStationData.totalLength) ||
        targetStationData.totalLength <= Precision::Confusion()) {
        errorMessage = "chain projection failure: target chain data is invalid";
        return false;
    }

    std::size_t bestEdgeIndex = invalidIndex;
    double bestParameter = 0.0;
    gp_Pnt bestPoint;
    double bestDistance = std::numeric_limits<double>::infinity();
    std::string nearestError;
    if (!findNearestPointOnChain(
            sourcePoint,
            targetChain,
            bestEdgeIndex,
            bestParameter,
            bestPoint,
            bestDistance,
            nearestError)) {
        errorMessage = "chain projection failure: " + nearestError;
        return false;
    }

    try {
        const double allowedDistance =
            params_.maxGapWidth + params_.midlineRepairTolerance;
        if (bestDistance > allowedDistance) {
            errorMessage = "chain projection failure: nearest point exceeds allowed gap width";
            return false;
        }

        const TopoDS_Edge& bestEdge =
            targetChain.edges[bestEdgeIndex];
        BRepAdaptor_Curve bestCurve(bestEdge);
        const double edgeStartParameter = BRep_Tool::Parameter(
            targetChain.vertices[bestEdgeIndex],
            bestEdge);
        double localLength = GCPnts_AbscissaPoint::Length(
            bestCurve,
            edgeStartParameter,
            bestParameter,
            Precision::Confusion());
        if (!std::isfinite(localLength) ||
            localLength < -Precision::Confusion() ||
            localLength >
                targetStationData.edgeLengths[bestEdgeIndex]
                    + Precision::Confusion()) {
            errorMessage = "chain projection failure: projected arc length is invalid";
            return false;
        }
        localLength = std::max(
            0.0,
            std::min(
                targetStationData.edgeLengths[bestEdgeIndex],
                localLength));

        double accumulatedLength = 0.0;
        for (std::size_t edgeIndex = 0;
             edgeIndex < bestEdgeIndex;
             ++edgeIndex) {
            accumulatedLength += targetStationData.edgeLengths[edgeIndex];
        }
        const double targetStation = GeometryUtils::Clamp01(
            (accumulatedLength + localLength)
            / targetStationData.totalLength);
        if (!std::isfinite(targetStation)) {
            errorMessage = "chain projection failure: projected chain station is invalid";
            return false;
        }

        projection.projectedPoint = bestPoint;
        projection.rawProjectedPoint = bestPoint;
        projection.targetEdgeIndex = bestEdgeIndex;
        projection.targetParameter = bestParameter;
        projection.targetStation = targetStation;
        projection.rawTargetStation = targetStation;
        projection.distance = bestDistance;
        projection.rawDistance = bestDistance;

        const double snapTolerance = std::min(
            params_.projectionVertexSnapTolerance,
            params_.maxGapWidth);
        const double stationSnapTolerance = std::min(
            0.25,
            snapTolerance / std::max(
                targetStationData.totalLength,
                Precision::Confusion()));

        std::size_t snappedVertexIndex = invalidIndex;
        double snappedVertexDistance =
            std::numeric_limits<double>::infinity();
        const std::size_t endpointIndices[2] = {
            bestEdgeIndex,
            bestEdgeIndex + 1
        };
        for (const std::size_t vertexIndex : endpointIndices) {
            const gp_Pnt vertexPoint =
                BRep_Tool::Pnt(targetChain.vertices[vertexIndex]);
            const double geometricDistance =
                bestPoint.Distance(vertexPoint);
            const double stationDistance = std::abs(
                targetStation
                - targetStationData.vertexStations[vertexIndex]);
            if (!std::isfinite(geometricDistance) ||
                geometricDistance > snapTolerance ||
                stationDistance > stationSnapTolerance) {
                continue;
            }

            if (snappedVertexIndex == invalidIndex ||
                geometricDistance + Precision::Confusion()
                    < snappedVertexDistance ||
                (std::abs(geometricDistance - snappedVertexDistance)
                        <= Precision::Confusion() &&
                    vertexIndex < snappedVertexIndex)) {
                snappedVertexIndex = vertexIndex;
                snappedVertexDistance = geometricDistance;
            }
        }

        // 吸附容差独立于细缝宽度，并同时检查链站点，避免吸到空间近但拓扑远的顶点。
        if (snappedVertexIndex != invalidIndex) {
            projection.projectedPoint =
                BRep_Tool::Pnt(targetChain.vertices[snappedVertexIndex]);
            projection.targetStation =
                targetStationData.vertexStations[snappedVertexIndex];
            projection.distance =
                sourcePoint.Distance(projection.projectedPoint);
            projection.snappedTargetVertexIndex =
                snappedVertexIndex;
            projection.snappedToVertex = true;
            if (!std::isfinite(projection.distance) ||
                projection.distance > allowedDistance) {
                errorMessage = "chain projection failure: snapped vertex exceeds allowed gap width";
                return false;
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "chain projection failure: " + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::selectSnappedVertexPairs(
    const EdgeChain& chainA,
    const EdgeChain& chainB,
    const std::vector<ChainProjectionResult>& projectionsFromA,
    const std::vector<ChainProjectionResult>& projectionsFromB,
    std::vector<SnappedVertexPair>& selectedPairs,
    std::string& errorMessage) const
{
    selectedPairs.clear();
    if (projectionsFromA.size() != chainA.vertices.size() ||
        projectionsFromB.size() != chainB.vertices.size()) {
        errorMessage = "projection conflict failure: projection counts are inconsistent";
        return false;
    }

    struct PairProposal {
        std::size_t vertexIndexA = 0;
        std::size_t vertexIndexB = 0;
        double distance = 0.0;
    };

    std::vector<PairProposal> proposals;
    auto addProposal = [&](std::size_t vertexIndexA,
                           std::size_t vertexIndexB,
                           double distance) {
        if (vertexIndexA == 0 ||
            vertexIndexA + 1 >= chainA.vertices.size() ||
            vertexIndexB == 0 ||
            vertexIndexB + 1 >= chainB.vertices.size() ||
            !std::isfinite(distance)) {
            return;
        }

        for (PairProposal& proposal : proposals) {
            if (proposal.vertexIndexA == vertexIndexA &&
                proposal.vertexIndexB == vertexIndexB) {
                proposal.distance = std::min(proposal.distance, distance);
                return;
            }
        }

        PairProposal proposal;
        proposal.vertexIndexA = vertexIndexA;
        proposal.vertexIndexB = vertexIndexB;
        proposal.distance = distance;
        proposals.push_back(proposal);
    };

    for (std::size_t vertexIndexA = 1;
         vertexIndexA + 1 < chainA.vertices.size();
         ++vertexIndexA) {
        const ChainProjectionResult& projection =
            projectionsFromA[vertexIndexA];
        if (projection.snappedToVertex) {
            addProposal(
                vertexIndexA,
                projection.snappedTargetVertexIndex,
                projection.distance);
        }
    }
    for (std::size_t vertexIndexB = 1;
         vertexIndexB + 1 < chainB.vertices.size();
         ++vertexIndexB) {
        const ChainProjectionResult& projection =
            projectionsFromB[vertexIndexB];
        if (projection.snappedToVertex) {
            addProposal(
                projection.snappedTargetVertexIndex,
                vertexIndexB,
                projection.distance);
        }
    }

    std::sort(
        proposals.begin(),
        proposals.end(),
        [](const PairProposal& lhs, const PairProposal& rhs) {
            if (lhs.distance != rhs.distance) {
                return lhs.distance < rhs.distance;
            }
            if (lhs.vertexIndexA != rhs.vertexIndexA) {
                return lhs.vertexIndexA < rhs.vertexIndexA;
            }
            return lhs.vertexIndexB < rhs.vertexIndexB;
        });

    std::vector<bool> pairedVerticesA(chainA.vertices.size(), false);
    std::vector<bool> pairedVerticesB(chainB.vertices.size(), false);
    for (const PairProposal& proposal : proposals) {
        if (pairedVerticesA[proposal.vertexIndexA] ||
            pairedVerticesB[proposal.vertexIndexB]) {
            continue;
        }

        // 一个已有顶点只能配对一次；按距离和顶点下标排序保证冲突处理确定性。
        pairedVerticesA[proposal.vertexIndexA] = true;
        pairedVerticesB[proposal.vertexIndexB] = true;
        SnappedVertexPair selectedPair;
        selectedPair.sourceVertexIndexA = proposal.vertexIndexA;
        selectedPair.sourceVertexIndexB = proposal.vertexIndexB;
        selectedPair.distance = proposal.distance;
        selectedPairs.push_back(selectedPair);
    }

    std::sort(
        selectedPairs.begin(),
        selectedPairs.end(),
        [](const SnappedVertexPair& lhs,
           const SnappedVertexPair& rhs) {
            if (lhs.sourceVertexIndexA != rhs.sourceVertexIndexA) {
                return lhs.sourceVertexIndexA < rhs.sourceVertexIndexA;
            }
            return lhs.sourceVertexIndexB < rhs.sourceVertexIndexB;
        });
    return true;
}

bool HealClosedNarrowGap::buildProjectionCorrespondenceNodes(
    const ClosedGapCandidate& candidate,
    const ChainStationData& chainAStations,
    const ChainStationData& chainBStations,
    std::vector<ProjectionCorrespondenceNode>& nodes,
    std::string& errorMessage) const
{
    nodes.clear();
    if (candidate.chainA.vertices.size() !=
            chainAStations.vertexStations.size() ||
        candidate.chainB.vertices.size() !=
            chainBStations.vertexStations.size()) {
        errorMessage = "correspondence construction failure: source station counts are inconsistent";
        return false;
    }

    try {
        for (const TopoDS_Vertex& vertex : candidate.chainA.vertices) {
            if (vertex.IsNull() || !GeometryUtils::IsFinitePoint(BRep_Tool::Pnt(vertex))) {
                errorMessage = "correspondence construction failure: chainA contains an invalid vertex";
                return false;
            }
        }
        for (const TopoDS_Vertex& vertex : candidate.chainB.vertices) {
            if (vertex.IsNull() || !GeometryUtils::IsFinitePoint(BRep_Tool::Pnt(vertex))) {
                errorMessage = "correspondence construction failure: chainB contains an invalid vertex";
                return false;
            }
        }

        std::vector<ChainProjectionResult> projectionsFromA(
            candidate.chainA.vertices.size());
        std::vector<ChainProjectionResult> projectionsFromB(
            candidate.chainB.vertices.size());

        // 必须双向投影，否则只存在于另一侧的内部拓扑顶点会被遗漏。
        for (std::size_t vertexIndexA = 1;
             vertexIndexA + 1 < candidate.chainA.vertices.size();
             ++vertexIndexA) {
            if (!projectPointToChain(
                    BRep_Tool::Pnt(
                        candidate.chainA.vertices[vertexIndexA]),
                    candidate.chainB,
                    chainBStations,
                    projectionsFromA[vertexIndexA],
                    errorMessage)) {
                errorMessage = "chainA vertex projection failed: " + errorMessage;
                return false;
            }
        }
        for (std::size_t vertexIndexB = 1;
             vertexIndexB + 1 < candidate.chainB.vertices.size();
             ++vertexIndexB) {
            if (!projectPointToChain(
                    BRep_Tool::Pnt(
                        candidate.chainB.vertices[vertexIndexB]),
                    candidate.chainA,
                    chainAStations,
                    projectionsFromB[vertexIndexB],
                    errorMessage)) {
                errorMessage = "chainB vertex projection failed: " + errorMessage;
                return false;
            }
        }

        std::vector<SnappedVertexPair> selectedPairs;
        if (!selectSnappedVertexPairs(
                candidate.chainA,
                candidate.chainB,
                projectionsFromA,
                projectionsFromB,
                selectedPairs,
                errorMessage)) {
            return false;
        }

        std::vector<bool> pairedVerticesA(
            candidate.chainA.vertices.size(),
            false);
        std::vector<bool> pairedVerticesB(
            candidate.chainB.vertices.size(),
            false);

        ProjectionCorrespondenceNode startNode;
        startNode.stationA = 0.0;
        startNode.stationB = 0.0;
        startNode.pointA =
            BRep_Tool::Pnt(candidate.chainA.vertices.front());
        startNode.pointB =
            BRep_Tool::Pnt(candidate.chainB.vertices.front());
        startNode.endpoint = true;
        nodes.push_back(startNode);

        for (const SnappedVertexPair& pair : selectedPairs) {
            if (pair.sourceVertexIndexA == 0 ||
                pair.sourceVertexIndexA + 1 >=
                    candidate.chainA.vertices.size() ||
                pair.sourceVertexIndexB == 0 ||
                pair.sourceVertexIndexB + 1 >=
                    candidate.chainB.vertices.size()) {
                errorMessage = "correspondence construction failure: snapped pair index is invalid";
                return false;
            }

            pairedVerticesA[pair.sourceVertexIndexA] = true;
            pairedVerticesB[pair.sourceVertexIndexB] = true;

            ProjectionCorrespondenceNode node;
            node.stationA =
                chainAStations.vertexStations[
                    pair.sourceVertexIndexA];
            node.stationB =
                chainBStations.vertexStations[
                    pair.sourceVertexIndexB];
            node.pointA = BRep_Tool::Pnt(
                candidate.chainA.vertices[
                    pair.sourceVertexIndexA]);
            node.pointB = BRep_Tool::Pnt(
                candidate.chainB.vertices[
                    pair.sourceVertexIndexB]);
            node.sourceVertexIndexA =
                pair.sourceVertexIndexA;
            node.sourceVertexIndexB =
                pair.sourceVertexIndexB;
            node.projectionDistance = pair.distance;
            node.hasSourceVertexA = true;
            node.hasSourceVertexB = true;
            node.snappedPair = true;
            nodes.push_back(node);
        }

        for (std::size_t vertexIndexA = 1;
             vertexIndexA + 1 < candidate.chainA.vertices.size();
             ++vertexIndexA) {
            if (pairedVerticesA[vertexIndexA]) {
                continue;
            }

            const ChainProjectionResult& projection =
                projectionsFromA[vertexIndexA];
            ProjectionCorrespondenceNode node;
            node.stationA =
                chainAStations.vertexStations[vertexIndexA];
            node.stationB = projection.rawTargetStation;
            node.pointA = BRep_Tool::Pnt(
                candidate.chainA.vertices[vertexIndexA]);
            node.pointB = projection.rawProjectedPoint;
            node.sourceVertexIndexA = vertexIndexA;
            node.projectionDistance = projection.rawDistance;
            node.hasSourceVertexA = true;
            nodes.push_back(node);
        }

        for (std::size_t vertexIndexB = 1;
             vertexIndexB + 1 < candidate.chainB.vertices.size();
             ++vertexIndexB) {
            if (pairedVerticesB[vertexIndexB]) {
                continue;
            }

            const ChainProjectionResult& projection =
                projectionsFromB[vertexIndexB];
            ProjectionCorrespondenceNode node;
            node.stationA = projection.rawTargetStation;
            node.stationB =
                chainBStations.vertexStations[vertexIndexB];
            node.pointA = projection.rawProjectedPoint;
            node.pointB = BRep_Tool::Pnt(
                candidate.chainB.vertices[vertexIndexB]);
            node.sourceVertexIndexB = vertexIndexB;
            node.projectionDistance = projection.rawDistance;
            node.hasSourceVertexB = true;
            nodes.push_back(node);
        }

        ProjectionCorrespondenceNode endNode;
        endNode.stationA = 1.0;
        endNode.stationB = 1.0;
        endNode.pointA =
            BRep_Tool::Pnt(candidate.chainA.vertices.back());
        endNode.pointB =
            BRep_Tool::Pnt(candidate.chainB.vertices.back());
        endNode.endpoint = true;
        nodes.push_back(endNode);
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "correspondence construction failure: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return mergeProjectionCorrespondenceNodes(
        nodes,
        candidate.chainA,
        candidate.chainB,
        chainAStations,
        chainBStations,
        errorMessage);
}

bool HealClosedNarrowGap::mergeProjectionCorrespondenceNodes(
    std::vector<ProjectionCorrespondenceNode>& nodes,
    const EdgeChain& chainA,
    const EdgeChain& chainB,
    const ChainStationData& chainAStations,
    const ChainStationData& chainBStations,
    std::string& errorMessage) const
{
    if (nodes.size() < 2) {
        errorMessage = "correspondence merge failure: node sequence is incomplete";
        return false;
    }

    const double snapTolerance = std::min(
        params_.projectionVertexSnapTolerance,
        params_.maxGapWidth);
    const double stationMergeToleranceA = std::min(
        0.05,
        std::max(
            kStationMergeTolerance,
            snapTolerance / std::max(
                chainAStations.totalLength,
                Precision::Confusion())));
    const double stationMergeToleranceB = std::min(
        0.05,
        std::max(
            kStationMergeTolerance,
            snapTolerance / std::max(
                chainBStations.totalLength,
                Precision::Confusion())));

    auto midpointOfNode =
        [](const ProjectionCorrespondenceNode& node) {
            return gp_Pnt(
                0.5 * (node.pointA.X() + node.pointB.X()),
                0.5 * (node.pointA.Y() + node.pointB.Y()),
                0.5 * (node.pointA.Z() + node.pointB.Z()));
        };

    std::sort(
        nodes.begin(),
        nodes.end(),
        [](const ProjectionCorrespondenceNode& lhs,
           const ProjectionCorrespondenceNode& rhs) {
            if (lhs.stationA != rhs.stationA) {
                return lhs.stationA < rhs.stationA;
            }
            if (lhs.stationB != rhs.stationB) {
                return lhs.stationB < rhs.stationB;
            }
            if (lhs.sourceVertexIndexA != rhs.sourceVertexIndexA) {
                return lhs.sourceVertexIndexA < rhs.sourceVertexIndexA;
            }
            return lhs.sourceVertexIndexB < rhs.sourceVertexIndexB;
        });

    std::vector<ProjectionCorrespondenceNode> mergedNodes;
    for (const ProjectionCorrespondenceNode& incomingNode : nodes) {
        if (!std::isfinite(incomingNode.stationA) ||
            !std::isfinite(incomingNode.stationB) ||
            incomingNode.stationA < 0.0 ||
            incomingNode.stationA > 1.0 ||
            incomingNode.stationB < 0.0 ||
            incomingNode.stationB > 1.0 ||
            !GeometryUtils::IsFinitePoint(incomingNode.pointA) ||
            !GeometryUtils::IsFinitePoint(incomingNode.pointB)) {
            errorMessage = "correspondence merge failure: node contains invalid geometry";
            return false;
        }

        std::size_t mergeIndex =
            std::numeric_limits<std::size_t>::max();
        for (std::size_t existingIndex = mergedNodes.size();
             existingIndex > 0;
             --existingIndex) {
            const ProjectionCorrespondenceNode& existingNode =
                mergedNodes[existingIndex - 1];
            if (incomingNode.stationA - existingNode.stationA >
                stationMergeToleranceA) {
                break;
            }
            if (std::abs(incomingNode.stationA - existingNode.stationA) >
                    stationMergeToleranceA ||
                std::abs(incomingNode.stationB - existingNode.stationB) >
                    stationMergeToleranceB ||
                midpointOfNode(incomingNode).Distance(
                    midpointOfNode(existingNode)) > snapTolerance) {
                continue;
            }
            if (incomingNode.endpoint != existingNode.endpoint) {
                continue;
            }
            mergeIndex = existingIndex - 1;
            break;
        }

        if (mergeIndex ==
            std::numeric_limits<std::size_t>::max()) {
            mergedNodes.push_back(incomingNode);
            continue;
        }

        ProjectionCorrespondenceNode& existingNode =
            mergedNodes[mergeIndex];
        if (existingNode.hasSourceVertexA &&
            incomingNode.hasSourceVertexA &&
            existingNode.sourceVertexIndexA !=
                incomingNode.sourceVertexIndexA) {
            errorMessage = "correspondence merge failure: different chainA vertices collapse to one node";
            return false;
        }
        if (existingNode.hasSourceVertexB &&
            incomingNode.hasSourceVertexB &&
            existingNode.sourceVertexIndexB !=
                incomingNode.sourceVertexIndexB) {
            errorMessage = "correspondence merge failure: different chainB vertices collapse to one node";
            return false;
        }

        const bool keepExistingGeometry =
            existingNode.snappedPair ||
            (!incomingNode.snappedPair &&
                existingNode.projectionDistance <=
                    incomingNode.projectionDistance);
        if (incomingNode.hasSourceVertexA &&
            !existingNode.hasSourceVertexA) {
            existingNode.stationA = incomingNode.stationA;
            existingNode.pointA = incomingNode.pointA;
            existingNode.sourceVertexIndexA =
                incomingNode.sourceVertexIndexA;
            existingNode.hasSourceVertexA = true;
        }
        else if (!existingNode.hasSourceVertexA &&
                 !incomingNode.hasSourceVertexA &&
                 !keepExistingGeometry) {
            existingNode.stationA = incomingNode.stationA;
            existingNode.pointA = incomingNode.pointA;
        }

        if (incomingNode.hasSourceVertexB &&
            !existingNode.hasSourceVertexB) {
            existingNode.stationB = incomingNode.stationB;
            existingNode.pointB = incomingNode.pointB;
            existingNode.sourceVertexIndexB =
                incomingNode.sourceVertexIndexB;
            existingNode.hasSourceVertexB = true;
        }
        else if (!existingNode.hasSourceVertexB &&
                 !incomingNode.hasSourceVertexB &&
                 !keepExistingGeometry) {
            existingNode.stationB = incomingNode.stationB;
            existingNode.pointB = incomingNode.pointB;
        }

        existingNode.projectionDistance = std::min(
            existingNode.projectionDistance,
            incomingNode.projectionDistance);
        existingNode.snappedPair =
            existingNode.snappedPair || incomingNode.snappedPair;
    }

    std::sort(
        mergedNodes.begin(),
        mergedNodes.end(),
        [](const ProjectionCorrespondenceNode& lhs,
           const ProjectionCorrespondenceNode& rhs) {
            if (lhs.stationA != rhs.stationA) {
                return lhs.stationA < rhs.stationA;
            }
            return lhs.stationB < rhs.stationB;
        });

    if (mergedNodes.size() < 2 ||
        !mergedNodes.front().endpoint ||
        !mergedNodes.back().endpoint ||
        std::abs(mergedNodes.front().stationA) >
            kStationMergeTolerance ||
        std::abs(mergedNodes.front().stationB) >
            kStationMergeTolerance ||
        std::abs(mergedNodes.back().stationA - 1.0) >
            kStationMergeTolerance ||
        std::abs(mergedNodes.back().stationB - 1.0) >
            kStationMergeTolerance) {
        errorMessage = "correspondence validation failure: endpoint nodes are invalid";
        return false;
    }
    mergedNodes.front().stationA = 0.0;
    mergedNodes.front().stationB = 0.0;
    mergedNodes.back().stationA = 1.0;
    mergedNodes.back().stationB = 1.0;

    // 对应关系必须在两条 chain 上同时前进；投影回退不能靠重新排序掩盖。
    for (std::size_t nodeIndex = 1;
         nodeIndex < mergedNodes.size();
         ++nodeIndex) {
        const ProjectionCorrespondenceNode& previous =
            mergedNodes[nodeIndex - 1];
        const ProjectionCorrespondenceNode& current =
            mergedNodes[nodeIndex];
        if (current.stationA - previous.stationA <=
                stationMergeToleranceA ||
            current.stationB - previous.stationB <=
                stationMergeToleranceB) {
            errorMessage = "correspondence validation failure: nodes are not strictly monotonic on both chains";
            return false;
        }
    }

    std::vector<int> sourceVertexUseCountA(
        chainA.vertices.size(),
        0);
    std::vector<int> sourceVertexUseCountB(
        chainB.vertices.size(),
        0);
    for (const ProjectionCorrespondenceNode& node : mergedNodes) {
        if (node.pointA.Distance(node.pointB) >
            params_.maxGapWidth + params_.midlineRepairTolerance) {
            errorMessage = "correspondence validation failure: node exceeds allowed gap width";
            return false;
        }
        if (node.hasSourceVertexA) {
            if (node.sourceVertexIndexA == 0 ||
                node.sourceVertexIndexA + 1 >=
                    chainA.vertices.size()) {
                errorMessage = "correspondence validation failure: chainA source vertex index is invalid";
                return false;
            }
            ++sourceVertexUseCountA[node.sourceVertexIndexA];
        }
        if (node.hasSourceVertexB) {
            if (node.sourceVertexIndexB == 0 ||
                node.sourceVertexIndexB + 1 >=
                    chainB.vertices.size()) {
                errorMessage = "correspondence validation failure: chainB source vertex index is invalid";
                return false;
            }
            ++sourceVertexUseCountB[node.sourceVertexIndexB];
        }
    }

    for (std::size_t vertexIndexA = 1;
         vertexIndexA + 1 < sourceVertexUseCountA.size();
         ++vertexIndexA) {
        if (sourceVertexUseCountA[vertexIndexA] != 1) {
            errorMessage = "correspondence validation failure: a chainA internal vertex is not mapped exactly once";
            return false;
        }
    }
    for (std::size_t vertexIndexB = 1;
         vertexIndexB + 1 < sourceVertexUseCountB.size();
         ++vertexIndexB) {
        if (sourceVertexUseCountB[vertexIndexB] != 1) {
            errorMessage = "correspondence validation failure: a chainB internal vertex is not mapped exactly once";
            return false;
        }
    }

    nodes = mergedNodes;
    return true;
}

bool HealClosedNarrowGap::buildProjectionTargetVertices(
    const ClosedGapCandidate& candidate,
    const std::vector<ProjectionCorrespondenceNode>& nodes,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    if (nodes.size() < 2) {
        errorMessage = "target vertex construction failure: correspondence nodes are incomplete";
        return false;
    }

    const std::size_t invalidIndex =
        std::numeric_limits<std::size_t>::max();
    midlineData.chain.vertices.clear();
    midlineData.chain.vertices.reserve(nodes.size());
    midlineData.targetVertexIndicesForChainA.assign(
        candidate.chainA.vertices.size(),
        invalidIndex);
    midlineData.targetVertexIndicesForChainB.assign(
        candidate.chainB.vertices.size(),
        invalidIndex);

    try {
        BRep_Builder builder;
        for (std::size_t nodeIndex = 0;
             nodeIndex < nodes.size();
             ++nodeIndex) {
            if (nodeIndex == 0) {
                TopoDS_Vertex targetStart;
                if (!buildTargetEndVertex(
                        candidate.startEnd,
                        targetStart,
                        errorMessage)) {
                    return false;
                }
                midlineData.chain.vertices.push_back(targetStart);
                midlineData.targetVertexIndicesForChainA.front() = 0;
                midlineData.targetVertexIndicesForChainB.front() = 0;
                continue;
            }
            if (nodeIndex + 1 == nodes.size()) {
                TopoDS_Vertex targetFinish;
                if (!buildTargetEndVertex(
                        candidate.finishEnd,
                        targetFinish,
                        errorMessage)) {
                    return false;
                }
                midlineData.chain.vertices.push_back(targetFinish);
                midlineData.targetVertexIndicesForChainA.back() =
                    nodeIndex;
                midlineData.targetVertexIndicesForChainB.back() =
                    nodeIndex;
                continue;
            }

            const ProjectionCorrespondenceNode& node =
                nodes[nodeIndex];
            // snapped pair 使用两个真实顶点；普通节点使用源顶点和对侧投影点，
            // 两种情况都取二者中点，避免把整条缝隙强行压到任意一侧。
            const gp_Pnt targetPoint(
                0.5 * (node.pointA.X() + node.pointB.X()),
                0.5 * (node.pointA.Y() + node.pointB.Y()),
                0.5 * (node.pointA.Z() + node.pointB.Z()));
            if (!GeometryUtils::IsFinitePoint(targetPoint)) {
                errorMessage = "target vertex construction failure: midpoint is not finite";
                return false;
            }

            double vertexTolerance = std::max({
                Precision::Confusion(),
                params_.midlineRepairTolerance,
                node.pointA.Distance(targetPoint)
                    + Precision::Confusion(),
                node.pointB.Distance(targetPoint)
                    + Precision::Confusion()
            });
            if (node.hasSourceVertexA) {
                const TopoDS_Vertex& sourceVertex =
                    candidate.chainA.vertices[
                        node.sourceVertexIndexA];
                vertexTolerance = std::max({
                    vertexTolerance,
                    BRep_Tool::Tolerance(sourceVertex),
                    BRep_Tool::Pnt(sourceVertex).Distance(targetPoint)
                        + Precision::Confusion()
                });
            }
            if (node.hasSourceVertexB) {
                const TopoDS_Vertex& sourceVertex =
                    candidate.chainB.vertices[
                        node.sourceVertexIndexB];
                vertexTolerance = std::max({
                    vertexTolerance,
                    BRep_Tool::Tolerance(sourceVertex),
                    BRep_Tool::Pnt(sourceVertex).Distance(targetPoint)
                        + Precision::Confusion()
                });
            }

            TopoDS_Vertex targetVertex;
            builder.MakeVertex(
                targetVertex,
                targetPoint,
                vertexTolerance);
            if (targetVertex.IsNull()) {
                errorMessage = "target vertex construction failure: failed to create midpoint vertex";
                return false;
            }
            midlineData.chain.vertices.push_back(targetVertex);

            if (node.hasSourceVertexA) {
                std::size_t& mappedIndex =
                    midlineData.targetVertexIndicesForChainA[
                        node.sourceVertexIndexA];
                if (mappedIndex != invalidIndex) {
                    errorMessage = "target vertex construction failure: chainA vertex maps more than once";
                    return false;
                }
                mappedIndex = nodeIndex;
            }
            if (node.hasSourceVertexB) {
                std::size_t& mappedIndex =
                    midlineData.targetVertexIndicesForChainB[
                        node.sourceVertexIndexB];
                if (mappedIndex != invalidIndex) {
                    errorMessage = "target vertex construction failure: chainB vertex maps more than once";
                    return false;
                }
                mappedIndex = nodeIndex;
            }
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "target vertex construction failure: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    for (std::size_t vertexIndexA = 0;
         vertexIndexA <
            midlineData.targetVertexIndicesForChainA.size();
         ++vertexIndexA) {
        if (midlineData.targetVertexIndicesForChainA[vertexIndexA] ==
            invalidIndex) {
            errorMessage = "target vertex construction failure: chainA vertex mapping is incomplete";
            return false;
        }
    }
    for (std::size_t vertexIndexB = 0;
         vertexIndexB <
            midlineData.targetVertexIndicesForChainB.size();
         ++vertexIndexB) {
        if (midlineData.targetVertexIndicesForChainB[vertexIndexB] ==
            invalidIndex) {
            errorMessage = "target vertex construction failure: chainB vertex mapping is incomplete";
            return false;
        }
    }
    if (midlineData.chain.vertices.size() != nodes.size()) {
        errorMessage = "target vertex construction failure: target vertex count is inconsistent";
        return false;
    }
    return true;
}

bool HealClosedNarrowGap::buildProjectionTargetEdges(
    const ClosedGapCandidate& candidate,
    const ChainStationData& chainAStations,
    const ChainStationData& chainBStations,
    const std::vector<ProjectionCorrespondenceNode>& nodes,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    if (nodes.size() < 2 ||
        midlineData.chain.vertices.size() != nodes.size()) {
        errorMessage = "target edge construction failure: target vertices are inconsistent";
        return false;
    }

    midlineData.chain.edges.clear();
    midlineData.chain.length = 0.0;
    midlineData.chain.edges.reserve(nodes.size() - 1);
    try {
        for (std::size_t targetEdgeIndex = 0;
             targetEdgeIndex + 1 < nodes.size();
             ++targetEdgeIndex) {
            const ProjectionCorrespondenceNode& startNode =
                nodes[targetEdgeIndex];
            const ProjectionCorrespondenceNode& endNode =
                nodes[targetEdgeIndex + 1];
            const double stationSpanA =
                endNode.stationA - startNode.stationA;
            const double stationSpanB =
                endNode.stationB - startNode.stationB;
            if (!std::isfinite(stationSpanA) ||
                !std::isfinite(stationSpanB) ||
                stationSpanA <= 0.0 ||
                stationSpanB <= 0.0) {
                errorMessage = "target edge construction failure: correspondence interval is invalid";
                return false;
            }

            const double segmentLengthA =
                chainAStations.totalLength * stationSpanA;
            const double segmentLengthB =
                chainBStations.totalLength * stationSpanB;
            const int sampleCount = computeMidlineSampleCount(
                segmentLengthA,
                segmentLengthB);
            if (sampleCount < 2) {
                errorMessage = "target edge construction failure: sample count is invalid";
                return false;
            }

            const gp_Pnt startPoint = BRep_Tool::Pnt(
                midlineData.chain.vertices[targetEdgeIndex]);
            const gp_Pnt endPoint = BRep_Tool::Pnt(
                midlineData.chain.vertices[targetEdgeIndex + 1]);
            if (!GeometryUtils::IsFinitePoint(startPoint) ||
                !GeometryUtils::IsFinitePoint(endPoint) ||
                startPoint.Distance(endPoint) <=
                    Precision::Confusion()) {
                errorMessage = "target edge construction failure: target endpoints are invalid";
                return false;
            }

            std::vector<gp_Pnt> guidancePoints;
            guidancePoints.reserve(
                static_cast<std::size_t>(sampleCount));
            guidancePoints.push_back(startPoint);
            for (int sampleIndex = 1;
                 sampleIndex + 1 < sampleCount;
                 ++sampleIndex) {
                const double local =
                    static_cast<double>(sampleIndex)
                    / static_cast<double>(sampleCount - 1);
                const double stationA =
                    startNode.stationA + stationSpanA * local;
                const double stationB =
                    startNode.stationB + stationSpanB * local;

                gp_Pnt pointA;
                gp_Pnt pointB;
                if (!sampleChainPoint(
                        candidate.chainA,
                        stationA,
                        pointA) ||
                    !sampleChainPoint(
                        candidate.chainB,
                        stationB,
                        pointB)) {
                    errorMessage = "target edge construction failure: failed to sample correspondence interval";
                    return false;
                }
                if (!GeometryUtils::IsFinitePoint(pointA) ||
                    !GeometryUtils::IsFinitePoint(pointB) ||
                    pointA.Distance(pointB) >
                        params_.maxGapWidth
                            + params_.midlineRepairTolerance) {
                    errorMessage = "target edge construction failure: interval samples exceed allowed gap width";
                    return false;
                }

                // A、B 使用各自站点插值，表达投影对应而不是错误的相同 t 对应。
                const gp_Pnt midpoint(
                    0.5 * (pointA.X() + pointB.X()),
                    0.5 * (pointA.Y() + pointB.Y()),
                    0.5 * (pointA.Z() + pointB.Z()));
                if (!GeometryUtils::IsFinitePoint(midpoint)) {
                    errorMessage = "target edge construction failure: sample midpoint is invalid";
                    return false;
                }
                if (guidancePoints.back().Distance(midpoint) >
                        Precision::Confusion() &&
                    midpoint.Distance(endPoint) >
                        Precision::Confusion()) {
                    guidancePoints.push_back(midpoint);
                }
            }

            if (guidancePoints.back().Distance(endPoint) <=
                Precision::Confusion()) {
                errorMessage = "target edge construction failure: guidance points collapse at endpoint";
                return false;
            }
            guidancePoints.push_back(endPoint);

            Handle(TColgp_HArray1OfPnt) interpolationPoints =
                new TColgp_HArray1OfPnt(
                    1,
                    static_cast<Standard_Integer>(
                        guidancePoints.size()));
            for (std::size_t pointIndex = 0;
                 pointIndex < guidancePoints.size();
                 ++pointIndex) {
                interpolationPoints->SetValue(
                    static_cast<Standard_Integer>(pointIndex + 1),
                    guidancePoints[pointIndex]);
            }

            GeomAPI_Interpolate interpolator(
                interpolationPoints,
                Standard_False,
                Precision::Confusion());
            interpolator.Perform();
            if (!interpolator.IsDone() ||
                interpolator.Curve().IsNull()) {
                errorMessage = "target edge construction failure: interpolation failed";
                return false;
            }

            BRepBuilderAPI_MakeEdge edgeMaker(
                interpolator.Curve(),
                midlineData.chain.vertices[targetEdgeIndex],
                midlineData.chain.vertices[targetEdgeIndex + 1]);
            if (!edgeMaker.IsDone()) {
                errorMessage = "target edge construction failure: edge builder failed";
                return false;
            }

            TopoDS_Edge targetEdge = TopologyUtils::OrientEdgeFromTo(
                edgeMaker.Edge(),
                midlineData.chain.vertices[targetEdgeIndex],
                midlineData.chain.vertices[targetEdgeIndex + 1]);
            TopoDS_Vertex targetStart;
            TopoDS_Vertex targetEnd;
            if (targetEdge.IsNull() ||
                !TopologyUtils::GetEdgeVertices(
                    targetEdge,
                    targetStart,
                    targetEnd) ||
                !targetStart.IsSame(
                    midlineData.chain.vertices[targetEdgeIndex]) ||
                !targetEnd.IsSame(
                    midlineData.chain.vertices[
                        targetEdgeIndex + 1])) {
                errorMessage = "target edge construction failure: canonical edge orientation is invalid";
                return false;
            }

            const double targetEdgeLength =
                GeometryUtils::ComputeEdgeLength(targetEdge);
            if (!std::isfinite(targetEdgeLength) ||
                targetEdgeLength <= Precision::Confusion()) {
                errorMessage = "target edge construction failure: canonical edge is degenerate";
                return false;
            }

            // canonical target edge 只创建一次，后续由 A、B 两侧 face 共同引用。
            midlineData.chain.edges.push_back(targetEdge);
            midlineData.chain.length += targetEdgeLength;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "target edge construction failure: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    if (midlineData.chain.edges.size() + 1 !=
            midlineData.chain.vertices.size() ||
        !std::isfinite(midlineData.chain.length) ||
        midlineData.chain.length <= Precision::Confusion()) {
        errorMessage = "target edge construction failure: canonical target chain is invalid";
        return false;
    }
    return true;
}

bool HealClosedNarrowGap::mergeProjectionNodesForShortTargetEdges(
    std::vector<ProjectionCorrespondenceNode>& nodes,
    const MidlineChainData& midlineData,
    bool& nodesMerged,
    std::string& errorMessage) const
{
    nodesMerged = false;
    if (nodes.size() < 2 ||
        midlineData.chain.vertices.size() != nodes.size() ||
        midlineData.chain.edges.size() + 1 != nodes.size()) {
        errorMessage = "short target edge merge failure: target topology is inconsistent";
        return false;
    }

    try {
        for (std::size_t edgeIndex = 0;
             edgeIndex < midlineData.chain.edges.size();
             ++edgeIndex) {
            const TopoDS_Edge& targetEdge =
                midlineData.chain.edges[edgeIndex];
            const TopoDS_Vertex& startVertex =
                midlineData.chain.vertices[edgeIndex];
            const TopoDS_Vertex& endVertex =
                midlineData.chain.vertices[edgeIndex + 1];
            const double edgeLength =
                GeometryUtils::ComputeEdgeLength(targetEdge);
            const double endpointTolerance =
                std::max(Precision::Confusion(),
                    BRep_Tool::Tolerance(startVertex)
                        + BRep_Tool::Tolerance(endVertex));
            const bool mayDegenerate =
                BRep_Tool::Degenerated(targetEdge) == Standard_True ||
                !std::isfinite(edgeLength) ||
                edgeLength <= endpointTolerance;
            if (!mayDegenerate) {
                continue;
            }

            // 端点节点代表两侧 chain 的公共尖点，不能为了消除短边而吞并内部拓扑。
            if (edgeIndex == 0 ||
                edgeIndex + 2 == nodes.size()) {
                std::ostringstream os;
                os << "short target edge merge failure: near-degenerate edge "
                   << edgeIndex
                   << " touches a chain endpoint; length=" << edgeLength
                   << ", endpointTolerance=" << endpointTolerance;
                errorMessage = os.str();
                return false;
            }

            const ProjectionCorrespondenceNode& firstNode =
                nodes[edgeIndex];
            const ProjectionCorrespondenceNode& secondNode =
                nodes[edgeIndex + 1];
            if ((firstNode.hasSourceVertexA &&
                 secondNode.hasSourceVertexA &&
                 firstNode.sourceVertexIndexA !=
                     secondNode.sourceVertexIndexA) ||
                (firstNode.hasSourceVertexB &&
                 secondNode.hasSourceVertexB &&
                 firstNode.sourceVertexIndexB !=
                     secondNode.sourceVertexIndexB)) {
                std::ostringstream os;
                os << "short target edge merge failure: edge "
                   << edgeIndex
                   << " cannot be removed without collapsing distinct "
                      "source vertices on the same chain; length="
                   << edgeLength
                   << ", endpointTolerance=" << endpointTolerance;
                errorMessage = os.str();
                return false;
            }

            const bool preferFirst =
                firstNode.snappedPair ||
                (!secondNode.snappedPair &&
                 firstNode.projectionDistance <=
                     secondNode.projectionDistance);
            ProjectionCorrespondenceNode mergedNode =
                preferFirst ? firstNode : secondNode;

            // 每一侧都优先保留该侧真实 source vertex 的站点和坐标。
            // 对互补的 A-only/B-only 节点，这会恢复成一个完整的 A/B 对应节点。
            if (firstNode.hasSourceVertexA) {
                mergedNode.stationA = firstNode.stationA;
                mergedNode.pointA = firstNode.pointA;
                mergedNode.sourceVertexIndexA =
                    firstNode.sourceVertexIndexA;
                mergedNode.hasSourceVertexA = true;
            }
            else if (secondNode.hasSourceVertexA) {
                mergedNode.stationA = secondNode.stationA;
                mergedNode.pointA = secondNode.pointA;
                mergedNode.sourceVertexIndexA =
                    secondNode.sourceVertexIndexA;
                mergedNode.hasSourceVertexA = true;
            }

            if (firstNode.hasSourceVertexB) {
                mergedNode.stationB = firstNode.stationB;
                mergedNode.pointB = firstNode.pointB;
                mergedNode.sourceVertexIndexB =
                    firstNode.sourceVertexIndexB;
                mergedNode.hasSourceVertexB = true;
            }
            else if (secondNode.hasSourceVertexB) {
                mergedNode.stationB = secondNode.stationB;
                mergedNode.pointB = secondNode.pointB;
                mergedNode.sourceVertexIndexB =
                    secondNode.sourceVertexIndexB;
                mergedNode.hasSourceVertexB = true;
            }

            mergedNode.endpoint = false;
            mergedNode.snappedPair =
                firstNode.snappedPair || secondNode.snappedPair;
            mergedNode.mergedForShortEdge = true;
            mergedNode.projectionDistance =
                mergedNode.pointA.Distance(mergedNode.pointB);
            if (!std::isfinite(mergedNode.projectionDistance) ||
                mergedNode.projectionDistance >
                    params_.maxGapWidth
                        + params_.midlineRepairTolerance) {
                errorMessage = "short target edge merge failure: merged node exceeds allowed gap width";
                return false;
            }

            const ProjectionCorrespondenceNode& previousNode =
                nodes[edgeIndex - 1];
            const ProjectionCorrespondenceNode& nextNode =
                nodes[edgeIndex + 2];
            if (mergedNode.stationA - previousNode.stationA <=
                    kStationMergeTolerance ||
                mergedNode.stationB - previousNode.stationB <=
                    kStationMergeTolerance ||
                nextNode.stationA - mergedNode.stationA <=
                    kStationMergeTolerance ||
                nextNode.stationB - mergedNode.stationB <=
                    kStationMergeTolerance) {
                errorMessage = "short target edge merge failure: merged node breaks correspondence monotonicity";
                return false;
            }

            nodes[edgeIndex] = mergedNode;
            nodes.erase(nodes.begin()
                + static_cast<std::ptrdiff_t>(edgeIndex + 1));
            nodesMerged = true;
            return true;
        }
    }
    catch (const Standard_Failure& failure) {
        errorMessage = "short target edge merge failure: "
            + RepairDiagnostics::FailureText(failure);
        return false;
    }

    return true;
}

bool HealClosedNarrowGap::buildSourceToTargetEdgeMapping(
    const std::vector<double>& sourceVertexStations,
    const std::vector<double>& targetStations,
    const char* chainName,
    std::vector<std::vector<std::size_t>>& mapping,
    std::string& errorMessage) const
{
    mapping.clear();
    const std::string name = chainName != nullptr ? chainName : "chain";
    if (sourceVertexStations.size() < 2 || targetStations.size() < 2) {
        errorMessage = "mapping failure: " + name
            + " station sequence is incomplete";
        return false;
    }

    mapping.resize(sourceVertexStations.size() - 1);
    for (std::size_t targetEdgeIndex = 0;
         targetEdgeIndex + 1 < targetStations.size();
         ++targetEdgeIndex) {
        const double intervalStart = targetStations[targetEdgeIndex];
        const double intervalEnd = targetStations[targetEdgeIndex + 1];
        const double midpoint = 0.5 * (intervalStart + intervalEnd);
        if (!std::isfinite(midpoint) ||
            midpoint <= intervalStart ||
            midpoint >= intervalEnd) {
            errorMessage = "mapping failure: target edge interval is invalid";
            return false;
        }

        const auto upper = std::upper_bound(
            sourceVertexStations.begin(),
            sourceVertexStations.end(),
            midpoint);
        if (upper == sourceVertexStations.begin() ||
            upper == sourceVertexStations.end()) {
            errorMessage = "mapping failure: " + name
                + " cannot locate the owner source edge for a target edge";
            return false;
        }

        const std::size_t sourceEdgeIndex = static_cast<std::size_t>(
            std::distance(sourceVertexStations.begin(), upper) - 1);
        if (sourceEdgeIndex >= mapping.size() ||
            midpoint <= sourceVertexStations[sourceEdgeIndex] ||
            midpoint >= sourceVertexStations[sourceEdgeIndex + 1]) {
            errorMessage = "mapping failure: " + name
                + " target edge midpoint is outside its source interval";
            return false;
        }
        mapping[sourceEdgeIndex].push_back(targetEdgeIndex);
    }

    return true;
}

bool HealClosedNarrowGap::validateReplacementMapping(
    const std::vector<std::vector<std::size_t>>& mapping,
    std::size_t sourceEdgeCount,
    std::size_t targetEdgeCount,
    const char* chainName,
    std::string& errorMessage) const
{
    const std::string name = chainName != nullptr ? chainName : "chain";
    if (mapping.size() != sourceEdgeCount) {
        errorMessage = "mapping failure: " + name
            + " replacement mapping count does not match source edges";
        return false;
    }
    if (targetEdgeCount == 0) {
        errorMessage = "mapping failure: canonical target edge chain is empty";
        return false;
    }

    std::vector<int> targetOwnership(targetEdgeCount, 0);
    std::size_t nextExpectedTargetIndex = 0;
    for (std::size_t sourceEdgeIndex = 0;
         sourceEdgeIndex < mapping.size();
         ++sourceEdgeIndex) {
        const std::vector<std::size_t>& replacements =
            mapping[sourceEdgeIndex];
        if (replacements.empty()) {
            errorMessage = "mapping failure: " + name
                + " source edge has no target edge";
            return false;
        }

        for (std::size_t replacementIndex = 0;
             replacementIndex < replacements.size();
             ++replacementIndex) {
            const std::size_t targetEdgeIndex =
                replacements[replacementIndex];
            if (targetEdgeIndex >= targetEdgeCount) {
                errorMessage = "mapping failure: " + name
                    + " target edge index is out of range";
                return false;
            }
            if (targetEdgeIndex != nextExpectedTargetIndex) {
                errorMessage = "mapping failure: " + name
                    + " target edges are not ordered by source edge";
                return false;
            }
            if (replacementIndex > 0 &&
                targetEdgeIndex != replacements[replacementIndex - 1] + 1) {
                errorMessage = "mapping failure: " + name
                    + " replacement indices are not continuous and increasing";
                return false;
            }
            ++targetOwnership[targetEdgeIndex];
            ++nextExpectedTargetIndex;
        }
    }

    if (nextExpectedTargetIndex != targetEdgeCount) {
        errorMessage = "mapping failure: " + name
            + " does not cover the full canonical target chain";
        return false;
    }
    for (std::size_t targetEdgeIndex = 0;
         targetEdgeIndex < targetOwnership.size();
         ++targetEdgeIndex) {
        if (targetOwnership[targetEdgeIndex] != 1) {
            errorMessage = "mapping failure: " + name
                + " does not cover every target edge exactly once";
            return false;
        }
    }
    return true;
}

bool HealClosedNarrowGap::buildMidlineChainByProjectionCorrespondence(
    const ClosedGapCandidate& candidate,
    MidlineChainData& midlineData,
    std::string& errorMessage) const
{
    midlineData = MidlineChainData();
    midlineData.constructionMode =
        MidlineConstructionMode::ProjectionCorrespondence;

    ChainStationData chainAStations;
    ChainStationData chainBStations;
    if (!buildNormalizedVertexStations(
            candidate.chainA,
            "chainA",
            chainAStations,
            errorMessage) ||
        !buildNormalizedVertexStations(
            candidate.chainB,
            "chainB",
            chainBStations,
            errorMessage)) {
        return false;
    }

    std::vector<ProjectionCorrespondenceNode> nodes;
    if (!buildProjectionCorrespondenceNodes(
            candidate,
            chainAStations,
            chainBStations,
            nodes,
            errorMessage)) {
        return false;
    }

    const std::size_t initialNodeCount = nodes.size();
    while (true) {
        if (!buildProjectionTargetVertices(
                candidate,
                nodes,
                midlineData,
                errorMessage) ||
            !buildProjectionTargetEdges(
                candidate,
                chainAStations,
                chainBStations,
                nodes,
                midlineData,
                errorMessage)) {
            return false;
        }

        bool nodesMerged = false;
        if (!mergeProjectionNodesForShortTargetEdges(
                nodes,
                midlineData,
                nodesMerged,
                errorMessage)) {
            return false;
        }
        if (!nodesMerged) {
            break;
        }

        ++midlineData.shortEdgeMergeCount;
        if (nodes.size() < 2 ||
            midlineData.shortEdgeMergeCount >
                static_cast<int>(initialNodeCount - 2)) {
            errorMessage = "short target edge merge failure: merge iteration limit exceeded";
            return false;
        }
    }
    midlineData.correspondenceNodes = nodes;

    std::vector<double> targetStationsA;
    std::vector<double> targetStationsB;
    targetStationsA.reserve(nodes.size());
    targetStationsB.reserve(nodes.size());
    for (const ProjectionCorrespondenceNode& node : nodes) {
        targetStationsA.push_back(node.stationA);
        targetStationsB.push_back(node.stationB);
    }

    // 一条 source edge 可以跨越多个对应区间，但每个 target edge
    // 在同一侧只能归属于一条 source edge。
    if (!buildSourceToTargetEdgeMapping(
            chainAStations.vertexStations,
            targetStationsA,
            "chainA",
            midlineData.replacementsForChainA,
            errorMessage) ||
        !buildSourceToTargetEdgeMapping(
            chainBStations.vertexStations,
            targetStationsB,
            "chainB",
            midlineData.replacementsForChainB,
            errorMessage)) {
        return false;
    }

    if (!validateReplacementMapping(
            midlineData.replacementsForChainA,
            candidate.chainA.edges.size(),
            midlineData.chain.edges.size(),
            "chainA",
            errorMessage) ||
        !validateReplacementMapping(
            midlineData.replacementsForChainB,
            candidate.chainB.edges.size(),
            midlineData.chain.edges.size(),
            "chainB",
            errorMessage)) {
        return false;
    }

    return true;
}

const char* HealClosedNarrowGap::midlineConstructionModeName(
    MidlineConstructionMode mode)
{
    switch (mode) {
    case MidlineConstructionMode::Auto:
        return "Auto";
    case MidlineConstructionMode::LegacyOneToOne:
        return "LegacyOneToOne";
    case MidlineConstructionMode::ProjectionCorrespondence:
        return "ProjectionCorrespondence";
    default:
        return "Unknown";
    }
}

const char* HealClosedNarrowGap::gapEndTypeName(GapEndType type)
{
    switch (type) {
    case GapEndType::SharedTip:
        return "SharedTip";
    case GapEndType::ShortEdge:
        return "ShortEdge";
    default:
        return "Unknown";
    }
}

