#include "MergeDuplicateVertex.h"

#include "GeometryUtils.h"
#include "TopologyUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListIteratorOfListOfStatus.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <BRepTools_ReShape.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>

namespace {
constexpr int kMaxAutomaticMergeIterations = 1000;

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
} // namespace

MergeDuplicateVertex::MergeDuplicateVertex(const Params& params)
    : params_(params)
{
    params_.zeroDistanceTolerance = std::max(params_.zeroDistanceTolerance, Precision::Confusion());
}

void MergeDuplicateVertex::setInputShape(const TopoDS_Shape& shape)
{
    inputShape_ = shape;
    resultShape_ = inputShape_;
    modified_ = false;
    mergedCount_ = 0;
    messages_.clear();
}

int MergeDuplicateVertex::perform()
{
    modified_ = false;
    mergedCount_ = 0;
    messages_.clear();
    resultShape_ = inputShape_;

    if (inputShape_.IsNull()) {
        appendMessage("Input shape is null");
        return 0;
    }

    TopoDS_Shape currentShape = inputShape_;
    const int iterationLimit = computeAutomaticIterationLimit(currentShape);
    {
        std::ostringstream os;
        os << "Automatic duplicate-vertex merge iteration limit=" << iterationLimit
           << ", zeroDistanceTolerance=" << params_.zeroDistanceTolerance;
        appendMessage(os.str());
    }

    for (int iter = 0; iter < iterationLimit; ++iter) {
        std::vector<Candidate> candidates = findDuplicateVertexCandidates(currentShape);
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
                return lhs.confidence > rhs.confidence;
            });

        {
            std::ostringstream os;
            os << "Iteration " << (iter + 1)
               << ": duplicate vertex candidate count=" << candidates.size();
            appendMessage(os.str());
        }

        if (candidates.empty()) {
            break;
        }

        bool mergedThisIteration = false;
        for (const Candidate& candidate : candidates) {
            TopoDS_Shape mergedShape;
            std::string errorMessage;
            if (!tryMergeCandidate(currentShape, candidate, mergedShape, errorMessage)) {
                appendMessage("Candidate failed: " + errorMessage);
                continue;
            }

            currentShape = mergedShape;
            resultShape_ = currentShape;
            modified_ = true;
            ++mergedCount_;
            mergedThisIteration = true;

            std::ostringstream os;
            os << "Merged duplicate vertex pair: distance=" << candidate.distance
               << ", confidence=" << candidate.confidence;
            appendMessage(os.str());
            break;
        }

        if (!mergedThisIteration) {
            break;
        }

        if (iter + 1 == iterationLimit) {
            appendMessage("Merge stopped: automatic iteration limit reached");
        }
    }

    resultShape_ = currentShape;
    return mergedCount_;
}

int MergeDuplicateVertex::perform(const TopoDS_Shape& shape)
{
    setInputShape(shape);
    return perform();
}

std::vector<MergeDuplicateVertex::Candidate> MergeDuplicateVertex::findDuplicateVertexCandidates(
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

    std::vector<int> nearZeroPartnerCount(static_cast<std::size_t>(vertexMap.Extent()) + 1, 0);
    std::vector<RawPair> rawPairs;
    int multiVertexPairCount = 0;

    for (int i = 1; i <= vertexMap.Extent(); ++i) {
        const TopoDS_Vertex firstVertex = TopoDS::Vertex(vertexMap.FindKey(i));
        if (firstVertex.IsNull()) {
            continue;
        }

        for (int j = i + 1; j <= vertexMap.Extent(); ++j) {
            const TopoDS_Vertex secondVertex = TopoDS::Vertex(vertexMap.FindKey(j));
            if (secondVertex.IsNull() || firstVertex.IsSame(secondVertex)) {
                continue;
            }

            double distance = 0.0;
            if (!isZeroDistanceVertexPair(firstVertex, secondVertex, distance)) {
                continue;
            }

            rawPairs.push_back({ i, j, distance });
            ++nearZeroPartnerCount[static_cast<std::size_t>(i)];
            ++nearZeroPartnerCount[static_cast<std::size_t>(j)];
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

    if (multiVertexPairCount > 0) {
        std::ostringstream os;
        os << "Skipped near-zero pairs involved in 3+ vertex groups: " << multiVertexPairCount;
        appendMessage(os.str());
    }

    return candidates;
}

bool MergeDuplicateVertex::isZeroDistanceVertexPair(
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

TopoDS_Vertex MergeDuplicateVertex::chooseKeptVertex(
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

bool MergeDuplicateVertex::isSafeToMergePair(
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

bool MergeDuplicateVertex::tryMergeCandidate(
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

bool MergeDuplicateVertex::repairAffectedTopology(
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

bool MergeDuplicateVertex::validateMergeResult(
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
        errorMessage = "shape type changed after duplicate vertex merge";
        return false;
    }

    const int freeEdgesBefore = TopologyUtils::CountFreeEdges(beforeShape);
    const int freeEdgesAfter = TopologyUtils::CountFreeEdges(afterShape);
    if (freeEdgesAfter > freeEdgesBefore) {
        errorMessage = "free edge count increased after duplicate vertex merge";
        return false;
    }

    const int collapsedEdgesBefore = countCollapsedEdges(beforeShape, collapsedEdgeTolerance);
    const int collapsedEdgesAfter = countCollapsedEdges(afterShape, collapsedEdgeTolerance);
    if (collapsedEdgesAfter > collapsedEdgesBefore) {
        errorMessage = "collapsed edge count increased after duplicate vertex merge";
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
        errorMessage = "new BRepCheck failures after duplicate vertex merge:";
        errorMessage += describeBRepCheckFailures(newFailures);
        if (!baselineFailures.empty()) {
            errorMessage += " baseline failures:";
            errorMessage += describeBRepCheckFailures(baselineFailures);
        }
        return false;
    }

    return true;
}

int MergeDuplicateVertex::computeAutomaticIterationLimit(const TopoDS_Shape& shape) const
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

int MergeDuplicateVertex::countCollapsedEdges(const TopoDS_Shape& shape, double tolerance) const
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

void MergeDuplicateVertex::appendMessage(const std::string& message)
{
    messages_.push_back(message);
}
