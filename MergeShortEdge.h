#pragma once

#include <vector>

#include <Precision.hxx>
#include <Standard_Handle.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

class Geom_Curve;

struct MergeCandidate {
    TopoDS_Edge edge1;//edge1通常是短边
    TopoDS_Edge edge2;
    TopoDS_Vertex commonVertex;
    double length1;
    double length2;
    double angle;

    MergeCandidate() : length1(0.0), length2(0.0), angle(0.0) {}
};

enum class ShortEdgeThresholdMode {
    ManualAbsolute,
    BoundingBoxRatio,
    HybridGlobalAndLocal
};

class MergeShortEdge {
public:
    MergeShortEdge(double bboxRatio = 0.001, double angleThresholdDeg = 175.0);
    ~MergeShortEdge() = default;

    MergeShortEdge(const MergeShortEdge&) = delete;
    MergeShortEdge& operator=(const MergeShortEdge&) = delete;
    MergeShortEdge(MergeShortEdge&&) = default;
    MergeShortEdge& operator=(MergeShortEdge&&) = default;

    void setInputShape(const TopoDS_Shape& shape);
    void setInputShapes(const std::vector<TopoDS_Shape>& shapes);

    void setLengthThreshold(double threshold) {
        lengthThreshold_ = threshold;
        thresholdMode_ = ShortEdgeThresholdMode::ManualAbsolute;
    }

    void setThresholdMode(ShortEdgeThresholdMode mode) { thresholdMode_ = mode; }
    ShortEdgeThresholdMode getThresholdMode() const { return thresholdMode_; }

    bool setLengthThresholdByBoundingBox(double ratio = 0.001);
    bool setLengthThresholdHybrid();

    void setAngleThresholdDeg(double degrees) { angleThresholdDeg_ = degrees; }

    double getLengthThreshold() const { return lengthThreshold_; }
    double getAngleThresholdDeg() const { return angleThresholdDeg_; }

    std::vector<TopoDS_Edge> getShortEdges();
    std::vector<MergeCandidate> findMergeCandidates();
    int performMerge();

    TopoDS_Shape getResultShape() const { return resultShape_; }
    const std::vector<TopoDS_Edge>& getShortEdges() const { return shortEdges_; }
    const std::vector<MergeCandidate>& getMergeCandidates() const { return mergeCandidates_; }
    const TopTools_IndexedDataMapOfShapeListOfShape& getVertexEdgeMap() const { return vertexEdgeMap_; }

private:
    double lengthThreshold_;
    double angleThresholdDeg_;
    double bboxRatio_;
    double localLengthRatio_; // 混合阈值模式下，短边还必须显著短于其局部邻边。
    ShortEdgeThresholdMode thresholdMode_; // 当前阈值来源：手工指定、包围盒比例或统计+局部混合模式。

    std::vector<TopoDS_Shape> inputShapes_;
    TopTools_IndexedDataMapOfShapeListOfShape vertexEdgeMap_;
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap_;

    std::vector<TopoDS_Edge> shortEdges_;
    std::vector<MergeCandidate> mergeCandidates_;
    TopoDS_Shape resultShape_;

    TopoDS_Shape preprocessWithShapeFixWireframe(const TopoDS_Shape& shape) const;
    //预处理
    TopoDS_Shape preprocessInputShape(const TopoDS_Shape& shape) const;
    TopoDS_Shape executeUnifySameDomain(const TopoDS_Shape& shape) const;
    // 按当前模式刷新有效阈值，避免手工阈值被后续流程隐式覆盖。
    bool refreshLengthThreshold();
    // 从全局边长分布的下尾部估计一个更稳健的自动短边阈值。
    double computeStatisticBasedLengthThreshold() const;
    // 混合模式下，要求该边在局部连接关系里也明显偏短，避免误伤正常小特征。
    bool satisfiesLocalLengthCriterion(const TopoDS_Edge& edge) const;
    // 校验新构造边的几何是否仍与两端顶点匹配，减少坏边进入后续拓扑修复。
    bool validateMergedEdgeGeometry(const TopoDS_Edge& edge) const;
    // 判断新边在指定面上是否已经拥有可用的 pcurve。
    bool hasValidPCurveOnFace(const TopoDS_Edge& edge, const TopoDS_Face& face) const;
    // 针对引用该边的所有面补建并校验 pcurve，同时修正 same-parameter。
    bool repairPCurvesForEdge(TopoDS_Edge& edge, const TopTools_ListOfShape& faces) const;
    // 对本轮新生成的边统一补建受影响面的 pcurve。
    bool repairAffectedPCurves(TopoDS_Shape& shape, std::vector<TopoDS_Edge>& newEdges) const;
    // 对合并后的整体结果做最终合法性检查，失败时用于回退。
    bool validateShape(const TopoDS_Shape& shape) const;
    // 取边的底层 3D 曲线并消除 Location 与裁剪曲线包装。
    Handle(Geom_Curve) getUnderlyingCurve(const TopoDS_Edge& edge) const;
    // 两条边顺序复合成一条 B 样条以尽量保留连接处细节。
    TopoDS_Edge createMergedEdgeByCompositeBSpline(const MergeCandidate& candidate) const;
    // 探测短边
    void detectShortEdges();
    TopoDS_Edge createMergedEdge(const MergeCandidate& candidate) const;
    TopoDS_Vertex getOtherVertex(const TopoDS_Edge& edge, const TopoDS_Vertex& commonVertex) const;
	// 根据当前候选列表构建合并批次，优先合并互不干扰的边以提升效率和结果质量。
    std::vector<std::vector<MergeCandidate>> buildCandidateBatches(const std::vector<MergeCandidate>& candidates) const;
    bool tryApplyBatch(const TopoDS_Shape& baseShape,
                       const std::vector<MergeCandidate>& batch,
                       TopoDS_Shape& committedShape,
                       int& committedCount) const;
    void buildTopologyMap();
    // 排除不适合合并的侯选边
    bool isHighRiskEdge(const TopoDS_Edge& edge) const;
    bool isShortEdge(const TopoDS_Edge& edge) const;
    bool findBestAdjacentEdge(const TopoDS_Edge& shortEdge, TopoDS_Vertex& bestCommonVertex, TopoDS_Edge& bestAdjacentEdge) const;
    std::vector<TopoDS_Edge> collectAllEdges() const;
};
