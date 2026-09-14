#pragma once
#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

// 两个原始支撑面上的局部重裁。保留等长短边的公共端点，删除错误边界链，
// 将 P-Q 连线投影到候选 A 的原 surface，使用同一共享 edge 重建两侧 face。
class HealIntersectingFaceGap {
public:
    struct Params {
        double maxGapWidth = 0.5;             // 投影边的缺陷带覆盖范围；不检查两侧原边界宽度。
        double maxEndpointCollapseDistance = 0.8; // 可删除的端部短边长度预算。
        double endpointLengthDifferenceTolerance = 1.0e-5; // 相邻短自由边长度差，模型单位。
        double pcurveTolerance = 1.0e-1;     // 重建面 pcurve 的计算精度，独立于修复容差预算。
        double maxRebuiltTolerance = 1.0e-1; // 与缺陷宽度独立，不用于强行连接。
        double maxRelativeAreaChange = 0.15;
        int maxBoundaryEdges = 8;            // 单侧缺陷链复杂度上限。
    };
    HealIntersectingFaceGap();
    explicit HealIntersectingFaceGap(const Params& params);
    int perform(const TopoDS_Shape& shape);
    TopoDS_Shape getResultShape() const { return result_; }
    bool isModified() const { return modified_; }
    const std::vector<std::string>& getMessages() const { return messages_; }
private:
    Params params_;
    TopoDS_Shape result_;
    bool modified_ = false;
    std::vector<std::string> messages_;
};
