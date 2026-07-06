#pragma once

#include <string>
#include <vector>

#include <Interface_Static.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelMap.hxx>
#include <TDocStd_Application.hxx>
#include <TDocStd_Document.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>

#include "MergeShortEdge.h"
#include "HealNarrowGap.h"
#include "HealClosedNarrowGap.h"
#include "MergeDuplicateVertex.h"
#include "SewShells.h"

class ApplicationInit {
public:
    ApplicationInit();
    ~ApplicationInit();

    Handle(TDocStd_Application) getApp() const { return anApp_; }
    Handle(TDocStd_Document) getDoc() const { return aDoc_; }
    TopoDS_Shape getImportedStepShape() const { return importedStepShape_; }

    bool importStep(const std::string& filePath);
    bool importIges(const std::string& filePath);
    // 使用 STEPControl_Reader 将 STEP 文件直接读成一个 TopoDS_Shape，不进入 XCAF 文档。
    bool importStepAsShape(const std::string& filePath);

    bool exportStep(const TopoDS_Shape& shape, const std::string& filePath) const;
    bool exportStepShape(const TopoDS_Shape& shape, const std::string& filePath) const;
    int mergeDuplicateVerticesOnStepShapeAndExport(
        const std::string& filePath,
        const MergeDuplicateVertex::Params& params = MergeDuplicateVertex::Params());
    int sewShellsOnStepShapeAndExport(
        const std::string& filePath,
        const SewShells::Params& params = SewShells::Params());
    int mergeShortEdgesOnLowestContainersAndExport(const std::string& filePath,
                                                   bool useManualLengthThreshold,
                                                   double manualLengthThreshold,
                                                   double angleThresholdDeg);
    int healNarrowGapsOnLowestContainersAndExport(const std::string& filePath,
                                                  const HealNarrowGap::Params& params = HealNarrowGap::Params());
    int healClosedGapsOnLowestContainersAndExport(
        const std::string& filePath,
        const HealClosedNarrowGap::Params& params = HealClosedNarrowGap::Params());
    int mergeDuplicateVerticesOnLowestContainersAndExport(
        const std::string& filePath,
        const MergeDuplicateVertex::Params& params = MergeDuplicateVertex::Params());
    int mergeDuplicateVerticesAndHealNarrowGapsOnLowestContainersAndExport(
        const std::string& filePath,
        const MergeDuplicateVertex::Params& duplicateVertexParams = MergeDuplicateVertex::Params(),
        const HealNarrowGap::Params& narrowGapParams = HealNarrowGap::Params());

    std::vector<TopoDS_Shape> getAllShapes() const;
    std::vector<TopoDS_Solid> getAllSolids() const;
    std::vector<TopoDS_Shell> getAllShells() const;
    std::vector<TopoDS_Edge> getAllEdges() const;

private:
    // 对 originalShape 整体执行短边合并，再一次性写回
    // ownerLabel，避免大型装配中按 shell 逐个 Replace 的高开销。
    struct EditableContainerTarget {
        TDF_Label ownerLabel;
        TopoDS_Shape originalShape;
    };

    bool transferAndCacheTopLevelShapes(Standard_Boolean isTransferred, const char* fileType);
    bool writeXcafDocumentStep(const Handle(TDocStd_Document)& doc,
                               const std::string& filePath,
                               bool preserveAssembly) const;
    void collectTopLevelTargets(std::vector<EditableContainerTarget>& targets) const;
    // 递归收集辅助函数：解析 reference，向下遍历嵌套装配，并在遇到
    // 第一个可整体处理的容器时停止继续拆分。
    // 判断 shape 自身是否是 solid/shell，或内部是否包含 solid/shell 拓扑。
    bool shapeContainsSolidOrShell(const TopoDS_Shape& shape) const;
    void setOutputStepUnitFromInput() const;
    void refreshShapeCache();
    void clearShapeCache();

    std::vector<TopoDS_Shape> aShapeVec_;
    std::vector<TopoDS_Shell> aShellVec_;
    std::vector<TopoDS_Solid> aSolidVec_;
    TopoDS_Shape importedStepShape_;
    std::string inputStepUnit_ = "MM";
    Handle(TDocStd_Application) anApp_;
    Handle(TDocStd_Document) aDoc_;
};
