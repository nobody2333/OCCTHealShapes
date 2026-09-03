#include <topo/Body.h>

#include <vector>
#include <format>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <TopExp_Explorer.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopAbs_ShapeEnum.hxx>

#include "ApplicationInit.h"
#include "GeometryUtils.h"
#include "ShapeRepairPipeline.h"
#include "ShapeRepairStages.h"
#include "SewShells.h"
#include "TopoTools.h"
#include "Viewer.h"

//连续模型转换为离散模型
void modelConvert(const std::string& fileName) {
    std::string filePath = "./assets/" + fileName + ".stp";
    TopoTools aTopoTools(0.5, fileName);
    auto startTime = std::chrono::system_clock::now();

    aTopoTools.importSTEP(filePath);
    std::vector<Body*> bodies = aTopoTools.getBodies();
    std::vector<Face*> faces;
    int faceNumber = 0;
    for (auto& body : bodies)
    {
        for (auto& face : body->faces)
        {
            faces.push_back(face);
            faceNumber++;
        }
    }

    //std::cout << std::format("Face number in bodies: {}", faceNumber) << std::endl;
    aTopoTools.exportOBJ("./results/" + fileName, faces);
    aTopoTools.exportVTK("./results/" + fileName, faces);

    auto endTime = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
    std::cout << "Time cost: " << duration.count() << " seconds" << std::endl;
}

void computeLength(const std::vector<TopoDS_Shape>& shapes) {
    for (std::size_t shapeIndex = 0; shapeIndex < shapes.size(); ++shapeIndex) {
        const TopoDS_Shape& shape = shapes[shapeIndex];
        if (shape.IsNull()) {
            std::cout << std::format("Shape {} is null, skip.", shapeIndex) << std::endl;
            continue;
        }

        TopTools_MapOfShape uniqueEdges;
        int edgeCount = 0;
        std::cout << std::format("Shape {} edge lengths:", shapeIndex) << std::endl;

        for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
            const TopoDS_Edge edge = TopoDS::Edge(exp.Current());
            if (!uniqueEdges.Add(edge)) {
                continue;
            }

            const double length = GeometryUtils::ComputeEdgeLength(edge);
            std::cout << std::format("  Edge {}: {:.9f}", edgeCount, length) << std::endl;
            ++edgeCount;
        }

        std::cout << std::format("Shape {} unique edge count: {}", shapeIndex, edgeCount) << std::endl;
    }
}

void mergeShortEdge(const std::string& fileName) {
    
    std::string inputFilePath = "./assets/" + fileName + ".stp";
    std::string outputFilePath = "./assets/fixedShortEdge/MergeEdge_" + fileName + ".stp";
    ApplicationInit appInit;

	auto t0 = std::chrono::system_clock::now();
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file fasiled" << inputFilePath << std::endl;
        return;
    }

    auto t1 = std::chrono::system_clock::now();
    std::cout << "ReadFile and Transfer cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
        << " s" << std::endl;

    std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return;
    }
    for (std::size_t i = 0; i < inputShapes.size(); ++i) {
        std::cout << std::format("Shape[{}] type: {}", i, static_cast<int>(inputShapes[i].ShapeType())) << std::endl;
    }


    const bool useManualLengthThreshold = false;
    const double manualLengthThreshold = 0.02;
    const double angleThresholdDeg = 160.0;


    
    std::cout << "Merging short edges......" << std::endl;
    int mergedCount = 0;
    mergedCount = appInit.mergeShortEdgesOnLowestContainersAndExport(outputFilePath,
                                                                               useManualLengthThreshold,
                                                                               manualLengthThreshold,
                                                                               angleThresholdDeg);
    if (mergedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return;
    }
    std::cout << "Merged short-edge candidate count: " << mergedCount << std::endl;


    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;
	auto t2 = std::chrono::system_clock::now();
    std::cout << "Total cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count()
        << " s" << std::endl;
    return;
}

int healNarrowGap(const std::string& fileName) {
    std::string inputFilePath = "./assets/" + fileName + ".stp";
    std::filesystem::create_directories("./assets/fixedNarrowGap");
    std::string outputFilePath = "./assets/fixedNarrowGap/HealGap_" + fileName + ".stp";
    ApplicationInit appInit;

    auto t0 = std::chrono::system_clock::now();
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file failed: " << inputFilePath << std::endl;
        return -1;
    }

    auto t1 = std::chrono::system_clock::now();
    std::cout << "ReadFile and Transfer cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
        << " s" << std::endl;

    std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return -1;
    }
    for (std::size_t i = 0; i < inputShapes.size(); ++i) {
        std::cout << std::format("Shape[{}] type: {}", i, static_cast<int>(inputShapes[i].ShapeType())) << std::endl;
    }

    HealNarrowGap::Params params;
    std::cout << "Healing narrow gaps......" << std::endl;
    const int repairedCount = appInit.healNarrowGapsOnLowestContainersAndExport(outputFilePath, params);
    if (repairedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Repaired narrow-gap candidate count: " << repairedCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    auto t2 = std::chrono::system_clock::now();
    std::cout << "Total cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count()
        << " s" << std::endl;

    return repairedCount;
}

int healClosedGap(const std::string& fileName) {
    std::string inputFilePath = "./assets/" + fileName + ".stp";
    std::filesystem::create_directories("./assets/fixedClosedGap");
    std::string outputFilePath = "./assets/fixedClosedGap/HealClosedGap_" + fileName + ".stp";
    ApplicationInit appInit;

    auto t0 = std::chrono::system_clock::now();
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file failed: " << inputFilePath << std::endl;
        return -1;
    }

    auto t1 = std::chrono::system_clock::now();
    std::cout << "ReadFile and Transfer cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
        << " s" << std::endl;

    std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return -1;
    }
    for (std::size_t i = 0; i < inputShapes.size(); ++i) {
        std::cout << std::format("Shape[{}] type: {}", i, static_cast<int>(inputShapes[i].ShapeType())) << std::endl;
    }

    HealClosedNarrowGap::Params params;
	params.midlineConstructionMode = HealClosedNarrowGap::MidlineConstructionMode::ProjectionCorrespondence;
    std::cout << "Healing closed narrow gaps......" << std::endl;
    const int repairedCount = appInit.healClosedGapsOnLowestContainersAndExport(outputFilePath, params);
    if (repairedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Repaired closed-gap candidate count: " << repairedCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    auto t2 = std::chrono::system_clock::now();
    std::cout << "Total cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count()
        << " s" << std::endl;

    return repairedCount;
}

// 单独运行简单孔洞修复；容器遍历和 STEP 导出复用通用流水线入口。
int healSimpleHole(const std::string& fileName) {
    const std::string inputFilePath = "./assets/" + fileName + ".stp";
    const std::filesystem::path outputDirectory = "./assets/fixedSimpleHole";
    std::filesystem::create_directories(outputDirectory);
    const std::string outputFilePath =
        (outputDirectory / ("HealSimpleHole_" + fileName + ".stp")).string();

    ApplicationInit appInit;
    const auto startTime = std::chrono::system_clock::now();
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file failed: " << inputFilePath << std::endl;
        return -1;
    }

    const std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return -1;
    }
    for (std::size_t i = 0; i < inputShapes.size(); ++i) {
        std::cout << std::format(
            "Shape[{}] type: {}",
            i,
            static_cast<int>(inputShapes[i].ShapeType())) << std::endl;
    }

    HealSimpleHole::Params holeParams;
    // 使用 HealSimpleHole::Params 中统一的已验证默认值。自动补孔无法仅凭拓扑
    // 区分设计开口；处理其他模型尺度时仍应按模型单位覆盖尺寸与误差参数。

    ShapeRepairPipeline pipeline;
    pipeline.addStage(
        ShapeRepairStages::MakeHealSimpleHoleStage(holeParams));

    std::cout << "Healing simple holes......" << std::endl;
    const int repairedCount = appInit.runRepairPipelineOnLowestContainersAndExport(
        outputFilePath,
        pipeline);
    if (repairedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Repaired simple-hole candidate count: "
              << repairedCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    const auto endTime = std::chrono::system_clock::now();
    std::cout << "Total cost: "
              << std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count()
              << " s" << std::endl;
    return repairedCount;
}

int mergeDuplicateVertex(const std::string& fileName) {
    std::string inputFilePath = "./assets/" + fileName + ".stp";
    std::filesystem::create_directories("./assets/fixedDuplicateVertex");
    std::string outputFilePath = "./assets/fixedDuplicateVertex/MergeVertex_" + fileName + ".stp";
    ApplicationInit appInit;

    auto t0 = std::chrono::system_clock::now();
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file failed: " << inputFilePath << std::endl;
        return -1;
    }

    auto t1 = std::chrono::system_clock::now();
    std::cout << "ReadFile and Transfer cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
        << " s" << std::endl;

    std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return -1;
    }
    for (std::size_t i = 0; i < inputShapes.size(); ++i) {
        std::cout << std::format("Shape[{}] type: {}", i, static_cast<int>(inputShapes[i].ShapeType())) << std::endl;
    }

    MergeDuplicateVertex::Params params;

    std::cout << "Merging duplicate vertices......" << std::endl;
    const int mergedCount = appInit.mergeDuplicateVerticesOnLowestContainersAndExport(outputFilePath, params);
    if (mergedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Merged duplicate-vertex pair count: " << mergedCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    auto t2 = std::chrono::system_clock::now();
    std::cout << "Total cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count()
        << " s" << std::endl;

    return mergedCount;
}

int mergeDuplicateVertexOnStepShape(const std::string& fileName) {
    std::string inputFilePath = "./assets/" + fileName + ".stp";
    std::filesystem::create_directories("./assets/fixedDuplicateVertexStepShape");
    std::string outputFilePath = "./assets/fixedDuplicateVertexStepShape/MergeVertex_" + fileName + ".stp";
    ApplicationInit appInit;

    auto t0 = std::chrono::system_clock::now();
    if (!appInit.importStepAsShape(inputFilePath)) {
        std::cerr << "Reading file by STEPControl_Reader failed: " << inputFilePath << std::endl;
        return -1;
    }

    auto t1 = std::chrono::system_clock::now();
    std::cout << "STEPControl_Reader cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
        << " s" << std::endl;

    MergeDuplicateVertex::Params params;
    params.forbidAcrossDifferentSolids = false;

    std::cout << "Merging duplicate vertices on imported STEP shape......" << std::endl;
    const int mergedCount = appInit.mergeDuplicateVerticesOnStepShapeAndExport(outputFilePath, params);
    if (mergedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Merged duplicate-vertex pair count: " << mergedCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    auto t2 = std::chrono::system_clock::now();
    std::cout << "Total cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count()
        << " s" << std::endl;

    return mergedCount;
}

int sewShellsOnStepShape(const std::string& fileName) {
    std::string inputFilePath = "./assets/" + fileName + ".stp";
    std::filesystem::create_directories("./assets/fixedSewShellsStepShape");
    std::string outputFilePath = "./assets/fixedSewShellsStepShape/SewShells_" + fileName + ".stp";
    ApplicationInit appInit;

    auto t0 = std::chrono::system_clock::now();
    if (!appInit.importStepAsShape(inputFilePath)) {
        std::cerr << "Reading file by STEPControl_Reader failed: " << inputFilePath << std::endl;
        return -1;
    }

    auto t1 = std::chrono::system_clock::now();
    std::cout << "STEPControl_Reader cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
        << " s" << std::endl;

    SewShells::Params params;
    params.nonManifoldMode = true;
    params.floatingEdgesMode = true;
    params.runShapeFixAfterSewing = true;
	params.rejectInvalidResult = false;
    params.sewingTolerance = 1.0e-2;
    params.maxTolerance = 1.0e-2;
	params.discardLooseTopLevelSubshapesAfterSewing = true;
    std::cout << "Sewing shells on imported STEP shape......" << std::endl;
    const int sewedGroupCount = appInit.sewShellsOnStepShapeAndExport(outputFilePath, params);
    if (sewedGroupCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Sewed shell group count: " << sewedGroupCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    auto t2 = std::chrono::system_clock::now();
    std::cout << "Total cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count()
        << " s" << std::endl;

    return sewedGroupCount;
}

int mergeDuplicateVertexAndHealNarrowGap(const std::string& fileName) {
    std::string inputFilePath = "./assets/" + fileName + ".stp";
    std::filesystem::create_directories("./assets/fixedDuplicateVertexNarrowGap");
    std::string outputFilePath = "./assets/fixedDuplicateVertexNarrowGap/dp&hn_" + fileName + ".stp";
    ApplicationInit appInit;

    auto t0 = std::chrono::system_clock::now();
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file failed: " << inputFilePath << std::endl;
        return -1;
    }

    auto t1 = std::chrono::system_clock::now();
    std::cout << "ReadFile and Transfer cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
        << " s" << std::endl;

    std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return -1;
    }
    for (std::size_t i = 0; i < inputShapes.size(); ++i) {
        std::cout << std::format("Shape[{}] type: {}", i, static_cast<int>(inputShapes[i].ShapeType())) << std::endl;
    }

    MergeDuplicateVertex::Params duplicateVertexParams;
    HealNarrowGap::Params narrowGapParams;

    std::cout << "Merging duplicate vertices, then healing narrow gaps......" << std::endl;
    const int repairedCount = appInit.mergeDuplicateVerticesAndHealNarrowGapsOnLowestContainersAndExport(
        outputFilePath,
        duplicateVertexParams,
        narrowGapParams);
    if (repairedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Duplicate-vertex/narrow-gap repair count: " << repairedCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    auto t2 = std::chrono::system_clock::now();
    std::cout << "Total cost: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count()
        << " s" << std::endl;

    return repairedCount;
}

// 通用修复流水线调用入口：通过 enabled 开关选择模块，addStage 顺序决定执行顺序。
int runShapeRepairPipeline(const std::string& fileName) {
    const std::string inputFilePath = "./assets/" + fileName + ".stp";
    const std::filesystem::path outputDirectory = "./assets/fixedRepairPipeline";
    std::filesystem::create_directories(outputDirectory);
    const std::string outputFilePath =
        (outputDirectory / ("RepairPipeline_" + fileName + ".stp")).string();

    ApplicationInit appInit;
    const auto startTime = std::chrono::system_clock::now();
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file failed: " << inputFilePath << std::endl;
        return -1;
    }

    // 调整这些开关即可指定参与本次修复的模块。
    const bool enableHealNarrowGap = true;
    const bool enableHealClosedNarrowGap = true;
    // 自动补孔需要结合模型单位设置尺寸范围；当前已启用。
    const bool enableHealSimpleHole = true;

    ShapeRepairPipeline::Params pipelineParams;
    pipelineParams.rejectShapeTypeChange = true;
    pipelineParams.validateAfterEachStage = true;
    pipelineParams.validateFinalResult = true;
    ShapeRepairPipeline pipeline(pipelineParams);

    HealNarrowGap::Params narrowGapParams;
    HealClosedNarrowGap::Params closedNarrowGapParams;
    HealSimpleHole::Params simpleHoleParams;

    // 性能对照：关闭 HealNarrowGap 的 crossed-face 两级兜底策略。
    narrowGapParams.enableCrossedFaceBoundaryRepair = false;

    // 与独立入口共用 HealSimpleHole::Params 中的默认配置。

	closedNarrowGapParams.midlineConstructionMode = HealClosedNarrowGap::MidlineConstructionMode::ProjectionCorrespondence;

    // 添加顺序就是执行顺序：先补孔，避免细缝修复改变孔洞自由边界的拓扑关系。
    pipeline.addStage(ShapeRepairStages::MakeHealSimpleHoleStage(
        simpleHoleParams,
        enableHealSimpleHole));
    pipeline.addStage(ShapeRepairStages::MakeHealNarrowGapStage(
        narrowGapParams,
        enableHealNarrowGap));
    pipeline.addStage(ShapeRepairStages::MakeHealClosedNarrowGapStage(
        closedNarrowGapParams,
        enableHealClosedNarrowGap));

    std::cout << "Running configured shape repair pipeline......" << std::endl;
    const int repairedCount = appInit.runRepairPipelineOnLowestContainersAndExport(
        outputFilePath,
        pipeline);
    if (repairedCount < 0) {
        std::cerr << "STEP file write failed: " << outputFilePath << std::endl;
        return -1;
    }

    std::cout << "Repair pipeline committed operation count: "
              << repairedCount << std::endl;
    std::cout << "STEP file write successfully to " << outputFilePath << std::endl;

    const auto endTime = std::chrono::system_clock::now();
    std::cout << "Total cost: "
              << std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count()
              << " s" << std::endl;
    return repairedCount;
}

int main()
{
    std::string fileName = "h6k";
    //modelConvert(fileName);
	//mergeShortEdge(fileName);
    //healNarrowGap(fileName);
    //healClosedGap(fileName);
	//healSimpleHole(fileName);
	//sewShellsOnStepShape(filename);
	runShapeRepairPipeline(fileName);
    return 0;
}
