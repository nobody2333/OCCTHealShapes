#include <topo/Body.h>

#include <vector>
#include <format>
#include <chrono>
#include <iostream>
#include <TopExp_Explorer.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopAbs_ShapeEnum.hxx>

#include "ApplicationInit.h"
#include "GeometryUtils.h"
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
    std::string outputFilePath = "./assets/fixedShortEdge/fixed" + fileName + ".stp";
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



int main()
{
    std::string fileName = "dumper truck";
    //modelConvert(fileName);
	mergeShortEdge(fileName);
    return 0;
}
