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
#include "MergeShortEdge.h"
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
    size_t pos = fileName.find('.');
    std::string suffix = fileName.substr(pos + 1);
    std::string modleName = fileName.substr(0, pos);

    std::string inputFilePath = "./assets/" + fileName;
    std::string outputFilePath = "./assets/fixedShortEdge/fixed" + modleName + ".stp";
    ApplicationInit appInit;

    if (suffix == "iges") {
        if (!appInit.importIges(inputFilePath)) {
            std::cerr << "Reading file fasiled" << inputFilePath << std::endl;
            return;
        }
    }
    else {
        if (!appInit.importStep(inputFilePath)) {
            std::cerr << "Reading file fasiled" << inputFilePath << std::endl;
            return;
        }
    }
    std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return;
    }

    std::cout << "Top-level shape count: " << inputShapes.size() << std::endl;
    for (std::size_t i = 0; i < inputShapes.size(); ++i) {
        std::cout << std::format("  Shape[{}] type: {}", i, static_cast<int>(inputShapes[i].ShapeType())) << std::endl;
    }


    MergeShortEdge shortEdgeMerger;
    shortEdgeMerger.setInputShapes(inputShapes);
    shortEdgeMerger.setLengthThreshold(0.02);
	shortEdgeMerger.setAngleThresholdDeg(120.0);
    std::cout << std::format("Merging short edges with length threshold: {} and angle threshold: {} degrees.", shortEdgeMerger.getLengthThreshold(), shortEdgeMerger.getAngleThresholdDeg()) << std::endl;


    std::cout << "Merging short edges......" << std::endl;
    int mergedCount = shortEdgeMerger.performMerge();

    std::cout << std::format("Merged {} edge pairs.", mergedCount) << std::endl;
    TopoDS_Shape resultShape = shortEdgeMerger.getResultShape();

    if (appInit.exportStep(resultShape, outputFilePath)) {
        std::cout << "STEP file write successfully to " << outputFilePath << std::endl;
    }
    else {
        std::cerr << "Error: Failed to write STEP file." << std::endl;
        return;
    }

    return;
}



int main()
{
    std::string fileName = "bentDihedral.stp";
    //modelConvert(fileName);
	mergeShortEdge(fileName);
    return 0;
}
