#include <topo/Body.h>

#include <vector>
#include <format>
#include <chrono>

#include "ApplicationInit.h"
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

void mergeShortEdge(const std::string& fileName) {
    std::string inputFilePath = "./assets/" + fileName + ".stp";

    std::string outputFilePath = "./assets/fixedShortEdge/fixed" + fileName + ".stp";

    ApplicationInit appInit;
    if (!appInit.importStep(inputFilePath)) {
        std::cerr << "Reading file fasiled" << inputFilePath << std::endl;
        return;
    }

    std::vector<TopoDS_Shape> inputShapes = appInit.getAllShapes();
    if (inputShapes.empty()) {
        std::cerr << "Cannot find valid shapes!" << std::endl;
        return;
    }

    MergeShortEdge shortEdgeMerger(0.005, 170.0);
    shortEdgeMerger.setInputShapes(inputShapes);

    std::cout << "Merging short edges......" << std::endl;
    int mergedCount = shortEdgeMerger.performMerge();

    // 获取合并后的最终结果
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
    std::string fileName = "dlr-f4";
    //modelConvert(fileName);
	mergeShortEdge(fileName);

    return 0;
}
