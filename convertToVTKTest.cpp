#include <topo/Body.h>

#include <vector>
#include <format>
#include <chrono>


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
    std::string filePath = "./assets/shortEdge/" + fileName + ".stp";

}



int main()
{
    std::string fileName = "F007-CarAircon-model3";
    //modelConvert(fileName);

    return 0;
}
