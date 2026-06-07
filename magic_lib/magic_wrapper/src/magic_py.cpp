//
// Created by max on 2021-07-05.
//

#include <pybind11/pybind11.h>
namespace py=pybind11;

#include "GrayImage.h"
#include "PyRAG.h"
#include "PyGMM.h"
#include "GradientWrapper.h"
#include "WatershedWrapper.h"
//#include "permuto_wrapper.h"
#include "R2ProductInfo.h"
#include "TextureWrapper.h"
// #include "ImageConversion.h"

// int get_randmax(){
//     return RAND_MAX;
// }


PYBIND11_MODULE(magic_py, m) {

    m.doc() = "MAGIC_PY python module";

    // standalone functions included in the magic_lib module
    m.def("Gradient", &Gradient, "Compute the gradient magnitude of an image");
    m.def("Watershed", &Watershed, "Compute the watershed transform from the gradient magnitude of an image");
    m.def("BlockAverage", &BlockAverage, "block average");
    // m.def("permuto", &permuto_gauss_transform, "permutohedral lattice gauss transform!!!");
    // m.def("get_randmax", &get_randmax, "get the value of RAND_MAX for this system");
    m.def("extract_glcitr", &extract_glcitr, "extract GLCM features");
    m.def("test_array_capsule", &test_array_capsule, "memory test");

    // reads product.xml and returns relevant data in a python dictionary
    // currently only returns pass direction and ground control points because that's all
    // I need right now and I am lazy. But any other relevant data can easily be added to
    // the returned dictionary.
    m.def("ReadR2Product", &GetProductInfo, py::arg("sProductPath")=" ", py::arg("blockaverage")=1, "extract information from a radarsat 2 product.xml file!!!");
    m.def("ReadR2Lut", &GetCalibrationLut,  "get the calibration lut from a radarsat-2 calibration file!!!");

    // the basic region adjacency graph where edges represent pairwise interactions
    py::class_<PyRAG>(m, "PyRAG")
        .def(py::init<py::array_t<float>, py::array_t<float>, py::array_t<int>, py::array_t<int>, int>())

        // getters for graph structure properties
        .def("GetEdges", &PyRAG::GetEdges)
        .def("GetVertexCount", &PyRAG::GetNodeCount)
        .def("GetEdgeCount", &PyRAG::GetEdgeCount)
        .def("GetOriginalWatershed", &PyRAG::GetOriginalWatershed)
        .def("GetWatershed", &PyRAG::GetWatershed)
        .def("GetChainCodes", &PyRAG::GetChainCodes)
        .def("GetBoundaryLengths", &PyRAG::GetBoundaryLengths)
        .def("GetBoundaryMap", &PyRAG::GetBoundaryMap)

        // getters for region attributes
        .def("GetRegionClassLabels", &PyRAG::GetRegionClassLabels)
        .def("GetRegionCentroids", &PyRAG::GetRegionCentroids)
        .def("GetRegionMeans", &PyRAG::GetRegionMeans)
        .def("GetRegionCoMeans", &PyRAG::GetRegionCoMeans)
        .def("GetRegionPixelCounts", &PyRAG::GetRegionPixelCounts)

        // getters for classification related attributes
        .def("GetConnectedClassMap", &PyRAG::GetConnectedClassMap)
        .def("GetConnectedClassMapGreedy", &PyRAG::GetConnectedClassMapGreedy)
        .def("GetConnectedClassMapBoundary", &PyRAG::GetConnectedClassMapBoundary)
        .def("GetMergingLUT", &PyRAG::GetMergingLUT)
        .def("GetMergingLUTFromMaps", &PyRAG::GetMergingLUTFromMaps)
        .def("GetGMMParams", &PyRAG::GetGMMParams)
        .def("GetGMMUnary", &PyRAG::GetGMMUnary)
        .def("GetEdgePenalty", &PyRAG::GetEdgePenalty)
        .def("GetBeta", &PyRAG::GetBeta)
        .def("GetBoundaryPopulationRatio", &PyRAG::GetBoundaryPopulationRatio)

        // methods related to updating classification parameters
        .def("UpdateGMMParam", &PyRAG::UpdateGMMParam)
        .def("UpdateBoundaryParam", &PyRAG::UpdateBoundaryParam)
        .def("UpdateBoundaryParamFisher", &PyRAG::UpdateBoundaryParamFisher)
        .def("SetCurrentIteration", &PyRAG::SetCurrentIteration)

        // graph operations
        .def("IRGSStep", &PyRAG::IRGSStep)
        .def("KmeansInitialization", &PyRAG::KmeansInitialization)
        .def("SetClassLabelsByIndex", &PyRAG::SetClassLabelsByIndex)
        .def("MergeRegions", &PyRAG::MergeRegions)
        ;

//    // the curvature graph / partial enumeration graph, where vertices represent clusters / cliques / patches of regions
//    // I really need to get my terminology sorted out
//    py::class_<PyCurvatureGraph>(m, "PyCurvatureGraph")
//        .def(py::init<py::array_t<float>, py::array_t<float>, py::array_t<int>, py::array_t<int>, int>())
//
//        // getters for graph structure properties
//        .def("GetEdges", &PyCurvatureGraph::GetEdges)
//        .def("GetVertexCount", &PyCurvatureGraph::GetNodeCount)
//        .def("GetEdgeCount", &PyCurvatureGraph::GetEdgeCount)
//        .def("GetOriginalWatershed", &PyCurvatureGraph::GetOriginalWatershed)
//        .def("GetWatershed", &PyCurvatureGraph::GetWatershed)
//        .def("GetChainCodes", &PyCurvatureGraph::GetChainCodes)
//        .def("GetBoundaryLengths", &PyCurvatureGraph::GetBoundaryLengths)
//        .def("GetBoundaryMap", &PyCurvatureGraph::GetBoundaryMap)
//
//        // getters for region attributes
//        .def("GetRegionClassLabels", &PyCurvatureGraph::GetRegionClassLabels)
//        .def("GetRegionCentroids", &PyCurvatureGraph::GetRegionCentroids)
//        .def("GetRegionMeans", &PyCurvatureGraph::GetRegionMeans)
//        .def("GetRegionCoMeans", &PyCurvatureGraph::GetRegionCoMeans)
//        .def("GetRegionPixelCounts", &PyCurvatureGraph::GetRegionPixelCounts)
//
//        // getters for classification related attributes
//        .def("GetConnectedClassMap", &PyCurvatureGraph::GetConnectedClassMap)
//        .def("GetMergingLUT", &PyCurvatureGraph::GetMergingLUT)
//        .def("GetMergingLUTFromMaps", &PyCurvatureGraph::GetMergingLUTFromMaps)
//        .def("GetGMMParams", &PyCurvatureGraph::GetGMMParams)
//        .def("GetGMMUnary", &PyCurvatureGraph::GetGMMUnary)
//        .def("GetEdgePenalty", &PyCurvatureGraph::GetEdgePenalty)
//        .def("GetBeta", &PyCurvatureGraph::GetBeta)
//        .def("GetBoundaryPopulationRatio", &PyCurvatureGraph::GetBoundaryPopulationRatio)
//
//        // methods related to updating classification parameters
//        .def("UpdateGMMParam", &PyCurvatureGraph::UpdateGMMParam)
//        .def("UpdateBoundaryParam", &PyCurvatureGraph::UpdateBoundaryParam)
//        .def("UpdateBoundaryParamFisher", &PyCurvatureGraph::UpdateBoundaryParamFisher)
//        .def("SetCurrentIteration", &PyCurvatureGraph::SetCurrentIteration)
//
//        // graph operations
//        //.def("IRGSStep", &PyCurvatureGraph::IRGSStep)
//        .def("KmeansInitialization", &PyCurvatureGraph::KmeansInitialization)
//        .def("SetClassLabelsByIndex", &PyCurvatureGraph::SetClassLabelsByIndex)
//        //.def("MergeRegions", &PyCurvatureGraph::MergeRegions) // I haven't adapted this method for PyCurvatureGraph objects yet
//
//        // extra methods for PyCurvatureGraph that don't exist for the regular PyRAG
//        .def("GetClusterCount", &PyCurvatureGraph::GetClusterCount)
//        .def("GetClusterOrdering", &PyCurvatureGraph::GetClusterOrdering)
//        .def("GetClusterEdgeCount", &PyCurvatureGraph::GetClusterEdgeCount)
//        .def("GetClusterCenterMap", &PyCurvatureGraph::GetClusterCenterMap)
//        .def("GetClusterBoundaries", &PyCurvatureGraph::GetClusterBoundaries)
//        .def("GetClusterEdges", &PyCurvatureGraph::GetClusterEdges)
//        .def("GetClusterRegions", &PyCurvatureGraph::GetClusterRegions)
//        .def("GetClusterDegrees", &PyCurvatureGraph::GetClusterDegrees)
//        .def("GetClusterSizes", &PyCurvatureGraph::GetClusterSizes)
//        .def("GetInnerBoundaryCounts", &PyCurvatureGraph::GetInnerBoundaryCounts)
//        .def("GetInnerBoundaryDepartureDirections", &PyCurvatureGraph::GetInnerBoundaryDepartureDirections) // isn't that a mouthful
//        .def("GetInnerBoundaryDisplacements", &PyCurvatureGraph::GetInnerBoundaryDisplacements)
//        .def("GetCliqueEdgeSharedRegions", &PyCurvatureGraph::GetCliqueEdgeSharedRegions)
//        .def("GetCutVertices", &PyCurvatureGraph::GetCutVertices)
//        .def("GetBoundaryPairChainCodes", &PyCurvatureGraph::GetBoundaryPairChainCodes)
//        .def("GetBoundaryPairChainCode", &PyCurvatureGraph::GetBoundaryPairChainCode)
//        .def("GetBoundaryChainCode", &PyCurvatureGraph::GetBoundaryChainCode)
//        .def("fTest", &PyCurvatureGraph::fTest)
//
//        ;

    /// Wrapper class for gaussian mixture model parameters, making it easy to pass them between python and c++
    py::class_<PyGMM>(m, "PyGMM")

            // constructor from numpy arrays passed positionally
            .def(py::init<int, int, py::array_t<float>, py::array_t<float>, py::array_t<float>, py::array_t<float>, py::array_t<float>, py::array_t<float>, py::array_t<float>>())

            // constructor from dict of GMM params. No error checking so be careful lol
            .def(py::init<py::dict>())

            .def("GetMu", &PyGMM::GetMu)
            .def("GetCovInv", &PyGMM::GetCovInv)
            .def("GetPrior", &PyGMM::GetPrior)
            .def("GetLnSigma", &PyGMM::GetLnSigma)
            .def("GetCovDet", &PyGMM::GetCovDet)
            .def("GetLnPrior", &PyGMM::GetLnPrior)

            .def("_dict", &PyGMM::todict);

}
