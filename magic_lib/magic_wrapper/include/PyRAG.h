//
// Created by max on 2021-12-23.
//

#ifndef CROWS_PYRAG_H
#define CROWS_PYRAG_H

#include "RegionAdjacencyGraph.h"
#include "RegionMerger.h"
#include "ImageConversion.h"
#include "Watershed.h"
#include "DataType.h"
#include "EstimateUtil.h"
#include "RegionBasedClustering.h"
#include "PyGMM.h"
#include "StrengthPenaltyBoundary.h"
#include "Fisher.h"

class PyRAG{
public:
    PyRAG(py::array_t<float> srcImage, py::array_t<float> Gradient, py::array_t<int> Watershed, py::array_t<int> mask, int nClassCnt, int nBoundaryMode=0);
    ~PyRAG();

protected:

    int m_nClassCnt;
    int m_nFeatureCnt;
    int m_nBoundaryMode;

    // Image data
    TGrayImage<float>** m_ppSrcImg;
    TGrayImage<float>* m_pGradient;
    TMonoImage* m_pMask;

    // map region pixels -> region labels
    TGrayImage<int>* m_pWatershedMap;

    // map boundary pixels -> boundary labels
    TGrayImage<int>* m_pBoundaryLabelMap;

    // the region adjacency graph
    TRegionMerger* m_pGraph;

    // Mapping from the (possibly nonconsecutive) region labels to a consecutive labelling.
    // If you don't do any region merging, this is just the sequence (0,1,...,N_regions-1).
    // If you merge regions, the sequence becomes nonconsecutive (some regions are missing) so this mapping
    // can be used to keep track of which regions are still present.
    std::vector<int> m_vConsecutiveMapping;

    // some parameters related to classification / IRGS
    GMM* gmm;
    bool m_bIRGSsetup; // flag to indicate whether the classifier has been setup for running IRGS.
    int m_nIRGSIterElapsed;
    TRegionClassifier* m_pClassifier;

public:

    // getters for various graph structure attributes
    py::array_t<int> GetEdges() const;
    int GetNodeCount() const { return m_pGraph->GetVertexCnt();}
    int GetEdgeCount() const { return m_pGraph->GetEdgeCnt();}
    py::array_t<int> GetOriginalWatershed() const;
    std::tuple<py::array_t<int>, py::array_t<int>> GetWatershed() const;
    std::tuple<py::array_t<int>, py::array_t<int>, py::array_t<int>, py::array_t<int>> GetChainCodes() const;
    py::array_t<int> GetBoundaryMap();
    py::array_t<int> GetBoundaryLengths();

    // getters for various region attributes
    py::array_t<int> GetRegionClassLabels() const;
    py::array_t<float> GetRegionCentroids() const;
    py::array_t<float> GetRegionMeans() const;
    py::array_t<float> GetRegionCoMeans() const;
    py::array_t<int> GetRegionPixelCounts() const;

    // getters for classification-related attributes
    py::array_t<int> GetConnectedClassMap();
    py::array_t<int> GetConnectedClassMapGreedy();
    py::array_t<int> GetConnectedClassMapBoundary();
    bool IsBoundary(TGrayImage<int>* pMap, TMonoImage* pMask, int nRow, int nCol, int*& nClassLabel, bool& bIndepBoundary);
    int GetClassCnt() const {return m_nClassCnt;}
    py::array_t<int> GetMergingLUT();
    py::array_t<int> GetMergingLUTFromMaps(py::array_t<int> aMapOld, py::array_t<int> aMapNew);
    PyGMM GetGMMParams() const;
    py::array_t<float> GetGMMUnary() const;
    py::array_t<float> GetEdgePenalty() const;
    double GetBeta() {return m_pGraph->GetCurBeta();}
    double GetBoundaryPopulationRatio() {return m_pGraph->GetBoundaryPopulationRatio(); }

    // methods for updating classification-related parameters
    bool UpdateGMMParam();
    bool UpdateBoundaryParamFisher(double dbK, double dbBeta1, double dbBeta2, int it);
    bool UpdateBoundaryParam(double dbK, int it);
    void SetupIRGS(); // set up the graph to do IRGS iterations
    void SetCurrentIteration(int it) { m_nIRGSIterElapsed = it; }

    // Operations on the graph (classification, region merging, setting class labels etc)
    void IRGSStep(double dbK, double dbBeta1, double dbBeta2, bool bVerbose); // do one IRGS iteration.
    py::array_t<int> KmeansInitialization();
    void SetClassLabelsByIndex(py::array_t<int> ClassLabels, py::array_t<int> NodeIndices);
    void PreprocessWatershed() const;
    py::array_t<int> MergeRegions();

private:
    // Unary / Pairwise functions copied from the GaussianRegionClustering for convenience
    double GetUnaryEnergy(TRegion *pRegion, int nClassIndex) const;
    double GetBinaryEnergy(int nLabel1, int nLabel2, TBoundary *pBoundary) const;


};

#endif //CROWS_PYRAG_H
