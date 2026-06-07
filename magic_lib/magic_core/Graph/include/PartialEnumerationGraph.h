//
// Created by max on 2021-12-14.
//

#include "RegionAdjacencyGraph.h"
#include "Graph.h"
#include "EstimateUtil.h"

#ifndef ICETOOLS_PARTIALENUMERATIONGRAPH_H
#define ICETOOLS_PARTIALENUMERATIONGRAPH_H

///
/// CLUSTER NODE
///
class TClusterNode: public TVertex{
public:
    TClusterNode(TGraph* pGraph, int Label);
    ~TClusterNode();

    // methods for cluster construction
    void AddCenterPoint(POINTex Pt) {m_vCenterPoints.push_back(Pt); m_nCenterPoints++; }

public:

    int m_nMemberRegions;
    std::vector<int> m_vMemberRegions;

    int m_nCenterPoints;
    std::vector<POINTex> m_vCenterPoints;

    int m_nMemberBoundaries;
    std::vector<int> m_vMemberBoundaries;
    std::vector<int> m_vMemberBoundaryDirs;
    std::vector<int> m_vMemberBoundaryHead;
    std::vector<std::pair<int, int>> m_vMemberBoundaryDisplacement;
};

///
/// CLUSTER EDGE
///
class TClusterEdge: public TEdge, public TAVLTreeNode {
public:
    TClusterEdge(TGraph* pGraph, int nLabel1, int nLabel2, int nBoundaryIndex);
    ~TClusterEdge();
    void AddSharedRegion(int SharedRegionIdx);
    int m_nBoundaryIndex;
    std::vector<int> m_vSharedRegions;
};


///
/// PARTIAL ENUMERATION GRAPH
///
class TPartialEnumerationGraph: public TGraph{
public:

    TPartialEnumerationGraph(TRegionAdjacencyGraph* pRAG,
                             TGrayImage<int>* pWatershedMap,
                             TGrayImage<int>* pBoundaryLabelMap,
                             TGrayImage<int>* pClusterCenterMap,
                             TMonoImage* pMask);
    ~TPartialEnumerationGraph();


    // constructs the cluster graph based on the watershed map and the cluster center map
    void ConstructGraph();
    void LoadClusterCenters();

    std::vector<int> m_vRAGCutVertices;

    std::vector<int> m_vChainCodes;
    std::vector<int> m_vChainCodeLengths;
    std::vector<int> m_vChainCodeOffsets;
    std::vector<POINTex> m_vStartPoints;
    std::vector<POINTex> m_vEndPoints;

    void AppendBoundaryPairChainCode(std::vector<int> & codes, std::vector<int> & pair_idx, std::vector<int> & clique_idx, TClusterNode* pNode, int b1, int b2);

    std::vector<int> GetBoundaryPairChainCode(TClusterNode* pNode, int b1, int b2);
    std::vector<int> GetInnerBoundaryChainCode(TClusterNode* pNode, int b);

    // for debugging, do not use this please if you can avoid it.
    TClusterNode* GetCluster(int nLabel){ return (TClusterNode*) m_pClusterLUT[nLabel];}

protected:
    TRegionAdjacencyGraph* m_pRAG;
    TGrayImage<int>* m_pWatershedMap;
    TGrayImage<int>* m_pBoundaryLabelMap;
    TGrayImage<int>* m_pClusterCenterMap;
    TMonoImage* m_pMask;

    TClusterNode** m_pClusterLUT;

    // functions for graph construction
    TClusterEdge* IsEdgeInList(int nLabel1, int nLabel2);
    TAVLTree m_EdgeList;
    TVertex* FirstNeighbor(TVertex* pVertex) override;
    TVertex* NextNeighbor(TVertex* pVertex) override;
    virtual TEdge* ConstructEdge(TVertex* pNeighbor1, TVertex* pNeighbor2);

};

#endif //ICETOOLS_PARTIALENUMERATIONGRAPH_H
