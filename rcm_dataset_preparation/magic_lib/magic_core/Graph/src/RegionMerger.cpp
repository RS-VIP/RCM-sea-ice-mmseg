// RegionMerger.cpp : implementation of the TRegionMerger class
//
// Merge similar regions
/////////////////////////////////////////////////////////////////////////////

#include "RegionMerger.h"
#include "EdgeWeightKey.h"
#include "EstimateUtil.h"
#include "RegionClustering.h"
#include "StrengthPenaltyBoundary.h"
#define BETA_THRESHOLD 1.2
#define BETA_DECAY 8

//
// Constructor (region based)
// 
TRegionMerger::TRegionMerger(TRegionClassifier* pClassifier,
	TGrayImage<int>* pMap,
	TMonoImage* pMask, 
	TGrayImage<FLOAT>** ppSrcImg,
	int nFeatureCnt, TRegionAdjacencyGraph::statistics_mode enStatMode,
	TGrayImage<FLOAT>* pGrad,
    TRegionAdjacencyGraph::boundary_mode enBoundMode)
	: TGIEPGraph(pClassifier, pMap, pMask, ppSrcImg, nFeatureCnt, enStatMode, pGrad, enBoundMode)
{
}

//
// Destructor
// 
TRegionMerger::~TRegionMerger(){}

//
// Merge nearest vertice
// 
BOOL TRegionMerger::Merge()
{
	// Find the most similar two vertice and merge them
	BOOL bChanged = false;
	m_pNearestEdge = FindNearest();
	int iTotalMergedSofar = 0;
	int loopcount = 0;

	while (!ShouldStop())
	{
		//printf("%d\n", m_nVertexCnt);
		if ( m_pNearestEdge->m_pNeighbor1->m_nNeighborCnt >=  m_pNearestEdge->m_pNeighbor2->m_nNeighborCnt )
			TGraph::Merge(m_pNearestEdge->m_pNeighbor1, m_pNearestEdge->m_pNeighbor2);
		else
			TGraph::Merge(m_pNearestEdge->m_pNeighbor2, m_pNearestEdge->m_pNeighbor1);

		loopcount++;
		//if(CanSaveInterResult())
		//	this->SaveIntermediateResult(this->m_nIter+1, loopcount,this->m_nVertexCnt,-2);

		bChanged = true;
		m_pNearestEdge = FindNearest(); //QIN_MOD
	}

	return bChanged;
}

std::vector<std::pair<int,int>> TRegionMerger::MergeAndReturnIndex()
{
    std::vector<std::pair<int,int>> MergedIndices;

    // Find the most similar two vertices and merge them
    m_pNearestEdge = FindNearest();

    while (!ShouldStop())
    {
        if ( m_pNearestEdge->m_pNeighbor1->m_nNeighborCnt >=  m_pNearestEdge->m_pNeighbor2->m_nNeighborCnt ){
            MergedIndices.push_back(std::make_pair(m_pNearestEdge->m_pNeighbor1->m_nLabel, m_pNearestEdge->m_pNeighbor2->m_nLabel));
            TGraph::Merge(m_pNearestEdge->m_pNeighbor1, m_pNearestEdge->m_pNeighbor2);
        }
        else {
            MergedIndices.push_back(std::make_pair(m_pNearestEdge->m_pNeighbor2->m_nLabel, m_pNearestEdge->m_pNeighbor1->m_nLabel));
            TGraph::Merge(m_pNearestEdge->m_pNeighbor2, m_pNearestEdge->m_pNeighbor1);
        }

        m_pNearestEdge = FindNearest(); //QIN_MOD
    }

    return MergedIndices;
}

//
// Find the edge that corresponds to the two nearest vertice
// 
// Return
//  NULL: no more vertex; otherwise: the corresponding edge
// 
TEdge* TRegionMerger::FindNearest()
{
	return m_pFirstEdge;
}

//
// Determine if the merging should stop or not
// 
// return
//  true: stop; false: continue
// 
BOOL TRegionMerger::ShouldStop()
{
	return (!(m_pNearestEdge && ((TStrengthPenaltyBoundary*) m_pNearestEdge)->GetEdgeWeight() < 0));
}


//
// Iterative classification and region merging
// 
// Param in
//  nIterations: number of iterations
// 
// return
//  the segmentation map
//
TImage* TRegionMerger::Merge(double& dbFinalBeta, TRegionClassifier*& pRetClassifier, int nBetaMode, bool bOutIntermResult, int nIterations, double dbK, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2)
{
	TImage *pRet;

	// Construct the graph
	ConstructGraph();

	// Initialize the classifier
	if (!m_bNullClassifier)
		m_pClassifier->Init();

	if ( nBetaMode == FIXED )
		m_dbBeta = dbBeta1;

	double dbFisher;

	int nUnchangedCnt = 0;


	for ( m_nIter = 0; m_nIter < nIterations; m_nIter++ )
	{

		// Update the penalty function parameter K
		UpdateK(dbK, m_nIter);

		if ( ((TRegionClustering*) m_pClassifier)->GetType() & TRegionClustering::VAFM )
			m_dbAlpha = UpdateAlpha(m_nIter, dbAlpha1, dbAlpha2);

		if ( nBetaMode == MLFIXED ){
			m_dbBeta = UpdateBeta(GetBoundaryPopulationRatio());
			m_dbBeta *= dbBeta1;
		}

		if ( nBetaMode == MLFISH ){
            //std::cout << "------------------------------------------" << std::endl;
            //std::cout << "Iteration: " << m_nIter << std::endl;
			m_dbBeta = UpdateBeta(GetBoundaryPopulationRatio());
            //std::cout << "beta_0 : " << m_dbBeta << std::endl;
			dbFisher = m_pClassifier->ComputeMinFisher();
            //std::cout << "Fisher: " << dbFisher << std::endl;
			m_dbBeta *= ( dbBeta1 * dbFisher / dbBeta2 / (1 + dbFisher / dbBeta2) );
            //std::cout << "adjusted beta : " << m_dbBeta << std::endl;


		}

		// Compute the energy corresponding to edge information
		UpdatePenalty();

		// Classification
		if (!m_bNullClassifier)
			m_pClassifier->IntermediateClassify();

		SortEdge();

		// Merging
		if (!Merge()){
			nUnchangedCnt ++;
		}
		else{
			nUnchangedCnt = 0;
		}

		if (!m_bNullClassifier){
			if (!m_pClassifier->UpdateParam())
				goto END;
		}
	}

END:
	// Final classification
	if (!m_bNullClassifier)
		m_pClassifier->FinalUpdateParam();

	pRetClassifier = m_pClassifier;

	if (!(pRet = m_pClassifier->GetClassMap())){
		TEstimateUtil* pEstUtil = new TEstimateUtil();
		pRet = new TGrayImage<int>(m_pMap);
		pEstUtil->RelabelMap( (TGrayImage<int>*) pRet, m_pMask );
		delete pEstUtil;
	}

	dbFinalBeta = m_dbBeta;
	return pRet;
}

std::tuple<double, double, bool> TRegionMerger::MergingStep(double dbK, double dbBeta1, double dbBeta2, int nIt) {
    double dbFisher;

    // Update the penalty function parameter K
    UpdateK(dbK, nIt);

    // Update Beta
    m_dbBeta = UpdateBeta(GetBoundaryPopulationRatio());
    dbFisher = m_pClassifier->ComputeMinFisher();
    m_dbBeta *= ( dbBeta1 * dbFisher / dbBeta2 / (1 + dbFisher / dbBeta2) );

    // Compute the energy corresponding to edge information
    UpdatePenalty();

    // Classification
    bool success = m_pClassifier->IntermediateClassify();

    SortEdge();
    bool changed =  Merge();

    m_pClassifier->UpdateParam();

    return std::make_tuple(m_dbBeta, dbFisher, changed);
}

// 
// Update the parameters associated with the specified vertex
// 
void TRegionMerger::Update(TVertex* p){
	TStrengthPenaltyBoundary* pEdge;
	double dbVal,neighborVertexArea,curVertexArea;
	TEdge *pCur;
	TVertex* neighborVertex;
	for (int i = 0; i < p->m_nNeighborCnt; i ++){
		pEdge = (TStrengthPenaltyBoundary*) p->m_pEdges[i];

		// Since the parameters of vertex have changed, need to recompute. Here,
		// only need to recompute the region energy difference

		//Debug by Zhijie
		//When the neighbor region is too small, we don't do the recomputation to make the method fast
		neighborVertex = p->GetNeighborVertex(pEdge);
		neighborVertexArea = (double)(((TGaussianRegion*)neighborVertex)->m_nPixelCnt);
		curVertexArea = (double)(((TGaussianRegion*)p)->m_nPixelCnt);
		if (neighborVertexArea/curVertexArea <= 0.005)
			continue;
				
		pEdge->ComputeRegionDiffEnergy();

		if (m_nEdgeCnt <= 1)  // only one or null elements in the list
			continue;

		// Detach from the linked list
		if (pEdge->m_pPrev)
			pEdge->m_pPrev->m_pNext = pEdge->m_pNext;
		else
			m_pFirstEdge = pEdge->m_pNext;
		if (pEdge->m_pNext)
			pEdge->m_pNext->m_pPrev = pEdge->m_pPrev;
		else
			m_pLastEdge = pEdge->m_pPrev;
		pEdge->m_pPrev = 0;
		pEdge->m_pNext = 0;

		// Insert back into the linked list in the order of its weight (merging
		// criterion)
		if ((dbVal = pEdge->GetEdgeWeight()) >= 0){
			// No need to sort, as it will not be considered in merging. Directly
			// insert into the end of the linked list
			if (m_pLastEdge)
				m_pLastEdge->m_pNext = pEdge;
			pEdge->m_pPrev = m_pLastEdge;
			m_pLastEdge = pEdge;
		}else{
			// Find the position and insert
			if (dbVal <=  m_pFirstEdge->GetEdgeWeight()){
				// Put as the first in the list
				pEdge->m_pNext = m_pFirstEdge;
				m_pFirstEdge->m_pPrev = pEdge;
				m_pFirstEdge = pEdge;
			}else{
				pCur = m_pFirstEdge;
				while (pCur->m_pNext){
					if (dbVal <= (pCur->m_pNext)->GetEdgeWeight())
						break;
					pCur = pCur->m_pNext;
				}
				pEdge->m_pNext = pCur->m_pNext;
				if (pEdge->m_pNext)
					pEdge->m_pNext->m_pPrev = pEdge;
				else
					m_pLastEdge = pEdge;
				pCur->m_pNext = pEdge;
				pEdge->m_pPrev = pCur;
			}
		}
	}
}


// Sort the edge by its weight (merging criterion)
//
void TRegionMerger::SortEdge()                               
{
	//TGIEPGraph::UpdatePenalty();	//	compute edge penalty value
	// Sort the edge by its weight (merging criterion) to speed up the
	// searching for the lowest
	TBoundary* pCur = (TBoundary*) m_pFirstEdge;
	TBoundary* pNext;
	double dbVal;
	m_BoundaryList.Reset();  // clear the graph
	m_pFirstEdge = 0;  // empty the linked list for space of result
	m_pLastEdge = 0;
	while (pCur){
		dbVal = pCur->GetEdgeWeight();	//	combine regional statistics and edge penalty value
		pNext = (TBoundary*) pCur->m_pNext;
		if (dbVal >= 0){
			// No need to sort, as it will not be considered in merging. Directly
			// insert it into a linked list
			pCur->m_pNext = m_pFirstEdge;
			if (m_pFirstEdge)
				m_pFirstEdge->m_pPrev = pCur;
			else
				m_pLastEdge = pCur;  // mark current as the last one in the linked list
			m_pFirstEdge = pCur;
		}else{
			pCur->Reset(new TEdgeWeightKey(dbVal));
			m_BoundaryList.Insert(pCur);
		}
		pCur = pNext;
	}

	// Put the AVL tree nodes into a linked list and concatenate to the existing
	// positive weight node linked list
	TBoundary *pFirst = 0;
	if (m_BoundaryList.m_pRoot)
		pFirst = (TBoundary*) m_BoundaryList.GetMostLeftChild(m_BoundaryList.m_pRoot);
	pCur = pFirst;
	while (pCur){
		pNext = (TBoundary*) m_BoundaryList.GetNextNode(pCur);
		if (pNext){
			// Link the two consecutive nodes together
			pCur->m_pNext = pNext;
			pNext->m_pPrev = pCur;
			pCur = pNext;
		}else{
			// Link with the existing positive weight node linked list
			pCur->m_pNext = m_pFirstEdge;
			if (m_pFirstEdge)
				m_pFirstEdge->m_pPrev = pCur;
			break;
		}
	}
	if (pFirst)
		m_pFirstEdge = pFirst;
	if (!m_pLastEdge)
		m_pLastEdge = pCur;

	try{
		m_pFirstEdge->m_pPrev = 0;
	}
	catch(...){}

}

// 
// Get the map in which each color represents a connected class
// region of the same class
// 
TGrayImage<int>* TRegionMerger::GetConnectedClassRegionMap()
{
	// Get maximum number of remaining segments
	TEstimateUtil EstUtil;
	int nMaxCnt = EstUtil.GetMaxLabel(m_pMap, m_pMask) + 1;

	// Generate the label mapping
	int* nLabels = new int[nMaxCnt];
	for (int i = 0; i < nMaxCnt; i ++)
		nLabels[i] = -1;
	int nNewLabel = 0;
	TRegion* pVertex = (TRegion*) m_pFirstVertex;
	while (pVertex)
	{
		if (nLabels[pVertex->m_nLabel] < 0)  // mapping for this label not determined yet
			SetLabelMapping(nLabels, pVertex, nNewLabel ++);
		pVertex = (TRegion*) pVertex->m_pNext;
	}

	// Obtain the map
	int nWidth = m_pMap->Width();
	int nHeight = m_pMap->Height();
	TGrayImage<int>* pMap = new TGrayImage<int>(nWidth, nHeight);
	int l;
	int r, c;
	for (r = 0; r < nHeight; r ++)
		for (c = 0; c < nWidth; c ++)
		{
			if ((l = (*m_pMap)(r, c)) >= 0)  // class label
				(*pMap)(r, c) = nLabels[l];
			else  // boundary label
				(*pMap)(r, c) = l;
		}

		// Remove interior boundaries
		int nLabel1, nLabel2;
		POINTex Pt;
		for (r = 0; r < nHeight; r ++){
			for (c = 0; c < nWidth; c ++){
				if ((*pMap)(r, c) < 0){  // boundary label
					Pt.Row = r; Pt.Col = c;
					if (!GetNeighboringRegionLabels(nLabel1, nLabel2, Pt))  // interior boundary, the first stores the label
						(*pMap)(r, c) = nLabel1;
				}
			}
		}

		delete[] nLabels;  //QIN_MOD
		return pMap;
}

// 
// Get a new graph with each vertex representing a connected region
// of the same class
// 
TGraph* TRegionMerger::GetConnectedClassRegions()
{
	if (m_enConstructMode == PIXEL)
		return 0;

	TGrayImage<int>* pMap = GetConnectedClassRegionMap();
	TRegionAdjacencyGraph* pRet = new TRegionAdjacencyGraph(pMap, m_pMask, m_ppSrcImg, m_nFeatureCnt, 
		TRegionAdjacencyGraph::GAUSSIAN);
	pRet->ConstructGraph();
	delete pMap;
	return pRet;
}

//
// Set the label mapping for the specified vertex
// 
void TRegionMerger::SetLabelMapping(int* nLabels, TRegion* pCurrent, int nNewLabel){
	TRegion* pNeighbor;
	nLabels[pCurrent->m_nLabel] = nNewLabel;
	for (int i = 0; i < pCurrent->m_nNeighborCnt; i ++){
		pNeighbor = (TRegion*) pCurrent->GetNeighborVertex(pCurrent->m_pEdges[i]);
		if (nLabels[pNeighbor->m_nLabel] < 0 &&  // mapping for neighbor not determined yet
			pNeighbor->m_nClassLabel == pCurrent->m_nClassLabel)  // belong to same class
			SetLabelMapping(nLabels, pNeighbor, nNewLabel);
	}
}
