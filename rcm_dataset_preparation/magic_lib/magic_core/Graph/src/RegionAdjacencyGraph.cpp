// RegionAdjacencyGraph.cpp : implementation of the TRegionAdjacencyGraph class
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "RegionAdjacencyGraph.h"
#include "EstimateUtil.h"

//
// Constructor (Region based)
// 
TRegionAdjacencyGraph::TRegionAdjacencyGraph(TGrayImage<int>* pMap,
                                             TMonoImage* pMask,
                                             TGrayImage<FLOAT>** ppSrcImg,
                                             int nFeatureCnt,
                                             TRegionAdjacencyGraph::statistics_mode enStatMode,
                                             TRegionAdjacencyGraph::boundary_mode enBoundaryMode)
	: TGraph(false)
{
	m_pMap = (TGrayImage<int>*) pMap->Clone();  // provide a copy for processing on (use new operation in Clone)
	m_pMask = pMask;
	m_ppSrcImg = ppSrcImg;
	m_nFeatureCnt = nFeatureCnt;
    m_pVertexLUT = 0;
    m_pEdgeLUT = 0;
    m_enConstructMode = REGION;
	m_enStatisticsMode = enStatMode;
    m_enBoundaryMode = enBoundaryMode;
}


//
// Destructor
// 
TRegionAdjacencyGraph::~TRegionAdjacencyGraph()
{
	delete m_pMap;
	delete[] m_pVertexLUT;
    delete[] m_pEdgeLUT;
	m_BoundaryList.m_pRoot = 0;  // To avoid double deletion. The edges are deleted in the TGraph class
}

//
// Construct the graph
// 
void TRegionAdjacencyGraph::ConstructGraph()
{
    ConstructGraphRegion();

    // Compute the parameters for each region
    int nWidth = m_pMap->Width();
    int nHeight = m_pMap->Height();

    TGaussianRegion* pVertex;

    int i;
    for (int r = 0; r < nHeight; r ++)
        for (int c = 0; c < nWidth; c ++)
        {
            if (m_pMask && !(*m_pMask)(r, c))
                continue;

            i = (*m_pMap)(r, c);
            if (i < 0)  // boundary site
                continue;

            pVertex = (TGaussianRegion*) m_pVertexLUT[i];
            pVertex->AccumulateStatistics(m_ppSrcImg, r, c);

            // Check if it is a border region segment
            if (IsBorderSegment(r, c))
                pVertex->m_bIsBorder = true;
        }

    pVertex = (TGaussianRegion*) m_pFirstVertex;
    while (pVertex)
    {
        pVertex->FinalizeStatistics();
        pVertex = (TGaussianRegion*) pVertex->m_pNext;
    }
    UpdateMaxDegree();
}

//
// generate a struct for the specified region
//
TVertex* TRegionAdjacencyGraph::CreateNewRegion(int nLabel)
{
    if (m_enStatisticsMode == GAUSSIAN)
        return new TGaussianRegion(this, nLabel, m_nFeatureCnt);
    else throw std::runtime_error("Unknown statistics mode!");
}

// 
// Check if it is on the border (no matter grid or mask)
// 
BOOL TRegionAdjacencyGraph::IsBorderSegment(int nRow, int nCol)
{
	if (nRow <= 0 || nRow >= m_pMap->Height() - 1 ||
	    nCol <= 0 || nCol >= m_pMap->Width() - 1)
		return true;

	int r, c;
	for (int i = 0; i < 8; i ++)
	{
		r = nRow + neighbor8[i].Row;
		c = nCol + neighbor8[i].Col;
		if (m_pMask && !(*m_pMask)(r, c))
			return true;
	}

	return false;
}

// 
// Construct the graph in REGION mode
// 
void TRegionAdjacencyGraph::ConstructGraphRegion()
{
	// Allocate memory for the vertex of graph
	TEstimateUtil EstUtil;
	int nVertexCnt = EstUtil.RelabelMap(m_pMap, m_pMask);  // make the labels consecutive
    if(nVertexCnt<1)
		return;

	if (m_pVertexLUT)
		delete[] m_pVertexLUT;

    m_pVertexLUT = new TVertex*[nVertexCnt];
	int i;
	for (i = 0; i < nVertexCnt; i ++)
        m_pVertexLUT[i] = CreateNewRegion(i);

	// Load the boundaries into into m_BoundaryList
    if (m_enBoundaryMode == FULL)
        LoadBoundarySegments();
    else
        LoadBoundarySegments();
        //LoadInnerBoundarySegments();

	// Add vertice in graph
	for (i = 0; i < nVertexCnt; i ++)
		AddVertex(m_pVertexLUT[i]);

    if (m_pEdgeLUT) delete[] m_pEdgeLUT;

    // fill up the edge LUT
    i=0;
    m_pEdgeLUT = new TEdge*[GetEdgeCnt()];
    for (TEdge* pCur = GetEdges(); pCur; pCur = pCur->m_pNext) m_pEdgeLUT[i++]=pCur;
}

// 
// Load the boundaries into a list
// 
void TRegionAdjacencyGraph::LoadBoundarySegments()
{
	int nWidth = m_pMap->Width();
	int nHeight = m_pMap->Height();
	TGrayImage<BYTE>* pVisitRecord = new TGrayImage<BYTE>(nWidth, nHeight);
    pVisitRecord->SetValue(0); // initially all unvisited

    int r, c;
	POINTex Pt;
	for (r = 0; r < nHeight; r ++)
		for (c = 0; c < nWidth; c ++)
	{
		if (m_pMask && !(*m_pMask)(r, c))
			continue;

		if ((*m_pMap)(r, c) < 0 &&  // boundary pixel value is defined as -2 by the watershed
		    (*pVisitRecord)(r, c) == 0)  // not visited before
		{
			Pt.Row = r; Pt.Col = c;
			AddBoundarySegment(pVisitRecord, Pt);
		}
	}
	delete pVisitRecord;
}

//
// Same as the LoadBoundarySegments function except that it does not include the boundary intersection points
// as part of the boundary. This is because I want to handle them as a special case
//
void TRegionAdjacencyGraph::LoadInnerBoundarySegments() {
    int nWidth = m_pMap->Width();
    int nHeight = m_pMap->Height();
    TGrayImage<BYTE>* pVisitRecord = new TGrayImage<BYTE>(nWidth, nHeight);
    pVisitRecord->SetValue(0); // initially all unvisited

    int r, c;
    int rr, cc;
    int nadj;
    POINTex Pt;
    for (r = 0; r < nHeight; r ++)
        for (c = 0; c < nWidth; c ++)
        {
            if (m_pMask && !(*m_pMask)(r, c))
                continue;

            if ((*m_pMap)(r, c) < 0 &&  // boundary pixel value is defined as -2 by the watershed
                (*pVisitRecord)(r, c) == 0)  // not visited before
            {

                nadj = 0;

                // loop over 4-connected neighbours and count the number of non-region pixels
                for (int i = 0; i < 4; i++) {

                    rr = r + neighbor4[i].Row;
                    cc = c + neighbor4[i].Col;

                    if (rr < 0 || rr >= nHeight ||
                        cc < 0 || cc >= nWidth) { // neighbour out of valid range
                        // nadj+=3;
                    }
                    else if (m_pMask && !(*m_pMask)(rr, cc)) {
                        // nadj+=3;
                    } else if ((*m_pMap)(rr, cc) < 0) nadj++;

                }

                // not a boundary intersection if nbranch < 3, so add it to the boundary list
                if (nadj < 3){
                    Pt.Row = r; Pt.Col = c;
                    AddBoundarySegment(pVisitRecord, Pt);
                } else {
                    (*pVisitRecord)(r,c)=1;
                }
            }
        }
    delete pVisitRecord;
}

// 
// Add a boundary segment into the list
// 
// Param in
//  StartPt: the starting point of the boundary segment
// 
// Param out
//  pVisitRecord: the record of visited sites
// 
void TRegionAdjacencyGraph::AddBoundarySegment(TGrayImage<BYTE>* pVisitRecord,
                                               POINTex StartPt)
{
	// Find the two labels of neighboring segments on current boundary site
	int nLabel1, nLabel2;
	if (!GetNeighboringRegionLabels(nLabel1, nLabel2, StartPt))
		return;

	// Compute the chain corresponding to the specified site
	TChainCode* pChain = ComputeChainCode(StartPt, nLabel1, nLabel2,
	                                      pVisitRecord);

	// Add boundary segment into the list
	TBoundary* pCurrent = IsBoundarySegmentInList(nLabel1, nLabel2);
	if (pCurrent)  // Boundary segment with same labels found, concatenate
	{
		pCurrent->AttachChain(pChain);
	}
	else  // Boundary segment with same labels not found, add
	{
		// Insert a new boundary segment into the linked list at the beginning position
		pCurrent = CreateNewBoundary(nLabel1, nLabel2, pChain, true);
		m_BoundaryList.Insert(pCurrent);
	}
}


// 
// Compute the corresponding chain code
// 
TChainCode* TRegionAdjacencyGraph::ComputeChainCode(POINTex StartPt,
                                                    int nLabel1, int nLabel2,
                                                    TGrayImage<BYTE>* pVisitRecord)
{
	int nWidth = m_pMap->Width();
	int nHeight = m_pMap->Height();

	// Mark current site as visited
	POINTex CurPt = StartPt, NextPt;
	(*pVisitRecord)(CurPt.Row, CurPt.Col) = 1;
    BOOL isb;

	TChainCode* pChain = new TChainCode(StartPt, 100);
	for (int i = 0; i < 4; i ++)  // 8-neighbor system is used for region, and thus 4-neighbor system for boundary
	{
		NextPt.Row = CurPt.Row + neighbor4[i].Row;
		NextPt.Col = CurPt.Col + neighbor4[i].Col;
		if (NextPt.Row < 0 || NextPt.Row >= nHeight ||
		    NextPt.Col < 0 || NextPt.Col >= nWidth)
			continue;
		if (m_pMask && !(*m_pMask)(NextPt.Row, NextPt.Col))
			continue;
		if ((*m_pMap)(NextPt.Row, NextPt.Col) < 0 &&  // boundary site
		    (*pVisitRecord)(NextPt.Row, NextPt.Col) == 0) // not visited before
		{

            if (1) // m_enBoundaryMode == FULL
                isb = IsBoundaryOfRegions(nLabel1, nLabel2, NextPt);
            else
                isb =  IsInnerBoundaryOfRegions(nLabel1, nLabel2, NextPt);

            if (isb){
                pChain->AddSite(NextPt, 0);  // here the second argument is useless
                CurPt = NextPt;
                i = 0;
                // Mark current site as visited
                (*pVisitRecord)(CurPt.Row, CurPt.Col) = 1;
            }

		}
	}
	return pChain;
}

// 
// Get the labels of two neighboring segments on the specifed
// boundary site. If more than two different segments intersect,
// just choose the first two found.
// 
// Param in
//  Pt: position of the site
// 
// Param out
//  nLabel1, nLabel2: the labels of the two segments
// 
// Return
//  true: the specified site is a boundary site (with at least two
//   different neighbors; 
//  false: not a boundary site
// 
BOOL TRegionAdjacencyGraph::GetNeighboringRegionLabels(int& nLabel1, int& nLabel2,
                                                       POINTex Pt)
{
	int nWidth = m_pMap->Width();
	int nHeight = m_pMap->Height();
	int r, c, k;

    if (1){
        nLabel1 = -1; nLabel2 = -1;
        for (int i = 0; i < 8; i ++)
        {
            r = Pt.Row + neighbor8[i].Row;
            c = Pt.Col + neighbor8[i].Col;
            if (r < 0 || r >= nHeight ||
                c < 0 || c >= nWidth)  // out of valid range
                continue;
            if (m_pMask && !(*m_pMask)(r, c))
                continue;
            if ((k = (*m_pMap)(r, c)) < 0)  // boundary site
                continue;
            if (nLabel1 == -1)  // first label
                nLabel1 = k;
            else if (nLabel1 != k)
            {
                nLabel2 = k;
                break;
            }
        }
        return (nLabel2 >= 0);
    } else {
        nLabel1 = -1; nLabel2 = -1;
        for (int i = 0; i < 4; i ++)
        {
            r = Pt.Row + neighbor4[i].Row;
            c = Pt.Col + neighbor4[i].Col;
            if (r < 0 || r >= nHeight ||
                c < 0 || c >= nWidth)  // out of valid range
                continue;
            if (m_pMask && !(*m_pMask)(r, c))
                continue;
            if ((k = (*m_pMap)(r, c)) < 0)  // boundary site
                continue;
            if (nLabel1 == -1)  // first label
                nLabel1 = k;
            else if (nLabel1 != k)
            {
                nLabel2 = k;
                break;
            }
        }
        return (nLabel2 >= 0);
    }

}

// 
// Check if boundary segment with the same labels already 
// exists in the list
// 
// Param in
//  nLabel1, nLabel2: the labels of the two segments
// 
// Return
//  the boundary segment in the list found
// 
TBoundary* TRegionAdjacencyGraph::IsBoundarySegmentInList(int nLabel1, int nLabel2)
{
	TBoundaryKey Key(nLabel1, nLabel2);
	TAVLTreeNode* p = m_BoundaryList.Find(Key);
	if (p)
		return (TBoundary*) p;
	else
		return 0;
}

// 
// Check if the site is boundary of segments specified by the two labels
// 
// Param in
//  nLabel1, nLabel2: the labels of the two specified segments
//  Pt: position of the site
// 
BOOL TRegionAdjacencyGraph::IsBoundaryOfRegions(int nLabel1, int nLabel2, 
                                                POINTex Pt)
{
	int nWidth = m_pMap->Width();
	int nHeight = m_pMap->Height();
	int r, c, k;
	BOOL bFound1 = false, bFound2 = false;
	for (int i = 0; i < 8; i ++)
	{
		r = Pt.Row + neighbor8[i].Row;
		c = Pt.Col + neighbor8[i].Col;
		if (r < 0 || r >= nHeight ||
		    c < 0 || c >= nWidth)  // out of valid range
			continue;
		if (m_pMask && !(*m_pMask)(r, c))
			continue;
		if ((k = (*m_pMap)(r, c)) < 0)  // boundary site
			continue;
		if (nLabel1 == k)  // first label found
			bFound1 = true;
		else if (nLabel2 == k)
			bFound2 = true;
	}

	return (bFound1 && bFound2);
}

//
// Check if the site is boundary of segments specified by the two labels, *excluding boundary intersection points*
//
// Param in
//  nLabel1, nLabel2: the labels of the two specified segments
//  Pt: position of the site
//
BOOL TRegionAdjacencyGraph::IsInnerBoundaryOfRegions(int nLabel1, int nLabel2,
                                                POINTex Pt)
{
//    std::cerr << "yyyyyyyyyyyyyyyyyy" << std::endl;
//    int nWidth = m_pMap->Width();
//    int nHeight = m_pMap->Height();
//    int r, c, k;
//    int nbranch=0;
//    BOOL bFound1 = false, bFound2 = false;
//
//    for (int i = 0; i < 8; i ++)
//    {
//        r = Pt.Row + neighbor8[i].Row;
//        c = Pt.Col + neighbor8[i].Col;
//        if (r < 0 || r >= nHeight ||
//            c < 0 || c >= nWidth)  // out of valid range
//            continue;
//        if (m_pMask && !(*m_pMask)(r, c)) // pixel not masked
//            continue;
//        if ((k = (*m_pMap)(r, c)) < 0) // boundary site
//        {
//            if (i % 2 == 0) nbranch++; //count how many boundary sites touch this one
//            continue;
//        }
//        if (nLabel1 == k)  // first label found
//            bFound1 = true;
//        else if (nLabel2 == k)
//            bFound2 = true;
//    }
//
//    return (bFound1 && bFound2 && (nbranch < 3));
    int nWidth = m_pMap->Width();
    int nHeight = m_pMap->Height();
    int r, c, k;
    BOOL bFound1 = false, bFound2 = false;
    for (int i = 0; i < 4; i ++)
    {
        r = Pt.Row + neighbor4[i].Row;
        c = Pt.Col + neighbor4[i].Col;
        if (r < 0 || r >= nHeight ||
            c < 0 || c >= nWidth)  // out of valid range
            continue;
        if (m_pMask && !(*m_pMask)(r, c))
            continue;
        if ((k = (*m_pMap)(r, c)) < 0)  // boundary site
            continue;
        if (nLabel1 == k)  // first label found
            bFound1 = true;
        else if (nLabel2 == k)
            bFound2 = true;
    }

    return (bFound1 && bFound2);
}


// 
// Get the first neighbor for constructing the graph
// 
TVertex* TRegionAdjacencyGraph::FirstNeighbor(TVertex* pVertex)
{
	int nLabel = pVertex->m_nLabel;
	TBoundaryKey lKey(nLabel, 0);
	TBoundaryKey hKey(nLabel, nLabel - 1);  // since segments with label greater than
	                                        // current is not in graph yet, those
	                                        // segments are excluded.
	TBoundary* pCur = (TBoundary*) m_BoundaryList.FindFirstInRange(lKey, hKey);
	if (pCur)
	{
		if (pCur->m_nLabel1 == nLabel)
			return m_pVertexLUT[pCur->m_nLabel2];
		else
			return m_pVertexLUT[pCur->m_nLabel1];
	}

	return 0;
}

// 
// Get the next neighbor for constructing the graph
// 
TVertex* TRegionAdjacencyGraph::NextNeighbor(TVertex* pVertex)
{
	TBoundary* pCur = (TBoundary*) m_BoundaryList.FindNextInRange();
	if (pCur)
	{
		if (pCur->m_nLabel1 == pVertex->m_nLabel)
			return m_pVertexLUT[pCur->m_nLabel2];
		else
			return m_pVertexLUT[pCur->m_nLabel1];
	}

	return 0;
}

// 
// Construct edge for the two specified region segments
// 
TEdge* TRegionAdjacencyGraph::ConstructEdge(TVertex* pNeighbor1, TVertex* pNeighbor2)
{
	// In the graph construction stage, use the existing elements
	// in the list directly
	TBoundaryKey Key(pNeighbor1->m_nLabel, pNeighbor2->m_nLabel);
	TAVLTreeNode* p = m_BoundaryList.Find(Key);
	if (p)
	{
		((TBoundary*) p)->SetVertexNeighbor(pNeighbor1, pNeighbor2);
		return (TBoundary*) p;
	}
	else
		return 0;
}

//
// Set the label mapping for the specified vertex
//
void TRegionAdjacencyGraph::SetLabelMapping(int* nLabels, TRegion* pCurrent, int nNewLabel){
    TRegion* pNeighbor;
    nLabels[pCurrent->m_nLabel] = nNewLabel;
    for (int i = 0; i < pCurrent->m_nNeighborCnt; i ++){
        pNeighbor = (TRegion*) pCurrent->GetNeighborVertex(pCurrent->m_pEdges[i]);
        if (nLabels[pNeighbor->m_nLabel] < 0 &&  // mapping for neighbor not determined yet
            pNeighbor->m_nClassLabel == pCurrent->m_nClassLabel)  // belong to same class
            SetLabelMapping(nLabels, pNeighbor, nNewLabel);
    }
}


BOOL TRegionAdjacencyGraph::SetClassLabelByIndex(int ClassLabel, int Index)
{
    TRegion* pRegion = (TRegion*) m_pVertexLUT[Index];

    if (!pRegion) return false; // region does not exist

    pRegion->m_nClassLabel = ClassLabel;
    return true;
}



