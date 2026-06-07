// Region.cpp : implementation of the TRegion class
//
// Class for region segment struct
/////////////////////////////////////////////////////////////////////////////

#include "Region.h"
#include "RegionAdjacencyGraph.h"

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

// 
// Constructor
//
TRegion::TRegion(TGraph* pGraph, int nLabel)
	: TVertex(pGraph, nLabel)
{
	m_nClassLabel = -1;

	m_dbXMean = 0;
	m_dbYMean = 0;

	m_dbEdgeMean = 0;
	m_nPixelCnt = 0;

	m_nTop = 100000; m_nLeft = 100000;  // large enough
	m_nBottom = -1; m_nRight = -1;
	m_bIsBorder = false;
}

// 
// Destructor
// 
TRegion::~TRegion()
{

}

// 
// Get the boundary length
// 
int TRegion::GetBoundaryLength()
{
	int nBoundaryLength = 0;
	for (int i = 0; i < m_nNeighborCnt; i ++)
		nBoundaryLength += ((TBoundary*) m_pEdges[i])->m_nLength;

	return nBoundaryLength;
}

/// Accumulate pixel statistics during graph construction
void TRegion::AccumulateStatistics(TGrayImage<FLOAT>** ppSrcImg, int r, int c)
{
    // Position parameters
    m_dbXMean += c;
    m_dbYMean += r;
    m_nPixelCnt ++;

    // Bounding rectangle
    if (r < m_nTop) m_nTop = r;
    if (r > m_nBottom) m_nBottom = r;
    if (c < m_nLeft) m_nLeft = c;
    if (c > m_nRight) m_nRight = c;
}

/// Finalize region statistics once graph construction is complete
void TRegion::FinalizeStatistics()
{
    m_dbXMean /= m_nPixelCnt;
    m_dbYMean /= m_nPixelCnt;
}


// 
// Preserve the site to a boundary chain associated with the specified
// region segment since it is still boundary site after merging
// 
// Return
//  the associated boundary found
//  NULL: the associated boundary not found for current region segment
// 
TBoundary* TRegion::PreserveNonDeletableSite(TRegion* pCurrent,
                                             int nLabel, POINTex Pt)
{
	TBoundary* pBoundary;
	int i;
	for (i = 0; i < pCurrent->m_nNeighborCnt; i ++)
	{
		pBoundary = (TBoundary*) pCurrent->m_pEdges[i];
		if (pBoundary->m_nLabel1 == nLabel || pBoundary->m_nLabel2 == nLabel)
			break;
	}
	if (i >= pCurrent->m_nNeighborCnt)  // not found
		return 0;

	pBoundary->AddSite(Pt, 0);  // the second argument is useless here
	return pBoundary;
}

// 
// Get a neighboring segment label on the specifed site. If more than
// two different segments intersect, just choose the first found.
// 
// Param in
//  Pt: position of the site
// 
// Param out
//  nLabel: the obtained label
// 
// Return
//  true: the specified site is a boundary site (with at least one
//   different neighbor; 
//  false: not boundary site
// 
BOOL TRegion::GetNeighboringRegionLabels(int& nLabel, POINTex Pt)
{
	TMonoImage* pMask = ((TRegionAdjacencyGraph*) m_pGraph)->GetMask();
	TGrayImage<int>* pMap = ((TRegionAdjacencyGraph*) m_pGraph)->GetMap();
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int r, c;
	for (int i = 0; i < 8; i ++)
	{
		r = Pt.Row + neighbor8[i].Row;
		c = Pt.Col + neighbor8[i].Col;
		if (r < 0 || r >= nHeight ||
		    c < 0 || c >= nWidth)  // out of valid range
			continue;
		if (pMask && !(*pMask)(r, c))
			continue;
		if ((nLabel = (*pMap)(r, c)) < 0)  // boundary site
			continue;
		if (nLabel != m_nLabel)
			return true;
	}

	return false;
}

