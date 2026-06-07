// RegionClustering.cpp : implementation of the TRegionClustering class
//
// Clustering the regions by using the MRF on GIEP graph
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "RegionClustering.h"
#include <time.h>
//#include "StrengthPenaltyBoundary.h"
//#include "EstimateUtil.h"
//#include "Matrix.h"  // For n-D minimum Fisher criterion
//#include "RegionBasedClustering.h"
//#include "PixelBasedClustering.h"
//#include "GIEPGraph.h"

#define REG_COEF 0.001

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

/// Constructor which sets up some parameters.
/**
 * Param in
 * nType: Type of clustering energy equation to use.
 * nInitMode: Initialization mode to use.
 * nInitInitMode: Initialization mode to use for to initialize the
 * initialization algorithm (which is a clustering algorithm in
 * itself).
 * dbRatio: Parameter for the initialization algorithm, controls the
 * stopping criterion.
 *
 * nInit1: Parameter for the initialization algorithm, the number of
 * iterations of the initialization clustering algorithm.
 * nInit2: Parameter for the initialization algorithm, how many times
 * the clustering algorithm is executed
 *
 * The initialization clustering algorithm will be run nInit2 times,
 * each time starting from a new random starting point. Each time, the
 * clustering algorithm will be run for nInit1 iterations. Out of the
 * nInit2 times, the result with the lowest energy is chosen.
 *
 * nIterFinal : After the lowest energy initilization is chosen, the
 * clustering algorithm is further run for this many interations. This
 * last result is used as the starting point for IRGS.
 *
 * nMode: mode for annealing
 * dbInitT: initial temperature
 * nIterations: number of iterations. If SINGLE_STEP has been enabled
 * in the mode, this parameter has no effect.
 * nClusterCnt: number of classes (clusters)
 */	 
TRegionClustering::TRegionClustering(int nType, int nInitMode, int nInitInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal,
									 int nMode, double dbInitT, int nIterations, int nClassCnt)
	: TMarkovNetClassifier(nMode, dbInitT, nIterations, nClassCnt)
{
	m_nInitMode = nInitMode;
	m_nInitInitMode = nInitInitMode;
	m_dbInitRatio = dbRatio;
	m_nInit1 = nInit1;
	m_nInit2 = nInit2;
	m_nInitIterFinal = nIterFinal;

	m_nType = nType;

	//srand((unsigned) time(NULL));
}

//
// Destructor
// 
TRegionClustering::~TRegionClustering()
{

}


// 
// Final classification, in which boundary penalty between neighboring segments is not used
// 
bool TRegionClustering::FinalClassify()
{
	double U, dbMin;
	int i, nMinIdx;
	int* nLabelFlag = new int[m_nClassCnt];
	memset(nLabelFlag, 0, m_nClassCnt * sizeof(int));
	TRegion* p = (TRegion*) m_pGraph->GetVertices();
	while (p)
	{
		dbMin = 999999999999;  // big enough
		for ( i = 0; i < m_nClassCnt; i++ )
		{
			if ( ( U = GetUnaryEnergy(p, i) ) < dbMin )
			{
				dbMin = U;
				nMinIdx = i;
			}
		}
		p->m_nClassLabel = nMinIdx;

		if ( nLabelFlag[nMinIdx] == 0 )
			nLabelFlag[nMinIdx] = 1;

		p = (TRegion*) p->m_pNext;
	}

	nMinIdx = 0;
	for ( i = 0; i < m_nClassCnt; i++ )
		nMinIdx += nLabelFlag[i];

	delete[] nLabelFlag;

	if ( nMinIdx == m_nClassCnt )
		return true;
	else
		return false;
}


// 
// Generate a random class label configuration
//
void TRegionClustering::GenerateRandomLabel()
{
	TRegion* pCur = (TRegion*) m_pGraph->GetVertices();
    int i, nLabel;
	int* nLabelFlag = new int[m_nClassCnt];
    bool bContinueFlag;
	do
	{
		memset(nLabelFlag, 0, m_nClassCnt * sizeof(int));
		while (pCur)
		{
			nLabel = rand() * m_nClassCnt / (RAND_MAX + 1.0);
			pCur->m_nClassLabel = nLabel;
		    if ( nLabelFlag[nLabel] == 0 )
				nLabelFlag[nLabel] = 1;
			pCur = (TRegion*) pCur->m_pNext;
		}
		
		nLabel = 0;
		for ( i = 0; i < m_nClassCnt; i++ )
			nLabel += nLabelFlag[i];

		if ( nLabel == m_nClassCnt )
			bContinueFlag = false;
		else
			bContinueFlag = true;
	} while (bContinueFlag);

	delete[] nLabelFlag;
}

//
// Check if the current site is at boundary
// Old method, kept in case needed.
// 
/*bool TRegionClustering::IsBoundary(TGrayImage<int>* pMap, int nRow, int nCol)
{
	TMonoImage* pMask = m_pGraph->GetMask();

	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int i, r, c;
	int l = (*pMap)(nRow, nCol);

	for ( i = 0; i < 8; i++ )
	{
		r = nRow + neighbor8[i].Row;
		c = nCol + neighbor8[i].Col;
		if ( r < 0 || r >= nHeight || c < 0 || c >= nWidth )
			continue;
		if ( pMask && !(*pMask)(r, c) )
			continue;

		if ( (*pMap)(r, c) != l )
		{
			return true;
		}

		
	}
	
	return false;
}*/

/// Checks a pixel location is a boundary point. 
/** Also, populates nClassLabel array
 * to indicate which classes are adjacent to pixel. nClassLabel must be
 * already allocated. bIndepBoundary is true when the pixel is a
 * boundary but has no neighbours with a valid class.
 */ 
bool TRegionClustering::IsBoundary(TGrayImage<int>* pMap, TMonoImage* pMask, int nRow, int nCol, int*& nClassLabel, bool& bIndepBoundary)
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int i, r, c;
	int l = (*pMap)(nRow, nCol);
	int l2;
	bool bRet = false;

	if (l < 0)
	{
		bIndepBoundary = true;

		for ( i = 0; i < 8; i++ )
		{
			r = nRow + neighbor8[i].Row;
			c = nCol + neighbor8[i].Col;

			if ( r < 0 || r >= nHeight || c < 0 || c >= nWidth )
				continue;
			if ( pMask && !(*pMask)(r, c) )
				continue;

			l2 = (*pMap)(r, c);

			if (l2 < 0)
				continue;

			nClassLabel[l2] = 1; // this class exists in the neighbourhood.
			bIndepBoundary = false; // at least one valid class beside this pixel.

		}

		return true;

	}

	bIndepBoundary = false;
	return false;
}

/// Gets map where each region has its proper label applied but boundary pixels are still not labeled.
TGrayImage<int>* TRegionClustering::GetLabeledGraphMap()
{	
	TGrayImage<int>* tpMap = m_pGraph->GetMap();
	TMonoImage* pMask = m_pGraph->GetMask();

	// Get maximum number of remaining segments
	TEstimateUtil* pEstUtil = new TEstimateUtil();
	int nMaxCnt = pEstUtil->GetMaxLabel(tpMap, pMask) + 1;
	delete pEstUtil;

	// Generate the label mapping
	int* nLabels = new int[nMaxCnt];
	memset(nLabels, 0, nMaxCnt * sizeof(int));
	TRegion* pVertex = (TRegion*) m_pGraph->GetVertices();
	while (pVertex)
	{
		nLabels[pVertex->m_nLabel] = pVertex->m_nClassLabel;
		pVertex = (TRegion*) pVertex->m_pNext;
	}

	// Obtain the classification map
	int nWidth = tpMap->Width();
	int nHeight = tpMap->Height();
	TGrayImage<int>* pMap = new TGrayImage<int>(tpMap);
	int l, r, c;
	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;
			if ( ( l = (*tpMap)(r, c) ) >= 0 )  // class label
				(*pMap)(r, c) = nLabels[l];
		}

	delete[] nLabels;

	int nRow, nCol, l2, i;
	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;
			if ( (*pMap)(r, c) < 0 )
			{
				// Checks pixels around current boundary pixel.
				// If all non-boundary pixels that touch the current one
				// belong to the same class, then the current one also
				// belongs to that class and is no longer a boundary point.
				l2 = -1;
				for ( i = 0; i < 8; i++ )
				{
					nRow = r + neighbor8[i].Row;
					nCol = c + neighbor8[i].Col;

					if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
						continue;
					if ( pMask && !(*pMask)(nRow, nCol) )  // not valid
						continue;
					if ( ( l = (*pMap)(nRow, nCol)) < 0 )  // also boundary
						continue;
					if ( l != l2 && l2 >= 0 )
						break;

					l2 = l;
				}

				if ( i >= 8 && l2 >= 0 )  // inside the same class
					(*pMap)(r, c) = l2;
			}
		}
	
	return pMap;
}
