// Region.cpp : implementation of the TGaussianRegion class
//
// Class for region segment struct with Gaussian statistics.
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "GaussianRegion.h"
#include "Boundary.h"
#include "GrayImage.h"
#include "MonoImage.h"
#include "RegionAdjacencyGraph.h"

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

// 
// Constructor
//
TGaussianRegion::TGaussianRegion(TGraph* pGraph, int nLabel, int nFeatureCnt)
	: TRegion(pGraph, nLabel)
{
	m_pMean = new TVector(nFeatureCnt, false);
	m_pMean->SetZero();
	m_pCoMean = new TMatrix(nFeatureCnt, nFeatureCnt);
	m_pCoMean->SetZero();

}

// 
// Destructor
// 
TGaussianRegion::~TGaussianRegion()
{
	if (m_pMean) delete m_pMean;
	if (m_pCoMean) delete m_pCoMean;

}

//
// Accumulate pixel statistics during graph construction
//
void TGaussianRegion::AccumulateStatistics(TGrayImage<FLOAT>** ppSrcImg, int r, int c)
{

    double dbVal;
    int nFeatureCnt = ((TRegionAdjacencyGraph*) m_pGraph)->GetFeatureCnt();

    // update mean and comean of the region
    for ( int i = 0; i < nFeatureCnt; i++ )
    {
        dbVal = (*ppSrcImg[i])(r, c);
        (*m_pMean)(i) += dbVal;
        m_pCoMean->m_ppData[i][i] += dbVal * dbVal;
        for ( int j = i + 1; j < nFeatureCnt; j++ )
            m_pCoMean->m_ppData[i][j] += dbVal * (*ppSrcImg[j])(r, c);
    }

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

//
// Finalize region statistics once graph construction is complete
//
void TGaussianRegion::FinalizeStatistics()
{
    int nFeatureCnt = ((TRegionAdjacencyGraph*) m_pGraph)->GetFeatureCnt();

    for (int j = 0; j < nFeatureCnt; j ++)
    {
        (*m_pMean)(j) /= m_nPixelCnt;
        for (int k = j; k < nFeatureCnt; k ++)
            m_pCoMean->m_ppData[j][k] /= m_nPixelCnt;
    }
    m_dbXMean /= m_nPixelCnt;
    m_dbYMean /= m_nPixelCnt;
}

//
// Compute the unary weight of the vertex
// 
double TGaussianRegion::GetVertexWeight()
{
	// Compute the covariance
	int nFeatureCnt = m_pMean->GetElementCnt();
	TMatrix* pCov = new TMatrix(nFeatureCnt, nFeatureCnt);
	double** ppCovData = pCov->m_ppData;
	double* pMean = m_pMean->m_ppData[0];
	double** pCoMean = m_pCoMean->m_ppData;
	for (int i = 0; i < nFeatureCnt; i ++)
	{
		ppCovData[i][i] = pCoMean[i][i] - pMean[i]*pMean[i];
		for (int j = i + 1; j < nFeatureCnt; j ++)
		{
			ppCovData[i][j] = pCoMean[i][j] - pMean[i]*pMean[j];
			ppCovData[j][i] = ppCovData[i][j];
		}
	}

	// Determinant
	double dbSigma;
	if (m_nPixelCnt <= 3)
		dbSigma = 1; // can't get a good estimate so just set it to 1
	else
		dbSigma = sqrt(pCov->Determinant());  //Note: it is the square root of the determinant value

	if (dbSigma <= sqrt(ZERO_THRESH) )  // if very small (near zero), let the std to be half of the intensity precision
		dbSigma = sqrt(ZERO_THRESH);

	delete pCov;

	return ( log(dbSigma) *  m_nPixelCnt / nFeatureCnt );
}

// 
// Merge with a specified region
// 
void TGaussianRegion::MergeWith(TVertex* pNeighbor)
{
	TMonoImage* pMask = ((TRegionAdjacencyGraph*) m_pGraph)->GetMask();
	TGrayImage<int>* pMap = ((TRegionAdjacencyGraph*) m_pGraph)->GetMap();
	TGrayImage<FLOAT>** ppSrcImg= ((TRegionAdjacencyGraph*) m_pGraph)->GetImageSource();

	int nFeatureCnt = m_pMean->GetElementCnt();
	TGaussianRegion* p = (TGaussianRegion*) pNeighbor;
	double* pMean = m_pMean->m_ppData[0];
	double** pCoMean = m_pCoMean->m_ppData;
	int i, j;
	for (i = 0; i < nFeatureCnt; i ++)
	{
		pMean[i] = pMean[i] * m_nPixelCnt + (*p->m_pMean)(i) * p->m_nPixelCnt;
		for (j = i; j < nFeatureCnt; j ++)
			pCoMean[i][j] = pCoMean[i][j] * m_nPixelCnt + 
			                (*p->m_pCoMean)(i, j) * p->m_nPixelCnt;
	}
	m_dbXMean = m_dbXMean * m_nPixelCnt + p->m_dbXMean * p->m_nPixelCnt;
	m_dbYMean = m_dbYMean * m_nPixelCnt + p->m_dbYMean * p->m_nPixelCnt;
	m_nPixelCnt = m_nPixelCnt + p->m_nPixelCnt;

	// Update the bounding rectangle
	if (m_nLeft > p->m_nLeft) m_nLeft = p->m_nLeft;
	if (m_nRight < p->m_nRight) m_nRight = p->m_nRight;
	if (m_nTop > p->m_nTop) m_nTop = p->m_nTop;
	if (m_nBottom < p->m_nBottom) m_nBottom = p->m_nBottom;

	// Flag if border segment
	m_bIsBorder = m_bIsBorder | p->m_bIsBorder;

	// Update the labels in the map
	for (i = p->m_nTop; i <= p->m_nBottom; i ++)
		for (j = p->m_nLeft; j <= p->m_nRight; j ++)
	{
		if (pMask && !(*pMask)(i, j))
			continue;

		if ((*pMap)(i, j) == p->m_nLabel)
			// Change labels in the map corresponding to the segment being merged
			(*pMap)(i, j) = m_nLabel;
	}

	// Update the labels of previous boundary sites
	double dbVal;
	TBoundary* pBoundary = (TBoundary*) GetEdge(pNeighbor);
	if (pBoundary && pBoundary->m_bOnGrid)  // the boundary sites are on grid, so region features need to account for those sites
	{
		TChainCode* pChain = pBoundary->m_pChain;
		POINTex Pt;
		int nLabel;
		while (pChain)
		{
			pChain->GetFirstSite(Pt);
			do {
				if (GetNeighboringRegionLabels(nLabel, Pt))
				{
					if (!PreserveNonDeletableSite(this, nLabel, Pt))  // associated boundary not found, look in neighbor
						PreserveNonDeletableSite((TRegion*) pNeighbor, nLabel, Pt);
				}
				else
				{
					// Change boundary labels between the two segments
					(*pMap)(Pt.Row, Pt.Col) = m_nLabel;

					// Update the bounding rectangle
					if (m_nLeft > Pt.Col) m_nLeft = Pt.Col;
					if (m_nRight < Pt.Col) m_nRight = Pt.Col;
					if (m_nTop > Pt.Row) m_nTop = Pt.Row;
					if (m_nBottom < Pt.Row) m_nBottom = Pt.Row;

					// Update the border indicating flag
				    if (Pt.Row <= 0 || Pt.Row >= pMap->Height() - 1 || Pt.Col <= 0 || Pt.Col >= pMap->Width() - 1)
		                m_bIsBorder = 1;
					else
					{
						int r, c;
						for (i = 0; i < 8; i ++)
						{
							r = Pt.Row + neighbor8[i].Row;
							c = Pt.Col + neighbor8[i].Col;
							if (pMask && !(*pMask)(r, c))
								m_bIsBorder = 1;
						}
					}
					// Include those boundary sites in computing the mean and deviation
					for (i = 0; i < nFeatureCnt; i ++)
					{
						dbVal = (*ppSrcImg[i])(Pt.Row, Pt.Col);
						pMean[i] += dbVal;
						pCoMean[i][i] += dbVal * dbVal;
						for (int j = i + 1; j < nFeatureCnt; j ++)
							pCoMean[i][j] += dbVal * (*ppSrcImg[j])(Pt.Row, Pt.Col);
					}

					m_dbXMean += Pt.Col;
					m_dbYMean += Pt.Row;

					m_nPixelCnt ++;
				}
			} while (pChain->GetNextSite(Pt));

			pChain = pChain->m_pNext;
		}
	}

	for (i = 0; i < nFeatureCnt; i ++)
	{
		pMean[i] = pMean[i] / m_nPixelCnt;
		for (j = i; j < nFeatureCnt; j ++)
            // Note: only the upper triangle of m_pCoMean is assigned values
			pCoMean[i][j] = pCoMean[i][j] / m_nPixelCnt;
	}
	m_dbXMean = m_dbXMean / m_nPixelCnt;
	m_dbYMean = m_dbYMean / m_nPixelCnt;
}

