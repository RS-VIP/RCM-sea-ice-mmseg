// GaussianStrengthPenaltyBoundary.cpp : implementation of the TGaussianStrengthPenaltyBoundary class
//
// Class for region boundary with strength penalty using Gaussian statistics.
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "GaussianStrengthPenaltyBoundary.h"
#include "GaussianRegion.h"
#include "GIEPGraph.h"

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

// 
// Constructor
//
TGaussianStrengthPenaltyBoundary::TGaussianStrengthPenaltyBoundary(TGraph* pGraph,
                                                   int nLabel1, int nLabel2,
                                                   TChainCode* pChain,
                                                   BOOL bOnGrid)
	: TStrengthPenaltyBoundary(pGraph, nLabel1, nLabel2, pChain, bOnGrid)
{

}

// 
// Destructor
// 
TGaussianStrengthPenaltyBoundary::~TGaussianStrengthPenaltyBoundary()
{
}

// 
// Compute the edge penalty energy at a given site
//
double TGaussianStrengthPenaltyBoundary::GetPenaltyEnergyAtSite(POINTex Pt, int nDir)
{
	/*TGrayImage<FLOAT>** pEdgeStrength = ((TGIEPGraph*) m_pGraph)->GetEdgePenalty();*/
	double K = ((TGIEPGraph*) m_pGraph)->GetCurK(), dbVal;
	int nFeatureCnt = ((TGIEPGraph*) m_pGraph)->GetFeatureCnt();
//	if (((TGIEPGraph*) m_pGraph)->GetMode() == TRegionAdjacencyGraph::PIXEL)
//	{
//	    TGrayImage<FLOAT>** ppSrcImg = ((TGIEPGraph*) m_pGraph)->GetImageSource();
//
//		// Since the coordinate of boundary sites are defined on a grid
//		// 2 times the original grid for the image, need to divide by 2
//		int r = (Pt.Row - neighbor8[nDir].Row) >> 1;
//		int c = (Pt.Col - neighbor8[nDir].Col) >> 1;
//		// Here the edge strength is the difference between neighboring
//		// sites
//		dbVal = 0;
//		for (int k = 0; k < nFeatureCnt; k++)
//			dbVal += ((*ppSrcImg[k])(r, c) - (*ppSrcImg[k])(r + neighbor8[nDir].Row, c + neighbor8[nDir].Col))*((*ppSrcImg[k])(r, c) - (*ppSrcImg[k])(r + neighbor8[nDir].Row, c + neighbor8[nDir].Col));
//		dbVal = sqrt(dbVal)/K;
//	}
//	else
//	{
//		TGrayImage<FLOAT>* pGradient =((TGIEPGraph*) m_pGraph)->GetGradient();
//		// Here the edge strength is the gradient magnitude at the site
//		dbVal = (*pGradient)(Pt.Row, Pt.Col) / K;
//	}

    TGrayImage<FLOAT>* pGradient =((TGIEPGraph*) m_pGraph)->GetGradient();
    // Here the edge strength is the gradient magnitude at the site
    dbVal = (*pGradient)(Pt.Row, Pt.Col) / K;

	return exp(- dbVal * dbVal);
}

// 
// Compute region energy difference between before merging and after merging
// 
void TGaussianStrengthPenaltyBoundary::ComputeRegionDiffEnergy()
{
	// Compute the covariance
	TGaussianRegion* p1 = (TGaussianRegion*) m_pNeighbor1;
	TGaussianRegion* p2 = (TGaussianRegion*) m_pNeighbor2;
	int nMergedPixelCnt = p1->m_nPixelCnt + p2->m_nPixelCnt;
	int nFeatureCnt = p1->m_pMean->GetElementCnt();
	TVector dbMergedMean(nFeatureCnt, false);
	TMatrix dbMergedCoMean(nFeatureCnt, nFeatureCnt);
	int i;
	for (i = 0; i < nFeatureCnt; i ++)
	{
		dbMergedMean(i) = ((*p1->m_pMean)(i) * p1->m_nPixelCnt +
		                   (*p2->m_pMean)(i) * p2->m_nPixelCnt) /
		                  nMergedPixelCnt;
		for (int j = i; j < nFeatureCnt; j ++)
			dbMergedCoMean(i, j) = ((*p1->m_pCoMean)(i, j) * p1->m_nPixelCnt +
			                        (*p2->m_pCoMean)(i, j) * p2->m_nPixelCnt) /
			                       nMergedPixelCnt;
	}
	for (i = 0; i < nFeatureCnt; i ++)
	{
		dbMergedCoMean(i, i) -= dbMergedMean(i) * dbMergedMean(i);
		for (int j = i + 1; j < nFeatureCnt; j ++)
		{
			dbMergedCoMean(i, j) -= dbMergedMean(i) * dbMergedMean(j);
			dbMergedCoMean(j, i) = dbMergedCoMean(i, j);
		}
	}

	double dbMergedSigma;
	if (nMergedPixelCnt <= 3)
		dbMergedSigma = 1;
	else
		// Compute the covariance determinant
		dbMergedSigma = sqrt(dbMergedCoMean.Determinant());
	if ( dbMergedSigma <= sqrt(ZERO_THRESH) )  // if very small (near zero), let the std to be half of the intensity precision
		dbMergedSigma = sqrt(ZERO_THRESH); // QIN_WHY 0.000001
	m_dbRegionEnergyDiff = ( log(dbMergedSigma) * nMergedPixelCnt / nFeatureCnt -  p1->GetVertexWeight() -  p2->GetVertexWeight() ); //QIN_nFeatureCnt
}
