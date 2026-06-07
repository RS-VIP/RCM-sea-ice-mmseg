// MarkovNetClassifier.cpp : implementation of the TMarkovNetClassifier class
//
// Markov network classifier. It is similar to Markov random field,
// but with irregular grid
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "MarkovNetClassifier.h"
#include "EstimateUtil.h"

//
// Constructor
// 
TMarkovNetClassifier::TMarkovNetClassifier(int nAnnealingMode, double dbInitT, int nIterCnt,
                                           int nClassCnt)
	: TRegionClassifier()
{
	m_pAnnealUtil = NULL;
	m_pScanUtil = NULL;
	m_nAnnealingMode = nAnnealingMode;
	m_dbInitT = dbInitT;
	m_nIterCnt = nIterCnt;
	m_nClassCnt = nClassCnt;
	m_nClassPixelCnt = NULL;

	//srand((unsigned) time(NULL));
}

//
// Destructor
// 
TMarkovNetClassifier::~TMarkovNetClassifier()
{
	if ( m_pAnnealUtil )
		delete m_pAnnealUtil;
	if ( m_pScanUtil ) 
		delete m_pScanUtil;

	if (m_nClassPixelCnt)
		delete [] m_nClassPixelCnt;

}

// 
// Initialize the classifier
//
bool TMarkovNetClassifier::Init()
{
	if ( m_nClassCnt <= 0 )
		return false;
	m_pScanUtil = new TSampler();
	if (m_pAnnealUtil)
		delete m_pAnnealUtil;
	m_pAnnealUtil = new TAnnealingUtil(m_dbInitT, 0.98); //Here the ratio for descreasing temperature is fixed at 0.9

	m_nClassPixelCnt = new int[m_nClassCnt];
	return true;
}

//
// Intermediate classification of the regions by using simulated annealing
// 
bool TMarkovNetClassifier::IntermediateClassify()
{
	if ( m_nClassCnt <= 0 )
		return false;

	if ( ! m_pAnnealUtil )
		m_pAnnealUtil = new TAnnealingUtil(m_dbInitT, 0.98);

	if (m_nAnnealingMode & SINGLE_STEP)  // one step of the annealing
	{
		Annealing();
		if (!(m_nAnnealingMode & NO_COOLING))
			// compute the next temparature
			m_pAnnealUtil->NextTemperature();
	}
	else  // a full annealing process
	{
		// Refine using simulated annealing
		for ( int i = 0; i < m_nIterCnt; i++ )
		{
			Annealing();
			if (!(m_nAnnealingMode & NO_COOLING))
				// compute the next temparature
				m_pAnnealUtil->NextTemperature();
		}

		delete m_pAnnealUtil;
		m_pAnnealUtil = 0;
	}

	return TRegionClassifier::IntermediateClassify();
}

// 
// Simulated annealing for intermediate classification
// 
// Parameter in
//  pList: the regions in the form of graph vertex linked list
//
// Parameter out
//  pList: the regions with class labels assigned
// 
void TMarkovNetClassifier::Annealing()
{
	if (!m_pGraph)
		return;

	// Get the number of vertice in graph
	int nRegionCnt = m_pGraph->GetVertexCnt();
	TRegion** ppList = new TRegion*[nRegionCnt];
	TRegion* pRegion = (TRegion*) m_pGraph->GetVertices();
	for (int i = nRegionCnt - 1; i >= 0; i --)
	{
		ppList[i] = pRegion;
		pRegion = (TRegion*) pRegion->m_pNext;
	}

	// Annealing
	double *dbEnergies = new double[m_nClassCnt];
	double U;
	m_pScanUtil->InitScan(nRegionCnt);
	int k = m_pScanUtil->GetNextScan();
	while (k >= 0)
	{
		pRegion = ppList[k];
		for (k = 0; k < m_nClassCnt; k ++)
		{
			// Unary property energy
			U = GetUnaryEnergy(pRegion, k);
			// Binary relation energy
			for ( int j = 0; j < pRegion->m_nNeighborCnt; j++ )
				U += GetBinaryEnergy(pRegion, (TBoundary*) pRegion->m_pEdges[j], k);
			dbEnergies[k] = U;  //QIN_nFeatureCnt

		}

		// Sampling
		pRegion->m_nClassLabel = m_pAnnealUtil->GibbsSampler(dbEnergies, m_nClassCnt);
		k = m_pScanUtil->GetNextScan();
	}

	delete[] dbEnergies;  //QIN_MOD	
	delete[] ppList;  //QIN_MOD
}

// 
// Compute the overall energy corresponding to the current configuration
// 
double TMarkovNetClassifier::OverallEnergy()
{	
	double U = 0;

	// Unary property energy
	TRegion* pRegion = (TRegion*) m_pGraph->GetVertices();
	while (pRegion)
	{
		U += GetUnaryEnergy(pRegion, pRegion->m_nClassLabel);
		pRegion = (TRegion*) pRegion->m_pNext;
	}

	// Binary relation energy, note that each boundary is only counted once
	TBoundary *pBoundary = (TBoundary*) m_pGraph->GetEdges();
	while (pBoundary)
	{
		pRegion = (TRegion*) pBoundary->m_pNeighbor1;
		U += GetBinaryEnergy(pRegion, pBoundary, pRegion->m_nClassLabel);
		pBoundary = (TBoundary*) pBoundary->m_pNext;
	}

	return U;
}
