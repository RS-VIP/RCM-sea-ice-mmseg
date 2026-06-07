/** \file
 *
 * \brief TRegionBasedClustering class.
 *
 * Class for algorithms that perform region-based clustering
 *
 * In this file, all the GMM and KMGMM methods are commented because they are
 * theoretically inappropriate for region-based processing (Alex Qin)
 *
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "RegionBasedClustering.h"
#include "EstimateUtil.h"
//#include <random>

#define REG_COEF 0.001

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

//
// Constructor
// 
TRegionBasedClustering::TRegionBasedClustering( TGrayImage<int>* pMap, TMonoImage* pMask,
	                                            TGrayImage<FLOAT>** ppSrcImg, int nClusterCnt, int nFeatureCnt)
	                  : TGaussianClustering(nClusterCnt)
{

	m_pGraph = new TRegionAdjacencyGraph(pMap, pMask, ppSrcImg, nFeatureCnt, 
										TRegionAdjacencyGraph::GAUSSIAN);
	m_pGraph->ConstructGraph(); 	// Construct the graph
	m_bGraphNew = true;

	m_nFeatureCnt = nFeatureCnt;

	// m_nFeatureCnt and m_nClusterCnt must be called before allocating.
	if (m_pGraph)
		Alloc();

	//srand((unsigned) time(NULL));
}

//
// Constructor
// 
TRegionBasedClustering::TRegionBasedClustering( TRegionAdjacencyGraph* pGraph, int nClusterCnt)
					   : TGaussianClustering(nClusterCnt)
{
	m_pGraph = pGraph;
	m_nFeatureCnt = pGraph->GetFeatureCnt();
	m_bGraphNew = false;

	// m_nFeatureCnt and m_nClusterCnt must be called before allocating.
	if (m_pGraph)
		Alloc();

	//srand((unsigned) time(NULL));
}
//
// Destructor
// 
TRegionBasedClustering::~TRegionBasedClustering()
{

    if ( m_bGraphNew )
		delete m_pGraph;
}


void TRegionBasedClustering::RegionInitParam( int nInitMode, double* VecMin, double* VecMax, double* Variance )
{
	if ( nInitMode == RAND1 )
	{
		int nFeatureCnt = m_pGraph->GetFeatureCnt();

		int i, j;
		for ( i = 0; i < m_nClusterCnt; i++ )
		{
			m_ppCovInv[i]->SetZero();
			m_dbPri[i] = 1.0 / m_nClusterCnt;  // assume the priors equal
		}
		
		double dbRange, dbMin;
		for ( i = 0; i < nFeatureCnt; i++ )
		{
			dbMin = VecMin[i];
			dbRange = VecMax[i] - dbMin;
			for ( j = 0; j < m_nClusterCnt; j++ )
			{
				// Choose means randomly
				(*m_ppMu[j])(i) = dbMin + rand() * dbRange / RAND_MAX;

				//// Evenly divide the range to m_nClusterCnt parts, and use half as the
				//// standard deviation. Different features are assumed indepent first
				//(*m_ppCovInv[j])(i, i) = 4 /
				//	((dbRange / m_nClusterCnt) * (dbRange / m_nClusterCnt));

				(*m_ppCovInv[j])(i, i) = 1 / Variance[i];
			}
		}

		double dbDeterm = 1;
		for ( i = 0; i < nFeatureCnt; i++ )
			dbDeterm *= ( 1 / (*m_ppCovInv[0])(i, i) );
		dbDeterm = sqrt(dbDeterm);   // Note it is the square root of the determinant
		for ( i = 0; i < m_nClusterCnt; i++ )
			m_dbDeterm[i] = dbDeterm;
	}

	else if ( nInitMode == RAND2 )
	{
		int nFeatureCnt = m_pGraph->GetFeatureCnt();

		double dbAvgVar = 0.0;
		for ( int i = 0; i < nFeatureCnt; i++ )
			dbAvgVar += Variance[i];
		dbAvgVar /= nFeatureCnt;

		GenerateRandomLabel();
		UpdateParam(dbAvgVar);
	}

	else
	{
		printf("Error: The specified intialization method does not exist!");
		//_getchar_nolock();
		exit(1);
	}
}

// 
// Generate a random class label configuration
//
void TRegionBasedClustering::GenerateRandomLabel()
{
	TRegion* pCur = (TRegion*) m_pGraph->GetVertices();
	int nVertexCnt = m_pGraph->GetVertexCnt();
    int i, nLabel;
	int* nLabelFlag = new int[m_nClusterCnt];
    bool bContinueFlag;
	int nCnt = 0;

	do
	{
		nCnt++;
		memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));
		while (pCur)
		{
            // Evil bug averted: if RAND_MAX is the same as the max integer value, overflow occurs and nearly all the
            // labels are set to zero. Oh the joy of porting old C++ code to linux... by the way there are probably other
            // bugs like this in the codebase so if you feel like your code is haunted that's something to check for.

			// nLabel = rand() * m_nClusterCnt / (RAND_MAX + 1.0);
            nLabel = rand() % m_nClusterCnt;

			pCur->m_nClassLabel = nLabel;
		    if ( nLabelFlag[nLabel] == 0 )
				nLabelFlag[nLabel] = 1;
			pCur = (TRegion*) pCur->m_pNext;
		}
		
		nLabel = 0;
		for ( i = 0; i < m_nClusterCnt; i++ )
			nLabel += nLabelFlag[i];

		if ( nLabel == m_nClusterCnt || nVertexCnt < m_nClusterCnt )
			bContinueFlag = false;
		else
			bContinueFlag = true;
	} while (bContinueFlag && nCnt < 10);

	delete[] nLabelFlag;
}

//
// Update the parameters for each cluster
//
bool TRegionBasedClustering::UpdateParam(double dbRegCoef)
{
	bool bUpdateFlag = true;

	int *nSampleCnt = new int[m_nClusterCnt];
	int i, j, k;
	for ( k = 0; k < m_nClusterCnt; k++ )
	{
		m_ppMu[k]->SetZero();
		m_ppCovInv[k]->SetZero();
		nSampleCnt[k] = 0;
	}

	// Mean and Variance
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	double *dbMu, **dbCoMu;
	TGaussianRegion* pCurrent;
	TVertex* pVertex = m_pGraph->GetVertices();
	while (pVertex)
	{
		pCurrent = (TGaussianRegion*) pVertex;
		k = pCurrent->m_nClassLabel;

		dbMu = m_ppMu[k]->m_ppData[0];
		dbCoMu = m_ppCovInv[k]->m_ppData;
		for ( i = 0; i < nFeatureCnt; i ++)
		{
			dbMu[i] += ( (*pCurrent->m_pMean)(i) * pCurrent->m_nPixelCnt );
			for ( j = i; j < nFeatureCnt; j ++)
				dbCoMu[i][j] += ( (*pCurrent->m_pCoMean)(i, j) * pCurrent->m_nPixelCnt );
		}
		nSampleCnt[k] += pCurrent->m_nPixelCnt;

		pVertex = pVertex->m_pNext;
	}

	int nTotalCnt = 0;
	for ( k = 0; k < m_nClusterCnt; k++ )
	{
		m_nClassPixelCnt[k] = nSampleCnt[k]; // store for this later.

		if (nSampleCnt[k])
		{
			dbMu = m_ppMu[k]->m_ppData[0];
			dbCoMu = m_ppCovInv[k]->m_ppData;
			for ( i = 0; i < nFeatureCnt; i++ )
			{
				dbMu[i] /= nSampleCnt[k];
				for ( j = i; j < nFeatureCnt; j++ )
					dbCoMu[i][j] /= nSampleCnt[k];
			}
			for ( i = 0; i < nFeatureCnt; i++ )
			{
				dbCoMu[i][i] -= dbMu[i] * dbMu[i];
				for ( int j = i + 1; j < nFeatureCnt; j++ )
				{
					dbCoMu[i][j] -= dbMu[i] * dbMu[j];
					dbCoMu[j][i] = dbCoMu[i][j];
				}
			}
		}
		else
		{
			bUpdateFlag = false;
			nSampleCnt[k]++;
		}
		
		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += dbRegCoef;

		m_dbDeterm[k] = sqrt(m_ppCovInv[k]->Determinant());
		m_ppCovInv[k]->Inv(m_ppCovInv[k]);

	    nTotalCnt += nSampleCnt[k];  // get total popuation
	}

	// Priors
	for ( k = 0; k < m_nClusterCnt; k++ )
		m_dbPri[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);

	delete[] nSampleCnt; //QIN_MOD

	return bUpdateFlag;
}

//
// Final update the parameters for each cluster
//
bool TRegionBasedClustering::FinalUpdateParam(TGrayImage<int>* pMap, double dbRegCoef)
{
	bool bUpdateFlag = true;

	int *nSampleCnt = new int[m_nClusterCnt];
	int i, j, k;
	for ( k = 0; k < m_nClusterCnt; k++ )
	{
		m_ppMu[k]->SetZero();
		m_ppCovInv[k]->SetZero();
		nSampleCnt[k] = 0;
	}

	// Mean and Variance
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	double *dbMu, **dbCoMu;
	TGaussianRegion* pCurRegion;
	TVertex* pVertex = m_pGraph->GetVertices();
	while (pVertex)
	{
		pCurRegion = (TGaussianRegion*) pVertex;
		k = pCurRegion->m_nClassLabel;

		dbMu = m_ppMu[k]->m_ppData[0];
		dbCoMu = m_ppCovInv[k]->m_ppData;
		for ( i = 0; i < nFeatureCnt; i++ )
		{
			dbMu[i] += ( (*pCurRegion->m_pMean)(i) * pCurRegion->m_nPixelCnt );
			for ( j = i; j < nFeatureCnt; j++ )
				dbCoMu[i][j] += ( (*pCurRegion->m_pCoMean)(i, j) * pCurRegion->m_nPixelCnt );
		}
		nSampleCnt[k] += pCurRegion->m_nPixelCnt;

		pVertex = pVertex->m_pNext;
	}

	// Assign boundary points to one of their neighboring regions
	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
    TChainCode *pChain;
	POINTex Pt;
	int nDir;
	TBoundary* pCurBoundary = (TBoundary*) m_pGraph->GetEdges();
	while (pCurBoundary)
	{
		pChain = pCurBoundary->m_pChain;
		while (pChain)
		{
			pChain->GetFirstSite(Pt, nDir);
			do {
				k = (*pMap)(Pt.Row, Pt.Col);
				dbMu = m_ppMu[k]->m_ppData[0];
		        dbCoMu = m_ppCovInv[k]->m_ppData;
		        for ( i = 0; i < nFeatureCnt; i++ )
		        {
					dbMu[i] += (*ppFeatures[i])(Pt.Row, Pt.Col);
			        for ( j = i; j < nFeatureCnt; j++ )
						dbCoMu[i][j] += ( (*ppFeatures[i])(Pt.Row, Pt.Col) * (*ppFeatures[j])(Pt.Row, Pt.Col) );
		         }
				  nSampleCnt[k]++;
			} while ( pChain->GetNextSite(Pt, nDir) );

			pChain = pChain->m_pNext;
		}
		pCurBoundary = (TBoundary*) pCurBoundary->m_pNext;
	}

	int nTotalCnt = 0;
	for ( k = 0; k < m_nClusterCnt; k++ )
	{
		m_nClassPixelCnt[k] = nSampleCnt[k]; // Store this for later.
		if (nSampleCnt[k])
		{
			dbMu = m_ppMu[k]->m_ppData[0];
			dbCoMu = m_ppCovInv[k]->m_ppData;
			for ( i = 0; i < nFeatureCnt; i++ )
			{
				dbMu[i] /= nSampleCnt[k];
				for ( j = i; j < nFeatureCnt; j++ )
					dbCoMu[i][j] /= nSampleCnt[k];
			}
			for ( i = 0; i < nFeatureCnt; i++ )
			{
				dbCoMu[i][i] -= dbMu[i] * dbMu[i];
				for ( j = i + 1; j < nFeatureCnt; j++ )
				{
					dbCoMu[i][j] -= dbMu[i] * dbMu[j];
					dbCoMu[j][i] = dbCoMu[i][j];
				}
			}
		}
		else
		{
			bUpdateFlag = false;
			nSampleCnt[k]++;
		}
		
		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += dbRegCoef;

		m_dbDeterm[k] = sqrt(m_ppCovInv[k]->Determinant());
		m_ppCovInv[k]->Inv(m_ppCovInv[k]);

	    nTotalCnt += nSampleCnt[k];  // get total popuation
	}

	// Priors
	for ( k = 0; k < m_nClusterCnt; k++ )
		m_dbPri[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);

	delete[] nSampleCnt; //QIN_MOD
	return bUpdateFlag;
}

// *** Region-based GMM is commented out because it is theoretically inappropriate for region-based processing (Alex Qin)
/// **************************************** Region-based GMM clustering **************************************** ///
// 
// Region-based GMM main routine
// 
// Parameter in
//
// Return
//  the clustering result map
//
//TGrayImage<int>* TRegionBasedClustering::RegionGMMClustering(int nInitMode, bool& bSucFlag, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal)
//{
//	TGrayImage<int>* pRet = 0;
//
//	// Initialize the classifier
//	Alloc(nClusterCnt);
//
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//	TVector **ppMuBest = AllocVectors(m_nClusterCnt, nFeatureCnt);  // for preserving the best
//	TMatrix **ppCovInvBest = AllocMatrice(m_nClusterCnt, nFeatureCnt);  // for preserving the best
//	double* dbPriBest = new double[m_nClusterCnt];  // for preserving the best
//	double* dbDetermBest = new double[m_nClusterCnt];  // for preserving the best
//
//	// Parameter estimation
//	double dbEnergy, dbMinEnergy = -1;
//	int i, nIter = 0;  // iteration number
//		
//	// Find the range and variance of each feature
//	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
//	TMonoImage* pMask = m_pGraph->GetMask();
//	double* VecMin = new double[nFeatureCnt];
//	double* VecMax = new double[nFeatureCnt];
//	double* Mu = new double[nFeatureCnt];
//	double* Variance = new double[nFeatureCnt];
//	for ( i = 0; i < nFeatureCnt; i++ )
//	{
//		((TImage*) ppFeatures[i])->GetIntensityRange(VecMax[i], VecMin[i], pMask);
//		Mu[i] = ppFeatures[i]->GetIntensityMean(pMask);
//	}
//
//	int nCnt = 0;
//	double dbAvgVar = 0.0;
//	for ( i = 0; i < nFeatureCnt; i++ )
//	{
//		Variance[i] = 0;
//		nCnt = 0;
//		for ( int j = 0; j < ppFeatures[i]->Height(); j++ )
//			for ( int k = 0; k < ppFeatures[i]->Width(); k++ )
//			{
//				if ( pMask && ! (*pMask)(j, k) )
//					continue;
//				Variance[i] += ( ( (*ppFeatures[i])(j, k) - Mu[i] ) * ( (*ppFeatures[i])(j, k) - Mu[i] ) );
//				nCnt++;
//			}
//		Variance[i] /= nCnt;
//		dbAvgVar += Variance[i];
//	}
//	dbAvgVar /= nFeatureCnt;
//
//	do {
//		RegionInitParam(nInitMode, VecMin, VecMax, Variance);
//			
//		dbEnergy = RegionGMMIteration(dbRatio, nInit1);   //QIN_20090504: 10 times
//
//		if (dbMinEnergy < 0 ||  // first iteration
//			dbEnergy < dbMinEnergy)  
//		{
//			dbMinEnergy = dbEnergy;
//			for ( i = 0; i < m_nClusterCnt; i++ )
//			{
//				ppMuBest[i]->CopyFrom(m_ppMu[i]);
//				ppCovInvBest[i]->CopyFrom(m_ppCovInv[i]);
//			}
//			memcpy(dbPriBest, m_dbPri, m_nClusterCnt * sizeof(double));
//			memcpy(dbDetermBest, m_dbDeterm, m_nClusterCnt * sizeof(double));
//		}
//		//printf("%d", nIter);
//	} while ( ++nIter < nInit2 );
//
//	// Set the parameters with the obtained best and continue with the EM
//	for ( i = 0; i < m_nClusterCnt; i++ )
//	{
//		m_ppMu[i]->CopyFrom(ppMuBest[i]);
//		m_ppCovInv[i]->CopyFrom(ppCovInvBest[i]);
//	}
//	memcpy(m_dbPri, dbPriBest, m_nClusterCnt * sizeof(double));
//	memcpy(m_dbDeterm, dbDetermBest, m_nClusterCnt * sizeof(double));
//
//	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);
//
//	// Further iterations to converge
//	dbMinEnergy = RegionGMMIteration(dbRatio, nIterFinal);
//
//	//printf("The final minimum energy is %.6f\n", dbMinEnergy);
//
//	for ( i = 0; i < m_nClusterCnt; i++ )
//	{
//		delete ppMuBest[i];
//		delete ppCovInvBest[i];
//	}
//	delete[] ppMuBest;
//	delete[] ppCovInvBest;
//	delete[] dbPriBest;
//	delete[] dbDetermBest;
//	delete[] VecMin;
//	delete[] VecMax;
//	delete[] Mu;
//	delete[] Variance;
//
//	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
//	//SortCluster();
//
//	//// Show the parameters
//	//ShowParam();
//		
//	RegionGMMFinalClassify(bSucFlag);
//
//	pRet = RegionGMMGetClassMap();
//
//	FinalUpdateParam(pRet, dbAvgVar);
//
//	return pRet;
//}

// 
// The iterations for region-based GMM by Expectation Maximization
// 
//double TRegionBasedClustering::RegionGMMIteration(double dbRatio, int nIterations)
//{
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//	int nIter = 0, i, j, k;
//	double dbOldEnergy, dbNewEnergy = -1, dbMax, dbTotal, dbPriTotal;
//	double* dbMu;
//	double** dbCoMu;
//	TVector* VecProb = new TVector(m_nClusterCnt, false);  // probability
//    
//	TVector **ppMuComp = AllocVectors(m_nClusterCnt, nFeatureCnt);
//	TMatrix **ppCovComp = AllocMatrice(m_nClusterCnt, nFeatureCnt);
//	double* dbPri = new double[m_nClusterCnt];  // for the current compuation
//	double* U = new double[m_nClusterCnt];
//	double* dbLnPrior = new double[m_nClusterCnt];
//	double* dbLnSigma = new double[m_nClusterCnt];
//
//	for ( i = 0; i < m_nClusterCnt; i++ )
//	{
//		dbLnPrior[i] = log(m_dbPri[i]);
//		dbLnSigma[i] = log(m_dbDeterm[i]);
//	}
//
//	do {
//		dbOldEnergy = dbNewEnergy;
//		dbNewEnergy = 0;
//		SetZeroVector(ppMuComp, m_nClusterCnt);
//		SetZeroMatrix(ppCovComp, m_nClusterCnt);
//		memset(dbPri, 0, sizeof(double) * m_nClusterCnt);
//
//	    TGaussianRegion* pRegion = (TGaussianRegion*) m_pGraph->GetVertices();
//
//		//if ( nIter == 77 )
//		//	nIter = 77;
//
//	    while (pRegion)
//	    {
//			// Expectation step
//			dbMax = 999999999999;
//			for ( i = 0; i < m_nClusterCnt; i++ )
//			{
//			    dbMu = m_ppMu[i]->m_ppData[0];
//				dbCoMu = m_ppCovInv[i]->m_ppData;
//				U[i] = 0;
//				for ( j = 0; j < nFeatureCnt; j++ )
//				{
//					U[i] += ( ( (*pRegion->m_pCoMean)(j, j) - 2 * (*pRegion->m_pMean)(j) * dbMu[j] + dbMu[j] * dbMu[j] ) * dbCoMu[j][j] );
//		            for ( k = j + 1; k < nFeatureCnt; k++ )
//						U[i] += ( ( (*pRegion->m_pCoMean)(j, k) - (*pRegion->m_pMean)(j) * dbMu[k] - (*pRegion->m_pMean)(k) * dbMu[j] + dbMu[j] * dbMu[k] ) * dbCoMu[j][k] * 2 );
//				}
//				U[i] = U[i] / 2 + dbLnSigma[i];
//				U[i] *= pRegion->m_nPixelCnt;
//				if ( U[i] < dbMax )                  // Here, the minimum value is taken instead of maximum since maximum may result in an infinite number
//					dbMax = U[i];
//			}
//
//			dbTotal = 0;
//			for ( i = 0; i < m_nClusterCnt; i++ )
//			{
//				U[i] -= dbMax;
//				(*VecProb)(i) = exp(-U[i]) * m_dbPri[i];
//				dbTotal += (*VecProb)(i);
//			}
//
//			for ( i = 0; i < m_nClusterCnt; i++ )
//				(*VecProb)(i) = (*VecProb)(i) / dbTotal;
//
//			dbNewEnergy += ( -log(dbTotal) + dbMax );
//
//			// Maximization step
//			for ( i = 0; i < m_nClusterCnt; i++ )
//			{
//				for ( j = 0; j < nFeatureCnt; j++ )
//				{
//					(*ppMuComp[i])(j) += ( (*pRegion->m_pMean)(j) * pRegion->m_nPixelCnt * (*VecProb)(i) );
//					for ( k = j; k < nFeatureCnt; k++ )
//						(*ppCovComp[i])(j, k) += ( (*pRegion->m_pCoMean)(j, k) * pRegion->m_nPixelCnt * (*VecProb)(i) );
//				}
//				dbPri[i] += ( pRegion->m_nPixelCnt * (*VecProb)(i) );
//			}
//
//			pRegion = (TGaussianRegion*) pRegion->m_pNext;
//		}
//
//		if ( nIter == 0 || ( ( nIter < nIterations ) &&( fabs(dbNewEnergy - dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) )
//		{
//			// Compute the total of the priors
//			dbPriTotal = 0;
//			for ( i = 0; i < m_nClusterCnt; i++ )
//			{
//				if ( dbPri[i] <= ZERO_THRESH )
//					dbPri[i] += ZERO_THRESH;
//				dbPriTotal += dbPri[i];
//			}
//
//			// Update the parameters
//			for ( i = 0; i < m_nClusterCnt; i++ )
//			{
//				for ( j = 0; j < nFeatureCnt; j++ )  // mean
//					(*m_ppMu[i])(j) = (*ppMuComp[i])(j) / dbPri[i];
//
//				for ( j = 0; j < nFeatureCnt; j++ )  // covariance
//				{
//					(*ppCovComp[i])(j, j) = (*ppCovComp[i])(j, j) / dbPri[i] - (*m_ppMu[i])(j) * (*m_ppMu[i])(j);
//					for ( k = j + 1; k < nFeatureCnt; k++ )
//					{
//						(*ppCovComp[i])(j, k) = (*ppCovComp[i])(j, k) / dbPri[i] - (*m_ppMu[i])(j) * (*m_ppMu[i])(k);
//                        (*ppCovComp[i])(k, j) = (*ppCovComp[i])(j, k);
//					}
//				}
//
//				if ( ppCovComp[i]->Eigen() <= ZERO_THRESH )
//					for ( j = 0; j < nFeatureCnt; j++ )
//						(*ppCovComp[i])(j, j) += REG_COEF;
//
//				ppCovComp[i]->Inv(m_ppCovInv[i]);  // inverse of covariance
//				m_dbPri[i] = dbPri[i] / dbPriTotal;  // prior
//				m_dbDeterm[i] = sqrt(ppCovComp[i]->Determinant());  // determinant
//			}
//	    }
//
//		nIter++;
//	} while ( ( nIter == 1 ) || ( ( nIter <= nIterations ) && ( fabs(dbNewEnergy - dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) );
//
//	for ( i = 0; i < m_nClusterCnt; i++ )
//	{
//		delete ppMuComp[i];
//		delete ppCovComp[i];
//	}
//	delete[] ppMuComp;
//	delete[] ppCovComp;
//	delete[] dbPri;
//	delete[] U;
//	delete[] dbLnPrior;
//	delete[] dbLnSigma;
//	delete VecProb;
//
//	return dbNewEnergy;
//}

//void TRegionBasedClustering::RegionGMMFinalClassify(bool& bSucFlag)
//{
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//	int i, j, k, nMaxIdx;
//	double dbMax;
//	double* dbMu;
//	double** dbCoMu;
//	double* U = new double[m_nClusterCnt];
//	double dbTmp;
//	double* dbLnSigma = new double[m_nClusterCnt];
//	for ( i = 0; i < m_nClusterCnt; i++ )
//		dbLnSigma[i] = log(m_dbDeterm[i]);
//
//	int* nLabelFlag = new int[m_nClusterCnt];
//	memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));
//
//	TGaussianRegion* pRegion = (TGaussianRegion*) m_pGraph->GetVertices();
//	while (pRegion)
//	{
//		dbMax = 999999999999;
//		for ( i = 0; i < m_nClusterCnt; i++ )
//		{
//			dbMu = m_ppMu[i]->m_ppData[0];
//			dbCoMu = m_ppCovInv[i]->m_ppData;
//			U[i] = 0;
//			for ( j = 0; j < nFeatureCnt; j++ )
//			{
//				U[i] += ( ( (*pRegion->m_pCoMean)(j, j) - 2 * (*pRegion->m_pMean)(j) * dbMu[j] + dbMu[j] * dbMu[j] ) * dbCoMu[j][j] );
//		        for ( k = j + 1; k < nFeatureCnt; k++ )
//					U[i] += ( ( (*pRegion->m_pCoMean)(j, k) - (*pRegion->m_pMean)(j) * dbMu[k] - (*pRegion->m_pMean)(k) * dbMu[j] + dbMu[j] * dbMu[k] ) * dbCoMu[j][k] * 2 );
//			}
//			U[i] = U[i] / 2 + dbLnSigma[i];
//			U[i] *= pRegion->m_nPixelCnt;
//			if ( U[i] < dbMax )
//				dbMax = U[i];
//		}
//
//		for ( i = 0; i < m_nClusterCnt; i++ )
//			U[i] -= dbMax;
//
//		dbMax = -999999999999;
//		for ( i = 0; i < m_nClusterCnt; i++ )
//		{
//			dbTmp = exp(-U[i]) * m_dbPri[i];
//			if ( dbTmp > dbMax )
//			{
//				dbMax = dbTmp;
//				nMaxIdx = i;
//			}
//		}
//
//		pRegion->m_nClassLabel = nMaxIdx;
//
//		if ( nLabelFlag[nMaxIdx] == 0 )
//			nLabelFlag[nMaxIdx] = 1;
//
//		pRegion = (TGaussianRegion*) pRegion->m_pNext;
//	}
//
//	nMaxIdx = 0;
//	for ( i = 0; i < m_nClusterCnt; i++ )
//		nMaxIdx += nLabelFlag[i];
//	if ( nMaxIdx == m_nClusterCnt )
//		bSucFlag = true;
//	else
//		bSucFlag = false;
//
//	delete[] U;
//	delete[] dbLnSigma;
//	delete[] nLabelFlag;
//}

// 
// Get the classification map for region-based GMM clustering
// by mapping each pixel of the image to its class category 
// 
//TGrayImage<int>* TRegionBasedClustering::RegionGMMGetClassMap()
//{
//	TGrayImage<int>* tpMap = m_pGraph->GetMap();
//	TMonoImage* pMask = m_pGraph->GetMask();
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
//
//	// Get maximum number of remaining segments
//	TEstimateUtil EstUtil(0, 0);
//	int nMaxCnt = EstUtil.GetMaxLabel(tpMap, pMask) + 1;
//
//	// Generate the label mapping
//	int* nLabels = new int[nMaxCnt];
//	memset(nLabels, 0, nMaxCnt * sizeof(int));
//	TRegion* pVertex = (TRegion*) m_pGraph->GetVertices();
//	while (pVertex)
//	{
//		nLabels[pVertex->m_nLabel] = pVertex->m_nClassLabel;
//		pVertex = (TRegion*) pVertex->m_pNext;
//	}
//
//	// Obtain the classification map
//	int nWidth = tpMap->Width();
//	int nHeight = tpMap->Height();
//	TGrayImage<int>* pMap = new TGrayImage<int>(tpMap);
//	int l, r, c;
//	for ( r = 0; r < nHeight; r++ )
//		for ( c = 0; c < nWidth; c++ )
//		{
//			if ( pMask && !(*pMask)(r, c) )
//				continue;
//			if ( ( l = (*tpMap)(r, c) ) >= 0 )  // class label
//				(*pMap)(r, c) = nLabels[l];
//		}
//
//	tpMap = new TGrayImage<int>(pMap);
//	// Assign boundary points to one of their neighboring regions
//	int i, j, nRow, nCol, nClassCnt;
//	int* nClass = new int[m_nClusterCnt];
//	TVector* VecFeature = new TVector(nFeatureCnt, false);
//	for ( r = 0; r < nHeight; r++ )
//		for ( c = 0; c < nWidth; c++ )
//		{
//			if ( pMask && !(*pMask)(r, c) )
//				continue;
//			if ( (*tpMap)(r, c) < 0 )
//			{
//				// Check if there are different non-boundary labels in the
//				// neighborhood
//				memset( nClass, 0, m_nClusterCnt * sizeof(int) );
//				nClassCnt = 0;
//				for ( j = 0; j < 8; j++ )
//				{
//					nRow = r + neighbor8[j].Row;
//					nCol = c + neighbor8[j].Col;
//					if ( nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth )
//						continue;
//					if ( pMask && !(*pMask)(nRow, nCol) )  // not valid
//						continue;
//					if ( ( l = (*tpMap)(nRow, nCol) ) < 0 )  // also boundary
//						continue;
//					if ( nClass[l] == 0 )
//					{
//						nClass[l]++;
//					    nClassCnt++;
//					}
//				}
//
//				if ( nClassCnt == 0 )
//				{
//					printf("Error: one boundary point does not belong to any neighoring regions!");
//					_getchar_nolock();
//					exit(1);
//				}
//
//				int nMinIdx;
//				double dbTmp, dbMin = 999999999999;
//				if ( nClassCnt == 1 )
//				{
//					for ( i = 0; i < m_nClusterCnt; i++ )
//					{
//						if ( nClass[i] > 0 )
//						{
//							(*pMap)(r, c) = i;
//							break;
//						}
//					}
//
//				}
//				else if ( nClassCnt > 1 )
//				{
//					for ( i = 0; i < nFeatureCnt; i++ )
//						(*VecFeature)(i) = (*ppFeatures[i])(r, c);
//					for ( i = 0; i < m_nClusterCnt; i++ )
//					{
//						if ( nClass[i] > 0 )
//						{
//							dbTmp = RegionGMMBoundaryClassify( VecFeature, i );
//							if ( dbTmp < dbMin )
//							{
//								dbMin = dbTmp;
//								nMinIdx = i;
//							}
//						}
//					}
//                    (*pMap)(r, c) = nMinIdx;
//				}
//			}
//		}
//
//	delete[] nLabels; //QIN_MOD
//	delete[] nClass;
//	delete tpMap;
//
// 	delete VecFeature;
//
//	return pMap;
//}
//
////
//// Boundary pixels classification, in which boundary pixels are classified into one of their neighoring regions
////
//int TRegionBasedClustering::RegionGMMBoundaryClassify(TVector* VecFeature, int nLabel1, int nLabel2)
//{
//	int i;
//	double dbTotal1, dbTotal2;
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//	TVector* pVecTemp1 = new TVector(nFeatureCnt, false);  // a temperoray vector
//	TVector* pVecTemp2 = new TVector(nFeatureCnt, false);  // a temperoray vector
//
//	for ( i = 0; i < nFeatureCnt; i++ )
//		(*pVecTemp1)(i) = (*VecFeature)(i) - (*m_ppMu[nLabel1])(i);
//	pVecTemp1->Mul(pVecTemp2, m_ppCovInv[nLabel1]);
//	dbTotal1 = m_dbPri[nLabel1] * exp( -pVecTemp2->VecMul(pVecTemp1) / 2) / m_dbDeterm[nLabel1];
//
//	for ( i = 0; i < nFeatureCnt; i++ )
//		(*pVecTemp1)(i) = (*VecFeature)(i) - (*m_ppMu[nLabel2])(i);
//		pVecTemp1->Mul(pVecTemp2, m_ppCovInv[nLabel2]);
//		dbTotal2 = m_dbPri[nLabel2] * exp( -pVecTemp2->VecMul(pVecTemp1) / 2) / m_dbDeterm[nLabel2];
//
//	delete pVecTemp1;
//	delete pVecTemp2;
//
//	if ( dbTotal1 >= dbTotal2 )
//		return nLabel1;
//	else
//		return nLabel2;
//}
//
////
//// Boundary pixels classification, in which boundary pixels are classified into one of their neighoring regions
////
//double TRegionBasedClustering::RegionGMMBoundaryClassify(TVector* VecFeature, int nLabel)
//{
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//	TVector* pVecTemp1 = new TVector(nFeatureCnt, false);  // a temperoray vector
//	TVector* pVecTemp2 = new TVector(nFeatureCnt, false);  // a temperoray vector
//
//	for ( int i = 0; i < nFeatureCnt; i++ )
//		(*pVecTemp1)(i) = (*VecFeature)(i) - (*m_ppMu[nLabel])(i);
//	pVecTemp1->Mul(pVecTemp2, m_ppCovInv[nLabel]);
//	
//	delete pVecTemp1;
//	delete pVecTemp2;
//
//	return -1 * m_dbPri[nLabel] * exp( -pVecTemp2->VecMul(pVecTemp1) / 2) / m_dbDeterm[nLabel];
//
//}

/// **************************************** Region-based KM clustering **************************************** ///
/**
 * 
 * Param in
 * nInitMode: Initialization mode to use.
 *
 * dbRatio: Parameter for the initialization algorithm, controls the
 * stopping criterion.
 *
 * nInit1: Parameter for the initialization algorithm, the number of
 * iterations of the initialization clustering algorithm.
 * nInit2: Parameter for the initialization algorithm, how many times
 * the clustering algorithm is executed
 *
 * The K-means clustering algorithm will be run nInit2 times,
 * each time starting from a new random starting point. Each time, the
 * K-means algorithm will be run for nInit1 iterations. Out of the
 * nInit2 times, the result with the lowest energy is chosen.
 *
 * nIterFinal : After the lowest energy initilization is chosen, the
 * K-means algorithm is further run for this many interations
 * for further convergence.
 *
 * Return
 * 	the clustering result map
 *
 */
TGrayImage<int>* TRegionBasedClustering::RegionKMClustering(int nInitMode, bool& bSucFlag, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	TGrayImage<int>* pRet = 0;

	// Initialize the classifier
	//Alloc(m_nClusterCnt);

	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	TVector **ppMuBest = AllocVectors(m_nClusterCnt, nFeatureCnt);  // for preserving the best

	// Parameter estimation
	double dbEnergy, dbMinEnergy = -1;
	int i, j, k, nIter = 0;  // iteration number
	
	// Find the range and variance of each feature
	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
	TMonoImage* pMask = m_pGraph->GetMask();
	double* VecMin = new double[nFeatureCnt];
	double* VecMax = new double[nFeatureCnt];
	double* Mu = new double[nFeatureCnt];
	double* Variance = new double[nFeatureCnt];

	for ( i = 0; i < nFeatureCnt; i++ )
	{
		((TImage*) ppFeatures[i])->GetIntensityRange(VecMax[i], VecMin[i], pMask);
		Mu[i] = ppFeatures[i]->GetIntensityMean(pMask);
	}

	int nCnt = 0;
	double dbAvgVar = 0.0;
	for ( i = 0; i < nFeatureCnt; i++ )
	{
		Variance[i] = 0;
		nCnt = 0;
		for ( j = 0; j < ppFeatures[i]->Height(); j++ )
			for ( k = 0; k < ppFeatures[i]->Width(); k++ )
			{
				if ( pMask && ! (*pMask)(j, k) )
					continue;
				Variance[i] += ( ( (*ppFeatures[i])(j, k) - Mu[i] ) * ( (*ppFeatures[i])(j, k) - Mu[i] ) );
				nCnt++;
			}
		Variance[i] /= nCnt;
		dbAvgVar += Variance[i];
		//System::Diagnostics::Debug::Writeline("Variance" + i +" = " + Variance[i]);
		//Debug::Writeline("dbAvgVar" + dbAvgVar);
	}
	dbAvgVar /= nFeatureCnt;

	do {
		RegionInitParam(nInitMode, VecMin, VecMax, Variance);

		dbEnergy = RegionKMIteration(dbRatio, nInit1); //QIN_20090504: 10 times?

		if (dbMinEnergy < 0 ||  // first iteration
			dbEnergy < dbMinEnergy)  
		{
			dbMinEnergy = dbEnergy;
			for ( i = 0; i < m_nClusterCnt; i++ )
				ppMuBest[i]->CopyFrom(m_ppMu[i]);
		}
		//printf("%d", nIter);
	} while ( ++nIter < nInit2 );

	// Set the parameters with the obtained best and continue with the EM
	for (i = 0; i < m_nClusterCnt; i ++)
		m_ppMu[i]->CopyFrom(ppMuBest[i]);

	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);

	// Further iterations to converge
	dbMinEnergy = RegionKMIteration(dbRatio, nIterFinal);
		
	//printf("The final minimum energy is %.6f\n", dbMinEnergy);

	for (i = 0; i < m_nClusterCnt; i ++)
		delete ppMuBest[i];
	delete[] ppMuBest;

	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	//SortCluster();

	//// Show the parameters
	//ShowParam();
		
	RegionKMFinalClassify(bSucFlag);

	pRet = RegionKMGetClassMap();       // This step is not needed in paper since boundary is not considered 20090408

	FinalUpdateParam(pRet, dbAvgVar);

	return pRet;
}

// 
// The iterations for region-based KM
// 
double TRegionBasedClustering::RegionKMIteration(double dbRatio, int nIterations)
{
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	int nIter = 0, i, j, nMinIdx;
	double dbOldEnergy, dbNewEnergy = -1, U, dbMin;
	double* dbMu;

	TVector **ppMuComp = AllocVectors(m_nClusterCnt, nFeatureCnt);
	int* nSampleCnt = new int[m_nClusterCnt];
	do {
		dbOldEnergy = dbNewEnergy;
		dbNewEnergy = 0;
		SetZeroVector( ppMuComp, m_nClusterCnt );
		memset( nSampleCnt, 0, m_nClusterCnt*sizeof(int) );

	    TGaussianRegion* pRegion = (TGaussianRegion*) m_pGraph->GetVertices();
	    while (pRegion)
	    {
			dbMin = 999999999999e9;
			for ( i = 0; i < m_nClusterCnt; i++ )
			{
			    dbMu = m_ppMu[i]->m_ppData[0];
				U = 0;
				for ( j = 0; j < nFeatureCnt; j++ )
					U += ( (*pRegion->m_pCoMean)(j, j) - 2 * (*pRegion->m_pMean)(j) * dbMu[j] + dbMu[j] * dbMu[j] );
				if ( U < dbMin )
				{
					dbMin = U;
					nMinIdx = i;
				}
			}

			for ( i = 0; i < nFeatureCnt; i++ )
				(*ppMuComp[nMinIdx])(i) += ( (*pRegion->m_pMean)(i) * pRegion->m_nPixelCnt );

			nSampleCnt[nMinIdx] += pRegion->m_nPixelCnt;

			dbNewEnergy += ( dbMin * pRegion->m_nPixelCnt );

			pRegion = (TGaussianRegion*) pRegion->m_pNext;
		}

		if ( nIter == 0 || ( ( nIter < nIterations ) &&( fabs(dbNewEnergy -dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) )
		{
			for ( i = 0; i < m_nClusterCnt; i++ )
				for ( j = 0; j < nFeatureCnt; j++ )  // mean
				{
					if ( nSampleCnt[i] > 0 )
						(*m_ppMu[i])(j) = (*ppMuComp[i])(j) / nSampleCnt[i];
				}
		}

		nIter++;

	} while ( ( nIter == 1 ) || ( ( nIter <= nIterations ) && ( fabs(dbNewEnergy -dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) );


	for ( i = 0; i < m_nClusterCnt; i++ )
		delete ppMuComp[i];
	delete[] ppMuComp;
	delete[] nSampleCnt;
	return dbNewEnergy;
}

void TRegionBasedClustering::RegionKMFinalClassify(bool& bSucFlag)
{
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	int i, j, nMinIdx;
	double U, dbMin;
	double* dbMu;
	int* nLabelFlag = new int[m_nClusterCnt];
	memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));

	TGaussianRegion* pRegion = (TGaussianRegion*) m_pGraph->GetVertices();
	while (pRegion)
	{
		dbMin = 999999999999;
		for (i = 0; i < m_nClusterCnt; i ++)
		{
			dbMu = m_ppMu[i]->m_ppData[0];
			U = 0;
			for ( j = 0; j < nFeatureCnt; j++ )
				U += ( ( (*pRegion->m_pCoMean)(j, j) - 2 * (*pRegion->m_pMean)(j) * dbMu[j] + dbMu[j] * dbMu[j] ) );
			if ( U < dbMin )
			{
				dbMin = U;
				nMinIdx = i;
			}
		}
		pRegion->m_nClassLabel = nMinIdx;

		if ( nLabelFlag[nMinIdx] == 0 )
			nLabelFlag[nMinIdx] = 1;

		pRegion = (TGaussianRegion*) pRegion->m_pNext;
	}
	
	nMinIdx = 0;
	for ( i = 0; i < m_nClusterCnt; i++ )
		nMinIdx += nLabelFlag[i];
	if ( nMinIdx == m_nClusterCnt )
		bSucFlag = true;
	else
		bSucFlag = false;

	delete[] nLabelFlag;
}

// 
// Get the classification map for region-based KM clustering
// by mapping each pixel of the image to its class category 
// 
TGrayImage<int>* TRegionBasedClustering::RegionKMGetClassMap()
{
	TGrayImage<int>* tpMap = m_pGraph->GetMap();
	TMonoImage* pMask = m_pGraph->GetMask();
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();

	// Get maximum number of remaining segments
	TEstimateUtil EstUtil;
	int nMaxCnt = EstUtil.GetMaxLabel(tpMap, pMask) + 1;

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

	tpMap = new TGrayImage<int>(pMap);
	// Assign boundary points to one of their neighboring regions
	int i, j, nRow, nCol, nClassCnt;
	int* nClass = new int[m_nClusterCnt];
	TVector* VecFeature = new TVector(nFeatureCnt, false);
	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;
			if ( (*tpMap)(r, c) < 0 )
			{
				// Check if there are different non-boundary labels in the
				// neighborhood
				memset( nClass, 0, m_nClusterCnt * sizeof(int) );
				nClassCnt = 0;
				for ( j = 0; j < 8; j++ )
				{
					nRow = r + neighbor8[j].Row;
					nCol = c + neighbor8[j].Col;
					if ( nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth )
						continue;
					if ( pMask && !(*pMask)(nRow, nCol) )  // not valid
						continue;
					if ( ( l = (*tpMap)(nRow, nCol) ) < 0 )  // also boundary
						continue;
					if ( nClass[l] == 0 )
					{
						nClass[l]++;
					    nClassCnt++;
					}
				}

				if ( nClassCnt == 0 )
					continue;
				//{
				//	printf("Error: one boundary point does not belong to any neighoring regions!");
				//	_getchar_nolock();
				//	exit(1);
				//}

				int nMinIdx;
				double dbTmp, dbMin = 999999999999;
				if ( nClassCnt == 1 )
				{
					for ( i = 0; i < m_nClusterCnt; i++ )
					{
						if ( nClass[i] > 0 )
						{
							(*pMap)(r, c) = i;
							break;
						}
					}

				}
				else if ( nClassCnt > 1 )
				{
					for ( i = 0; i < nFeatureCnt; i++ )
						(*VecFeature)(i) = (*ppFeatures[i])(r, c);
					for ( i = 0; i < m_nClusterCnt; i++ )
					{
						if ( nClass[i] > 0 )
						{
							dbTmp = RegionKMBoundaryClassify( VecFeature, i );
							if ( dbTmp < dbMin )
							{
								dbMin = dbTmp;
								nMinIdx = i;
							}
						}
					}
                    (*pMap)(r, c) = nMinIdx;
				}
			}
		}

	delete[] nLabels;
	delete[] nClass;
	delete tpMap;

 //   TChainCode *pChain;
	//POINT Pt;
	//int nDir, nLabel1, nLabel2, i;
	//TVector* VecFeature = new TVector(nFeatureCnt, false);

	//TBoundary* pCur = (TBoundary*) m_pGraph->GetEdges();
	//while (pCur)
	//{
	//	nLabel1 = ((TRegion*) pCur->m_pNeighbor1)->m_nClassLabel;
	//	nLabel2 = ((TRegion*) pCur->m_pNeighbor2)->m_nClassLabel;
	//	pChain = pCur->m_pChain;
	//	while (pChain)
	//	{
	//		pChain->GetFirstSite(Pt, nDir);
	//		do {
 //               if ( nLabel1 == nLabel2 )
	//				(*pMap)(Pt.Row, Pt.Col) = nLabel1;
	//			else
	//			{
	//				for ( i = 0; i < nFeatureCnt; i++ )
	//					(*VecFeature)(i) = (*ppFeatures[i])(Pt.Row, Pt.Col);
	//				
	//				(*pMap)(Pt.Row, Pt.Col) = RegionKMBoundaryClassify( VecFeature, nLabel1, nLabel2 );
	//			}

	//		} while (pChain->GetNextSite(Pt, nDir));

	//		pChain = pChain->m_pNext;
	//	}
	//	pCur = (TBoundary*) pCur->m_pNext;
	//}

	//TGrayImage<int>* pDispMap = new TGrayImage<int>(nWidth, nHeight);
	//for (int r = 0; r < nHeight; r ++)
	//	for (int c = 0; c < nWidth; c ++)
	//	{
	//		if ( ( l = (*pMap)(r, c) ) >= 0 )  // class label
	//		{
	//			if (l == 0)
	//				(*pDispMap)(r, c) = 128;
	//			if (l == 1)
	//				(*pDispMap)(r, c) = 0;
	//			if (l == 2)
	//				(*pDispMap)(r, c) = 255;
	//			if (l == 3)
	//				(*pDispMap)(r, c) = 64;
	//			if (l == 4)
	//				(*pDispMap)(r, c) = 192;
	//			if (l == 5)
	//				(*pDispMap)(r, c) = 32;
	//			if (l == 6)
	//				(*pDispMap)(r, c) = 224;
	//		}
	//		else  // boundary label
	//			(*pDispMap)(r, c) = l;
	//	}
	//
	//pDispMap->Save("KM_Region_Result.bmp");
	//delete pDispMap;

	delete VecFeature;

	return pMap;
}

//
// Boundary pixels classification, in which boundary pixels are classified into one of their neighoring regions
//
//int TRegionBasedClustering::RegionKMBoundaryClassify(TVector* VecFeature, int nLabel1, int nLabel2)
//{
//	int i;
//	double dbTotal1 = 0, dbTotal2 = 0;
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//
//	for ( i = 0; i < nFeatureCnt; i++ )
//		dbTotal1 += ( ( (*VecFeature)(i) - (*m_ppMu[nLabel1])(i) ) * ( (*VecFeature)(i) - (*m_ppMu[nLabel1])(i) ) );
//
//	for ( i = 0; i < nFeatureCnt; i++ )
//		dbTotal2 += ( ( (*VecFeature)(i) - (*m_ppMu[nLabel2])(i) ) * ( (*VecFeature)(i) - (*m_ppMu[nLabel2])(i) ) );
//
//	if ( dbTotal1 >= dbTotal2 )
//		return nLabel1;
//	else
//		return nLabel2;
//}

//
// Boundary pixels classification, in which boundary pixels are classified into one of their neighoring regions
//
double TRegionBasedClustering::RegionKMBoundaryClassify(TVector* VecFeature, int nLabel)
{
	double dbTotal = 0;
	int nFeatureCnt = m_pGraph->GetFeatureCnt();

	for ( int i = 0; i < nFeatureCnt; i++ )
		dbTotal += ( ( (*VecFeature)(i) - (*m_ppMu[nLabel])(i) ) * ( (*VecFeature)(i) - (*m_ppMu[nLabel])(i) ) );

	return dbTotal;
}

/// **************************************** Region-based KM + GMM clustering **************************************** ///
// 
// Region-based KM + GMM main routine
// 
// Parameter in
//
// Return
//  the clustering result map
//
//TGrayImage<int>* TRegionBasedClustering::RegionKMGMMClustering(int nInitMode, bool& bSucFlag, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal)
//{
//	TGrayImage<int>* pRet = 0;
//
//	// Initialize the classifier
//	bool sucFlag = Alloc(nClusterCnt);
//
//	int nFeatureCnt = m_pGraph->GetFeatureCnt();
//	TVector **ppMuBest = AllocVectors(m_nClusterCnt, nFeatureCnt);  // for preserving the best
//
//	// Parameter estimation
//	double dbEnergy, dbMinEnergy = -1;
//	int i, nIter = 0;  // iteration number
//
//	// Find the range and variance of each feature
//	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
//	TMonoImage* pMask = m_pGraph->GetMask();
//	double* VecMin = new double[nFeatureCnt];
//	double* VecMax = new double[nFeatureCnt];
//	double* Mu = new double[nFeatureCnt];
//	double* Variance = new double[nFeatureCnt];
//	for ( i = 0; i < nFeatureCnt; i++ )
//	{
//		((TImage*) ppFeatures[i])->GetIntensityRange(VecMax[i], VecMin[i], pMask);
//		Mu[i] = ppFeatures[i]->GetIntensityMean(pMask);
//	}
//
//	int nCnt = 0;
//	double dbAvgVar = 0.0;
//	for ( i = 0; i < nFeatureCnt; i++ )
//	{
//		Variance[i] = 0;
//		nCnt = 0;
//		for ( int j = 0; j < ppFeatures[i]->Height(); j++ )
//			for ( int k = 0; k < ppFeatures[i]->Width(); k++ )
//			{
//				if ( pMask && ! (*pMask)(j, k) )
//					continue;
//				Variance[i] += ( ( (*ppFeatures[i])(j, k) - Mu[i] ) * ( (*ppFeatures[i])(j, k) - Mu[i] ) );
//				nCnt++;
//			}
//		Variance[i] /= nCnt;
//		dbAvgVar += Variance[i];
//	}
//	dbAvgVar /= nFeatureCnt;
//
//	do {
//		RegionInitParam(nInitMode, VecMin, VecMax, Variance);
//			
//		dbEnergy = RegionKMIteration(dbRatio, nInit1);  //QIN_20090504: 10 times?
//
//		if (dbMinEnergy < 0 ||  // first iteration
//			dbEnergy < dbMinEnergy)  
//		{
//			dbMinEnergy = dbEnergy;
//			for ( i = 0; i < m_nClusterCnt; i++ )
//				ppMuBest[i]->CopyFrom(m_ppMu[i]);
//		}
//		//printf("%d", nIter);
//	} while ( ++nIter < nInit2 );
//
//	// Set the parameters with the obtained best and continue with the EM
//	for (i = 0; i < m_nClusterCnt; i ++)
//	{
//		m_ppMu[i]->CopyFrom(ppMuBest[i]);
//		delete ppMuBest[i];
//	}	
//	delete[] ppMuBest;
//
//	RegionKMFinalClassify(bSucFlag);
//
//	//RegionKMGetClassMap();
//
//	// Update class parameters
//	// *****************************************************************//
//	int *nSampleCnt = new int[m_nClusterCnt];
//	int k;
//	for (k = 0; k < m_nClusterCnt; k ++)
//	{
//		m_ppMu[k]->SetZero();
//		m_ppCovInv[k]->SetZero();
//		nSampleCnt[k] = 0;
//	}
//
//	// Mean and Variance
//	int j;
//	double *dbMu, **dbCoMu;
//	TGaussianRegion* pCurRegion;
//	TVertex* pVertex = m_pGraph->GetVertices();
//	while (pVertex)
//	{
//		pCurRegion = (TGaussianRegion*) pVertex;
//		k = pCurRegion->m_nClassLabel;
//
//		dbMu = m_ppMu[k]->m_ppData[0];
//		dbCoMu = m_ppCovInv[k]->m_ppData;
//		for ( i = 0; i < nFeatureCnt; i++ )
//		{
//			dbMu[i] += ( (*pCurRegion->m_pMean)(i) * pCurRegion->m_nPixelCnt );
//			for ( j = i; j < nFeatureCnt; j++ )
//				dbCoMu[i][j] += ( (*pCurRegion->m_pCoMean)(i, j) * pCurRegion->m_nPixelCnt );
//		}
//		nSampleCnt[k] += pCurRegion->m_nPixelCnt;
//
//		pVertex = pVertex->m_pNext;
//	}
//
//	int nTotalCnt = 0;
//	for ( k = 0; k < m_nClusterCnt; k++ )
//	{
//		if ( nSampleCnt[k] > 0)
//		{
//			dbMu = m_ppMu[k]->m_ppData[0];
//			dbCoMu = m_ppCovInv[k]->m_ppData;
//			for ( i = 0; i < nFeatureCnt; i++ )
//			{
//				dbMu[i] /= nSampleCnt[k];
//				for ( j = i; j < nFeatureCnt; j++ )
//					dbCoMu[i][j] /= nSampleCnt[k];
//			}
//			for ( i = 0; i < nFeatureCnt; i++ )
//			{
//				dbCoMu[i][i] -= dbMu[i] * dbMu[i];
//				for ( j = i + 1; j < nFeatureCnt; j++ )
//				{
//					dbCoMu[i][j] -= dbMu[i] * dbMu[j];
//					dbCoMu[j][i] = dbCoMu[i][j];
//				}
//			}
//		}
//		else
//			nSampleCnt[k]++;
//		
//		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
//			for ( i = 0; i < nFeatureCnt; i++ )
//				(*m_ppCovInv[k])(i, i) += dbAvgVar;
//		m_dbDeterm[k] = sqrt(m_ppCovInv[k]->Determinant());
//		m_ppCovInv[k]->Inv(m_ppCovInv[k]);
//
//		nTotalCnt += nSampleCnt[k];  // get total popuation
//	}
//
//	// Priors
//	for ( k = 0; k < m_nClusterCnt; k++ )
//		m_dbPri[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);
//
//	delete[] nSampleCnt;
//	// *****************************************************************//
//
//	// Further iterations to converge
//	dbMinEnergy = RegionGMMIteration(dbRatio, nIterFinal);
//
//	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
//	//SortCluster();
//
//	//// Show the parameters
//	//ShowParam();
//		
//	RegionGMMFinalClassify(bSucFlag);
//
//	pRet = RegionGMMGetClassMap();
//
//	FinalUpdateParam(pRet, dbAvgVar);
//
//	return pRet;
//}

/// Get segmentation statistics object that lists all the class statistics.
TSegmentationStatistics* TRegionBasedClustering::GetSegStats(TGrayImage<int>* pMap)
{

	

	if (pMap)
	{
		// do an update based on the map.
		FinalUpdateParam(pMap, ZERO_THRESH);
	}

/*	TGaussianSegmentationStatistics tscc(m_nClusterCnt);
	
	TMatrix* temp = 0;

	for (int k = 0; k < m_nClusterCnt; k++)
	{
		tscc.m_vMu[k] = *m_ppMu[k];
		temp = m_ppCovInv[k]->Inv();
		tscc.m_vCov[k] = *temp;
		delete temp;
		tscc.m_vPixelCounts[k] = m_nClassPixelCnt[k];
	}
	
	//return new TGaussianSegmentationStatistics(tscc);
	return tscc.Clone();*/

	return TGaussianClustering::GetSegStats();
}
