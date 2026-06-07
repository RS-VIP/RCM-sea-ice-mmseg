// PixelBasedClustering.cpp:  Implementation of the TPixelBasedClustering class
//
// Class for algorithms that perform pixel-based clustering
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "EstimateUtil.h"
#include "PixelBasedClustering.h"

#define Pi 3.1415926535897
#define REG_COEF 0.001

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

#define BETA_CNT 23
const static double betaLUT[BETA_CNT][2] = {
	{0.6410440, 0.05},
	{0.5941100, 0.15},
	{0.5617630, 0.20}, 
	{0.5195670, 0.25}, 
	{0.4613010, 0.30}, 
	{0.3728720, 0.35}, 
	{0.2246080, 0.40}, 
	{0.1375060, 0.45}, 
	{0.0928173, 0.50}, 
	{0.0663682, 0.55}, 
	{0.0503680, 0.60}, 
	{0.0397756, 0.65}, 
	{0.0327342, 0.70}, 
	{0.0280129, 0.75}, 
	{0.0248621, 0.80}, 
	{0.0228758, 0.85}, 
	{0.0212809, 0.90}, 
	{0.0188605, 1.00}, 
	{0.0175266, 1.15}, 
	{0.0163295, 1.40}, 
	{0.0157941, 1.70}, 
	{0.0153300, 2.0}, 
	{0.0147133, 3.0}, 
//	{0.0135030, 5.0},
};

//
// Constructor for PixelBasedClustering
// 
TPixelBasedClustering::TPixelBasedClustering(TMonoImage* pMask, 
											 TGrayImage<FLOAT>** ppSrcImg, 
											 int nClusterCnt,
											 int nFeatureCnt)
               :       TGaussianClustering(nClusterCnt)
{
	//srand((unsigned) time(NULL));

	if (m_nClusterCnt > 0 && nFeatureCnt > 0)
	{
		m_nFeatureCnt = nFeatureCnt;
		Alloc();
	}

	m_pMask = pMask;
	m_ppSrcImg = ppSrcImg;

}


//
// Destructor
// 
TPixelBasedClustering::~TPixelBasedClustering()
{

}

////
//// Allocate memories needed
//// 
//bool TPixelBasedClustering::Alloc(int nClusterCnt, int nFeatureCnt)
//{
//	if ( nClusterCnt <= 0 || nFeatureCnt <= 0 )
//		return false;
//	//if ( m_nClusterCnt == nClusterCnt && m_nFeatureCnt == nFeatureCnt )
//	//	return true;
//
//	//Clear();
//
//	//m_nClusterCnt = nClusterCnt;
//	//m_nFeatureCnt = nFeatureCnt;
//
//	m_ppMu = AllocVectors(nClusterCnt, nFeatureCnt);
//	m_ppCovInv = AllocMatrice(nClusterCnt, nFeatureCnt);
//	m_dbDeterm = new double[nClusterCnt];
//	m_dbPri = new double[nClusterCnt];
//
//	return true;
//}

void TPixelBasedClustering::GenerateRandomMap(TGrayImage<int>* pMap)
{
	int i, j, nLabel;
	int* nLabelFlag = new int[m_nClusterCnt];
    bool bContinueFlag;
	int nCnt = 0;
	int nPixelCnt = 0;
	for ( i = 0; i < pMap->Height(); i++ )
		for ( j = 0; j <  pMap->Width(); j++ )
		{
			if ( m_pMask && !(*m_pMask)(i, j) )
				continue;
			nPixelCnt++;
		}
	do
	{
		memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));
		for ( i = 0; i < pMap->Height(); i++ )
			for ( j = 0; j <  pMap->Width(); j++ )
			{
				if ( m_pMask && !(*m_pMask)(i, j) )
					continue;
				nLabel = rand() * m_nClusterCnt / (RAND_MAX + 1.0);
				(*pMap)(i, j) = nLabel;
				if ( nLabelFlag[nLabel] == 0 )
					nLabelFlag[nLabel] = 1;
			}
		
		nLabel = 0;
		for ( i = 0; i < m_nClusterCnt; i++ )
			nLabel += nLabelFlag[i];

		if ( nLabel == m_nClusterCnt || nPixelCnt < m_nClusterCnt )
			bContinueFlag = false;
		else
			bContinueFlag = true;

	} while ( bContinueFlag && nCnt < 10 );
}

void TPixelBasedClustering::GenerateRandomHist(int* nHistLabel, int nHistBinCnt)
{
	int i, nLabel;
	int* nLabelFlag = new int[m_nClusterCnt];
    bool bContinueFlag;
	do
	{
		memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));
		for ( i = 0; i < nHistBinCnt - 1; i++ )
		{
			nLabel = rand() * m_nClusterCnt / (RAND_MAX + 1.0);
			nHistLabel[i] = nLabel;
			if ( nLabelFlag[nLabel] == 0 )
				nLabelFlag[nLabel] = 1;
		}
		nLabel = 0;
		for ( i = 0; i < m_nClusterCnt; i++ )
			nLabel += nLabelFlag[i];

		if ( nLabel == m_nClusterCnt )
			bContinueFlag = false;
		else
			bContinueFlag = true;

	} while ( bContinueFlag );
}

void TPixelBasedClustering::PixelInitParam( int nInitMode, double* VecMin, double* VecMax, double* Variance )
{
	// Alloc(m_nClusterCnt, m_nFeatureCnt);
	if ( nInitMode == RAND1 )
	{
		int i, j;
		for ( i = 0; i < m_nClusterCnt; i++ )
		{
			m_ppCovInv[i]->SetZero();
			m_dbPri[i] = 1.0 / m_nClusterCnt;  // assume the priors equal
		}
		
		double dbRange, dbMin;
		for ( i = 0; i < m_nFeatureCnt; i++ )
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
		for ( i = 0; i < m_nFeatureCnt; i++ )
			dbDeterm *= ( 1 / (*m_ppCovInv[0])(i, i) );
		dbDeterm = sqrt(dbDeterm);   // Note it is the square root of the determinant
		for ( i = 0; i < m_nClusterCnt; i++ )
			m_dbDeterm[i] = dbDeterm;
	}

	else if ( nInitMode == RAND2 )
	{
		double dbAvgVar = 0.0;
		for ( int i = 0; i < m_nFeatureCnt; i++ )
			dbAvgVar += Variance[i];
		dbAvgVar /= m_nFeatureCnt;

		int nHeight = m_ppSrcImg[0]->Height();
		int nWidth = m_ppSrcImg[0]->Width();

		TGrayImage<int>* pMap = new TGrayImage<int>( nWidth, nHeight );
		GenerateRandomMap(pMap);
		UpdateParam(pMap, dbAvgVar);
		delete pMap;
	}

	else
	{
		cout << "Error: The specified pixel-based clustering initializaiton mode does not exit!\n";
		//_getchar_nolock();
		exit(1);
	}
}

void TPixelBasedClustering::PixelInitParam( int nInitMode, int* pHist, int nHistBinCnt, double* Variance ) // For fast mode
{
	//Alloc(m_nClusterCnt, 1);
	if ( nInitMode == RAND1 )
	{
		int i, j;
		for ( i = 0; i < m_nClusterCnt; i++ )
		{
			m_ppCovInv[i]->SetZero();
			m_dbPri[i] = 1.0 / m_nClusterCnt;  // assume the priors equal
		}
		
		double dbRange;
		for ( i = 0; i < 1; i++ )
		{
			dbRange = nHistBinCnt - 1;
			for ( j = 0; j < m_nClusterCnt; j++ )
			{
				// Choose means randomly
				(*m_ppMu[j])(i) = 0 + rand() * dbRange / RAND_MAX;

				//// Evenly divide the range to m_nClusterCnt parts, and use half as the
				//// standard deviation. Different features are assumed indepent first
				//(*m_ppCovInv[j])(i, i) = 4 /
				//	((dbRange / m_nClusterCnt) * (dbRange / m_nClusterCnt));

				(*m_ppCovInv[j])(i, i) = 1 / Variance[i];
			}
		}

		double dbDeterm = 1;
		for ( i = 0; i < 1; i++ )
			dbDeterm *= ( 1 / (*m_ppCovInv[0])(i, i) );
		dbDeterm = sqrt(dbDeterm);   // Note it is the square root of the determinant
		for ( i = 0; i < m_nClusterCnt; i++ )
			m_dbDeterm[i] = dbDeterm;
	}

	else if ( nInitMode == RAND2 )
	{
		int* nHistLabel = new int[nHistBinCnt - 1];
		GenerateRandomHist(nHistLabel, nHistBinCnt);

		SetZeroVector(m_ppMu, m_nClusterCnt);
		SetZeroMatrix(m_ppCovInv, m_nClusterCnt);

		int i, k;
		double dbAvgVar = 0.0;
		for ( i = 0; i < 1; i++ )
			dbAvgVar += Variance[i];

		int* nSampleCnt = new int[m_nClusterCnt];
		memset(nSampleCnt, 0, m_nClusterCnt * sizeof(int));
		for ( i = 0; i < nHistBinCnt - 1; i++ )
		{
			k = nHistLabel[i];
			(*m_ppMu[k])(0) += ( pHist[i] * i );
			(*m_ppCovInv[k])(0, 0) += ( pHist[i] * i * i );
			nSampleCnt[k] += pHist[i];
		}

        int nTotalCnt = 0;
		for ( i = 0; i < m_nClusterCnt; i++ )
		{
			if ( nSampleCnt[i] > 0 )
			{
				(*m_ppMu[i])(0) /= nSampleCnt[i];
				(*m_ppCovInv[i])(0, 0) = (*m_ppCovInv[i])(0, 0) / nSampleCnt[i] - (*m_ppMu[i])(0) * (*m_ppMu[i])(0);
			}
			else
				nSampleCnt[i]++;

			if ( m_ppCovInv[i]->Eigen() <= ZERO_THRESH )
				for ( int j = 0; j < m_nClusterCnt; j++ )
					(*m_ppCovInv[i])(j, j) += dbAvgVar;
			m_dbDeterm[i] = sqrt((*m_ppCovInv[i])(0, 0));
			(*m_ppCovInv[i])(0, 0) = 1 / (*m_ppCovInv[i])(0, 0);

			nTotalCnt += nSampleCnt[i];
		}
        
		for ( i = 0; i < m_nClusterCnt; i++ )
			m_dbPri[i] = ((double) nSampleCnt[i]) / ((double) nTotalCnt);

		delete[] nHistLabel;
		delete[] nSampleCnt;
	}

	else
	{
		cout << "Error: The specified pixel-based clustering initializaiton mode does not exit!\n";
		//_getchar_nolock();
		exit(1);
	}
}


/// ********************************************************** GMM Clustering ************************************************************///
// 
// GMM clustering
// 
// Parameter in
//  nInitMode: the initialization method of the GMM
//  bSucFlag:  the success flag of the GMM, if the GMM returns the number of clusters smaller than the specified one, then it fails
//  nMode: fast mode or not (no effect when m_nFeatureCnt > 1)
//
// Return
//  the clustering result map
// 
TGrayImage<int>* TPixelBasedClustering::GaussMixClustering(int nInitMode, bool& bSucFlag, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	// univariate - possibility of fast mode.
	if (m_nFeatureCnt == 1)
	{
		if (!(nMode & FAST))
			GaussMixParamEst( nInitMode, dbRatio, nInit1, nInit2, nIterFinal );
		else
		{
			// For fast mode, the parameter is based on histogram. This significantly
			// improves the speed. However, it assumes that the features are univariate
			// and are of range [0 255], and thus are only useful for intensity features.
			int *pHist = m_ppSrcImg[0]->GetIntensityHisto( -0.4, 255.6, 256, m_pMask );  // histogram (interesting respresentation)
			GaussMixParamEst( nInitMode, pHist, 256, dbRatio, nInit1, nInit2, nIterFinal );
			delete[] pHist;
		}

		return GaussMixClassify( bSucFlag );
	}
	else if (m_nFeatureCnt > 1) // multivariate.
	{
		GaussMixParamEst( nInitMode, dbRatio, nInit1, nInit2, nIterFinal );
		return GaussMixClassify( bSucFlag );
	}
	else
		return NULL;

}

//
///// ********************************************************** GMM Clustering ************************************************************///
//// 
//// GMM clustering in univariate feature case.
//// 
//// Parameter in
////  nInitMode: the initialization method of the GMM
////  bSucFlag:  the success flag of the GMM, if the GMM returns the number of clusters smaller than the specified one, then it fails
////  pFeature:  the univariate input feature image
////  pMask: the mask image
////  nClusterCnt: number of classes
////  nMode: fast mode or not
////
//// Return
////  the clustering result map
//// 
//TGrayImage<int>* TPixelBasedClustering::GaussMixClustering( int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>* pFeature, TMonoImage* pMask, int nClusterCnt, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal )
//{
//	if (!(nMode & FAST))
//		GaussMixParamEst( nInitMode, &pFeature, 1, pMask, nClusterCnt, dbRatio, nInit1, nInit2, nIterFinal );
//	else
//	{
//		// For fast mode, the parameter is based on histogram. This significantly
//		// improves the speed. However, it assumes that the features are univariate
//		// and are of range [0 255], and thus are only useful for intensity features.
//		int *pHist = pFeature->GetIntensityHisto( -0.4, 255.6, 256, pMask );  // histogram (interesting respresentation)
//		GaussMixParamEst( nInitMode, pHist, 256, nClusterCnt, dbRatio, nInit1, nInit2, nIterFinal );
//		delete[] pHist;
//	}
//
//	return GaussMixClassify( bSucFlag, &pFeature, pMask );
//}
//
//// 
//// GMM clustering in multivariate feature case.
//// 
//// Parameter in
////  nInitMode: the initialization method of the GMM
////  bSucFlag:  the success flag of the GMM, if the GMM returns the number of clusters smaller than the specified one, then it fails
////  ppFeatures: the multivariate input feature image
////  nFeatureCnt: feature space dimensionality
////  pMask: the mask image
////  nClusterCnt: number of classes
////
//// Return
////  the clustering result map
//// 
//TGrayImage<int>* TPixelBasedClustering::GaussMixClustering( int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>** ppFeatures, int nFeatureCnt, TMonoImage* pMask, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal )
//{
//	GaussMixParamEst( nInitMode, ppFeatures, nFeatureCnt, pMask, nClusterCnt, dbRatio, nInit1, nInit2, nIterFinal );
//
//	return GaussMixClassify( bSucFlag, ppFeatures, pMask );
//}

// 
// Parameter estimation for GMM using expectation maximization based on histogram
// 
void TPixelBasedClustering::GaussMixParamEst(int nInitMode, int* pHist, int nHistBinCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	// Allocate memory
	// Alloc(m_nClusterCnt, 1);
	TVector **ppMuBest = AllocVectors(m_nClusterCnt, 1);  // for preserving the best
	TMatrix **ppCovInvBest = AllocMatrice(m_nClusterCnt, 1);  // for preserving the best
	double* dbPriBest = new double[m_nClusterCnt];  // for preserving the best
	double* dbDetermBest = new double[m_nClusterCnt];  // for preserving the best

	// Find the range and variance of each feature
	double* VecMin = new double[1];
	double* VecMax = new double[1];
	double* Variance = new double[1];
	VecMin[0] = 0;
	VecMax[0] = nHistBinCnt - 1;
	
	int i, nCnt = 0;
	double dbMu = 0;
	Variance[0] = 0;
	for ( i = 0; i < nHistBinCnt - 1; i++ )
	{
		dbMu += ( pHist[i] * i );
		nCnt += pHist[i];
	}
	dbMu /= nCnt;

	for ( i = 0; i < nHistBinCnt - 1; i++ )
		Variance[0] += ( pHist[i] * ( i - dbMu ) * ( i - dbMu ) );

	Variance[0] /= nCnt;

	// Estimation
	double dbEnergy, dbMinEnergy = -999999999999;
	int nIter = 0;  // iteration number
	do {
		PixelInitParam(nInitMode, pHist, nHistBinCnt, Variance);  // randomly initialize the parameters
		dbEnergy = GaussMixEMIterations( dbRatio, nInit1, pHist, nHistBinCnt ); // 10 times?

		if (dbMinEnergy < 0 ||  // first iteration
		    dbEnergy < dbMinEnergy)  
			{
				dbMinEnergy = dbEnergy;
				for ( i = 0; i < m_nClusterCnt; i++ )
				{
					ppMuBest[i]->CopyFrom(m_ppMu[i]);
					ppCovInvBest[i]->CopyFrom(m_ppCovInv[i]);
				}
				memcpy(dbPriBest, m_dbPri, m_nClusterCnt * sizeof(double));
				memcpy(dbDetermBest, m_dbDeterm, m_nClusterCnt * sizeof(double));
			}
	} while ( ++nIter < nInit2 );

	// Set the parameters with the obtained best and continue with the EM
	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		m_ppMu[i]->CopyFrom(ppMuBest[i]);
		m_ppCovInv[i]->CopyFrom(ppCovInvBest[i]);
	}
	memcpy(m_dbPri, dbPriBest, m_nClusterCnt * sizeof(double));
	memcpy(m_dbDeterm, dbDetermBest, m_nClusterCnt * sizeof(double));

	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);
	
	// Further iterations to converge
	dbMinEnergy = GaussMixEMIterations( dbRatio, nIterFinal, pHist, nHistBinCnt );

	//printf("The final minimum energy is %.6f\n", dbMinEnergy);

	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	//SortCluster();

	for (i = 0; i < m_nClusterCnt; i ++)
	{
		delete ppMuBest[i];
		delete ppCovInvBest[i];
	}
	delete[] ppMuBest;
	delete[] ppCovInvBest;
	delete[] dbPriBest;
	delete[] dbDetermBest;
	delete[] VecMin;
	delete[] VecMax;
	delete[] Variance;
}

// 
// Parameter estimation for GMM using expectation maximization based on multiple features
// 
void TPixelBasedClustering::GaussMixParamEst(int nInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	// Allocate memory
	// Alloc(m_nClusterCnt, m_nFeatureCnt);
	TVector **ppMuBest = AllocVectors(m_nClusterCnt, m_nFeatureCnt);  // for preserving the best
	TMatrix **ppCovInvBest = AllocMatrice(m_nClusterCnt, m_nFeatureCnt);  // for preserving the best
	double* dbPriBest = new double[m_nClusterCnt];  // for preserving the best
	double* dbDetermBest = new double[m_nClusterCnt];  // for preserving the best

	// Find the range and variance of each feature
	double* VecMin = new double[m_nFeatureCnt];
	double* VecMax = new double[m_nFeatureCnt];
	double* Mu = new double[m_nFeatureCnt];
	double* Variance = new double[m_nFeatureCnt];
	int i;
	for ( i = 0; i < m_nFeatureCnt; i++ )
	{
		((TImage*) m_ppSrcImg[i])->GetIntensityRange(VecMax[i], VecMin[i], m_pMask);
		Mu[i] = m_ppSrcImg[i]->GetIntensityMean(m_pMask);
	}

	int nCnt = 0;
	for ( i = 0; i < m_nFeatureCnt; i++ )
	{
		Variance[i] = 0;
		nCnt = 0;
		for ( int j = 0; j < m_ppSrcImg[i]->Height(); j++ )
			for ( int k = 0; k < m_ppSrcImg[i]->Width(); k++ )
			{
				if ( m_pMask && ! (*m_pMask)(j, k) )
					continue;
				Variance[i] += ( ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) * ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) );
				nCnt++;
			}
		Variance[i] /= nCnt;
	}

	// Parameter estimation
	double dbEnergy, dbMinEnergy = -999999999999;
	int nIter = 0;  // iteration number
	do {
		PixelInitParam(nInitMode, VecMin, VecMax, Variance);  // initialize randomly the parameters
		dbEnergy = GaussMixEMIterations( dbRatio, nInit1); //10 times?

		if (dbMinEnergy < 0 ||  // first iteration
		    dbEnergy < dbMinEnergy)  
		{
			dbMinEnergy = dbEnergy;
			for (i = 0; i < m_nClusterCnt; i ++)
			{
				ppMuBest[i]->CopyFrom(m_ppMu[i]);
				ppCovInvBest[i]->CopyFrom(m_ppCovInv[i]);
			}
			memcpy(dbPriBest, m_dbPri, m_nClusterCnt * sizeof(double));
			memcpy(dbDetermBest, m_dbDeterm, m_nClusterCnt * sizeof(double));
		}
		//printf("%d", nIter);
	} while ( ++nIter < nInit2 );

	// Set the parameters with the obtained best and continue with the EM
	for (i = 0; i < m_nClusterCnt; i ++)
	{
		m_ppMu[i]->CopyFrom(ppMuBest[i]);
		m_ppCovInv[i]->CopyFrom(ppCovInvBest[i]);
	}
	memcpy(m_dbPri, dbPriBest, m_nClusterCnt * sizeof(double));
	memcpy(m_dbDeterm, dbDetermBest, m_nClusterCnt * sizeof(double));

	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);
	
	// Further iterations to converge
	dbMinEnergy = GaussMixEMIterations( dbRatio, nIterFinal );

	//printf("The final minimum energy is %.6f\n", dbMinEnergy);

	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	//SortCluster();

	//// Show the parameters
	//ShowParam();

	for (i = 0; i < m_nClusterCnt; i ++)
	{
		delete ppMuBest[i];
		delete ppCovInvBest[i];
	}
	delete[] ppMuBest;
	delete[] ppCovInvBest;
	delete[] dbPriBest;
	delete[] dbDetermBest;
	delete[] VecMin;
	delete[] VecMax;
	delete[] Mu;
	delete[] Variance;
}


// 
// The EM for GMM in fast mode (i.e. using histogram)
// 
double TPixelBasedClustering::GaussMixEMIterations( double dbRatio, int nIterations, int* pHist, int nHistBinCnt )
{
	TVector **ppMuComp = AllocVectors(m_nClusterCnt, 1);
	TMatrix **ppCovComp = AllocMatrice(m_nClusterCnt, 1);
	double* dbPri = new double[m_nClusterCnt];  // for the current compuation
	int nIter = 0, i;
	double dbOldEnergy, dbNewEnergy = -999999999999, dbPriTotal;
	TVector VecFeature(1, false);  // feature vector
	TVector VecProb(m_nClusterCnt, false);  // probability
	do {
		dbOldEnergy = dbNewEnergy;
		dbNewEnergy = 0;
		SetZeroVector(ppMuComp, m_nClusterCnt);
		SetZeroMatrix(ppCovComp, m_nClusterCnt);
		memset(dbPri, 0, sizeof(double) * m_nClusterCnt);
		for ( i = 0; i < nHistBinCnt; i++ )
		{
			// E step
			VecFeature(0) = i;
			dbNewEnergy += ( pHist[i] * GaussMixExpectation(VecProb, VecFeature) );

			// M step
			VecProb *= pHist[i];
			GaussMixMaximization(ppMuComp, ppCovComp, dbPri, VecProb, VecFeature);
		}

		if ( nIter == 0 || ( ( nIter < nIterations ) &&( fabs(dbNewEnergy - dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) )
		{
			// Compute the total of the priors
			dbPriTotal = 0;
			for ( i = 0; i < m_nClusterCnt; i++ )
			{
				if ( dbPri[i] <= ZERO_THRESH )
					dbPri[i] += ZERO_THRESH;
				dbPriTotal += dbPri[i];
			}

			// Update the parameters
			for ( i = 0; i < m_nClusterCnt; i++ )
			{
				(*m_ppMu[i])(0) = (*ppMuComp[i])(0) / dbPri[i];  // mean
				(*ppCovComp[i])(0, 0) = (*ppCovComp[i])(0, 0) / dbPri[i] - (*m_ppMu[i])(0) * (*m_ppMu[i])(0);  // covariance
				if ( (*ppCovComp[i])(0, 0) <= ZERO_THRESH )
					(*ppCovComp[i])(0, 0) += REG_COEF;
				(*m_ppCovInv[i])(0, 0) = 1 / (*ppCovComp[i])(0, 0);  // inverse of covariance
				m_dbPri[i] = dbPri[i] / dbPriTotal;  // prior
				m_dbDeterm[i] = sqrt((*ppCovComp[i])(0, 0));  // determinant
			}
		}

		nIter++;

	} while ( ( nIter == 1 ) || ( ( nIter <= nIterations ) && ( fabs(dbNewEnergy - dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) );


	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		delete ppMuComp[i];
		delete ppCovComp[i];
	}
	delete[] ppMuComp;
	delete[] ppCovComp;
	delete[] dbPri;

	return dbNewEnergy;
}

// 
// The EM for GMM
// 
double TPixelBasedClustering::GaussMixEMIterations(double dbRatio, int nIterations)
{
	TVector **ppMuComp = AllocVectors(m_nClusterCnt, m_nFeatureCnt);
	TMatrix **ppCovComp = AllocMatrice(m_nClusterCnt, m_nFeatureCnt);
	double* dbPri = new double[m_nClusterCnt];  // for the current compuation
	int nIter = 0, i, k;
	double dbOldEnergy, dbNewEnergy = -999999999999, dbPriTotal;
	TVector VecFeature(m_nFeatureCnt, false);  // feature vector
	TVector VecProb(m_nClusterCnt, false);  // probability
	do {
		dbOldEnergy = dbNewEnergy;
		dbNewEnergy = 0;
		SetZeroVector(ppMuComp, m_nClusterCnt);
		SetZeroMatrix(ppCovComp, m_nClusterCnt);
		memset(dbPri, 0, sizeof(double) * m_nClusterCnt);
		for ( int r = m_ppSrcImg[0]->Height() - 1; r >= 0; r-- )
			for ( int c = m_ppSrcImg[0]->Width() - 1; c >= 0; c-- )
			{
				if (m_pMask && !(*m_pMask)(r, c))
					continue;

				// E step
				for (i = 0; i < m_nFeatureCnt; i ++)
					VecFeature(i) = (*m_ppSrcImg[i])(r, c);
				dbNewEnergy += GaussMixExpectation(VecProb, VecFeature);

				// M step
				GaussMixMaximization(ppMuComp, ppCovComp, dbPri, VecProb, VecFeature);
			}

		if ( nIter == 0 || ( ( nIter < nIterations ) &&( fabs(dbNewEnergy - dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) )
		{
			// Compute the total of the priors
			dbPriTotal = 0;
			for ( k = 0; k < m_nClusterCnt; k++ )
			{
				if ( dbPri[k] <= ZERO_THRESH )
					dbPri[k] += ZERO_THRESH;
				dbPriTotal += dbPri[k];
			}

			// Update the parameters
			for ( k = 0; k < m_nClusterCnt; k++ )
			{
				for ( i = 0; i < m_nFeatureCnt; i++ )  // mean
					(*m_ppMu[k])(i) = (*ppMuComp[k])(i) / dbPri[k];

				for ( i = 0; i < m_nFeatureCnt; i++ )  // covariance
				{
					(*ppCovComp[k])(i, i) = (*ppCovComp[k])(i, i) / dbPri[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(i);
					for ( int j = i + 1; j < m_nFeatureCnt; j++ )
					{
						(*ppCovComp[k])(i, j) = (*ppCovComp[k])(i, j) / dbPri[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
						(*ppCovComp[k])(j, i) = (*ppCovComp[k])(i, j);
					}
				}
					
				if ( ppCovComp[k]->Eigen() <= ZERO_THRESH )
					for ( i = 0; i < m_nFeatureCnt; i++ )  // covariance
						(*ppCovComp[k])(i, i) += REG_COEF;

				ppCovComp[k]->Inv(m_ppCovInv[k]);  // inverse of covariance
				m_dbPri[k] = dbPri[k] / dbPriTotal;  // prior
				m_dbDeterm[k] = sqrt(ppCovComp[k]->Determinant());  // determinant
			}
		}

		nIter++;

	} while ( ( nIter == 1 ) || ( ( nIter <= nIterations ) && ( fabs(dbNewEnergy - dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) );

	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		delete ppMuComp[i];
		delete ppCovComp[i];
	}
	delete[] ppMuComp;
	delete[] ppCovComp;
	delete[] dbPri;
	return dbNewEnergy;
}

// 
// The expectation step for GMM
// 
// Param in
//  VecFeature: the feature
//
// Param out
//  VecProb: the probability
// 
// Return
//  the associated energy
// 
double TPixelBasedClustering::GaussMixExpectation(TVector& VecProb, TVector& VecFeature)
{
 //   TVector* m_pVecTemp1 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation
	//TVector* m_pVecTemp2 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation
	//double dbTotal = 0;
	//int i;
	//for ( i = 0; i < m_nClusterCnt; i++ )
	//{
	//	for ( int j = 0; j < m_nFeatureCnt; j++ )
	//		(*m_pVecTemp1)(j) = VecFeature(j) - (*m_ppMu[i])(j);
	//	m_pVecTemp1->Mul(m_pVecTemp2, m_ppCovInv[i]);
	//	VecProb(i) = m_dbPri[i] * exp( -m_pVecTemp2->VecMul(m_pVecTemp1) / 2) / m_dbDeterm[i];
	//	dbTotal += VecProb(i);
	//}

	//delete m_pVecTemp1;
	//delete m_pVecTemp2;

	//for ( i = 0; i < m_nClusterCnt; i++ )
	//	VecProb(i) = VecProb(i) / dbTotal;

	//return -log(dbTotal);

    TVector* m_pVecTemp1 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation
	TVector* m_pVecTemp2 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation
	double dbTotal;

	double* U = new double[m_nClusterCnt];

	double* dbLnPrior = new double[m_nClusterCnt];
	double* dbLnSigma = new double[m_nClusterCnt];
	int i;
	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		dbLnPrior[i] = log(m_dbPri[i]);
		dbLnSigma[i] = log(m_dbDeterm[i]);
	}

	double dbMin = 999999999999;
	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		for ( int j = 0; j < m_nFeatureCnt; j++ )
			(*m_pVecTemp1)(j) = VecFeature(j) - (*m_ppMu[i])(j);
		m_pVecTemp1->Mul(m_pVecTemp2, m_ppCovInv[i]);
		U[i] = m_pVecTemp2->VecMul(m_pVecTemp1) / 2 + dbLnSigma[i];
		if ( U[i] < dbMin ) 
			dbMin = U[i];
	}

	delete m_pVecTemp1;
	delete m_pVecTemp2;
    
	dbTotal = 0;
	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		U[i] -= dbMin;
		VecProb(i) = exp(-U[i]) * m_dbPri[i];
		dbTotal += VecProb(i);
	}

	for ( i = 0; i < m_nClusterCnt; i++ )
		VecProb(i) = VecProb(i) / dbTotal;

	//dbTotal = 0;
	//for ( i = 0; i < m_nClusterCnt; i++ )
	//	dbTotal -= ( -U[i] - dbMin + dbLnPrior[i] ) * VecProb(i);

	dbTotal = -log(dbTotal) + dbMin;

	delete[] U;
	delete[] dbLnPrior;
	delete[] dbLnSigma;

	return dbTotal;
}

// 
// The maximization step for GMM
// 
// Param in
//  VecProb: the probability
//  VecFeature: the feature
//
// Param out
//  ppMu, ppCov, dbPri: the new sum of mean, covariance and prior
// 
void TPixelBasedClustering::GaussMixMaximization(TVector** ppMu, TMatrix** ppCov, double *dbPri, TVector& VecProb, TVector& VecFeature)
{
	double dbTmp;
	for ( int k = 0; k < m_nClusterCnt; k++ )
	{
		for ( int i = 0; i < m_nFeatureCnt; i++ )
		{
			dbTmp = VecFeature(i) * VecProb(k);
			(*ppMu[k])(i) += dbTmp;
			for ( int j = i; j < m_nFeatureCnt; j++ )
				(*ppCov[k])(i, j) += ( dbTmp * VecFeature(j) );
		}
		dbPri[k] += VecProb(k);
	}
}

// 
// Classification based on the obtained Guass mixture parameters
// 
// Parameter in
//  ppFeatures: the feature vector image
//  pMask: the mask image
//
// Return
//  the classification map
// 
TGrayImage<int>* TPixelBasedClustering::GaussMixClassify( bool& bSucFlag )
{
	////int nWidth = ppFeatures[0]->Width();
	////int nHeight = ppFeatures[0]->Height();
	////TGrayImage<int>* pMap = new TGrayImage<int>(nWidth, nHeight);

	////TVector* m_pVecTemp1 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation
	////TVector* m_pVecTemp2 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation

	////// Classification
	////TVector VecFeature(m_nFeatureCnt, false);  // feature vector
	////TVector VecProb(m_nClusterCnt, false);  // probability
	////double dbMax;
	////int nLabel, i;
	////int* nLabelFlag = new int[m_nClusterCnt];
	////memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));
	////for ( int r = 0; r < nHeight; r++ )
	////{
	////	for ( int c = 0; c < nWidth; c++ )
	////	{
	////		if (pMask && !(*pMask)(r, c))
	////		{
	////			(*pMap)(r, c) = -2;
	////			continue;
	////		}
	////
	////		for ( i = 0; i < m_nFeatureCnt; i++ )
	////			VecFeature(i) = (*ppFeatures[i])(r, c);
	////		for ( i = 0; i < m_nClusterCnt; i++ )
	////		{
	////			for ( int j = 0; j < m_nFeatureCnt; j++ )
	////				(*m_pVecTemp1)(j) = VecFeature(j) - (*m_ppMu[i])(j);
	////			m_pVecTemp1->Mul(m_pVecTemp2, m_ppCovInv[i]);
	////			VecProb(i) = m_dbPri[i] * exp( -m_pVecTemp2->VecMul(m_pVecTemp1) / 2 ) / m_dbDeterm[i];
	////		}

	////		dbMax = -999999999999;
	////		for ( i = 0; i < m_nClusterCnt; i++ )
	////		{
	////			if (VecProb(i) > dbMax)
	////			{
	////				dbMax = VecProb(i);
	////				nLabel = i;
	////			}
	////		}
	////		(*pMap)(r, c) = nLabel;
	////		if ( nLabelFlag[nLabel] == 0 )
	////			nLabelFlag[nLabel] = 1;
	////	}
	////}
	////
	////nLabel = 0;
	////for ( i = 0; i < m_nClusterCnt; i++ )
	////	nLabel += nLabelFlag[i];
	////if ( nLabel == m_nClusterCnt )
	////	bSucFlag = true;
	////else
	////	bSucFlag = false;

	////delete[] nLabelFlag;


	int nWidth = m_ppSrcImg[0]->Width();
	int nHeight = m_ppSrcImg[0]->Height();
	TGrayImage<int>* pMap = new TGrayImage<int>(nWidth, nHeight);

	TVector* m_pVecTemp1 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation
	TVector* m_pVecTemp2 = new TVector(m_nFeatureCnt, false);  // a temperoray vector during the computation

	// Classification
	double dbMax, dbTmp;
	int nLabel;
	int* nLabelFlag = new int[m_nClusterCnt];
	memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));
    
	double* U = new double[m_nClusterCnt];

	double* dbLnSigma = new double[m_nClusterCnt];
	int i;
	for ( i = 0; i < m_nClusterCnt; i++ )
		dbLnSigma[i] = log(m_dbDeterm[i]);

	dbMax = 999999999999;
	for ( int r = 0; r < nHeight; r++ )
		for ( int c = 0; c < nWidth; c++ )
		{
			if (m_pMask && !(*m_pMask)(r, c))
			{
				(*pMap)(r, c) = -2;
				continue;
			}

			for ( i = 0; i < m_nClusterCnt; i++ )
			{
				for ( int j = 0; j < m_nFeatureCnt; j++ )
					(*m_pVecTemp1)(j) = (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j);
				m_pVecTemp1->Mul(m_pVecTemp2, m_ppCovInv[i]);
				U[i] = m_pVecTemp2->VecMul(m_pVecTemp1) / 2 + dbLnSigma[i];
				if ( U[i] < dbMax )
					dbMax = U[i];
			}

			for ( i = 0; i < m_nClusterCnt; i++ )
				U[i] -= dbMax;

			dbMax = -999999999999;
			for ( i = 0; i < m_nClusterCnt; i++ )
			{
				dbTmp = exp(-U[i]) * m_dbPri[i];
				if ( dbTmp > dbMax )
				{
					dbMax = dbTmp;
					nLabel = i;
				}
			}

			(*pMap)(r, c) = nLabel;
			if ( nLabelFlag[nLabel] == 0 )
				nLabelFlag[nLabel] = 1;
		}
	
	nLabel = 0;
	for ( i = 0; i < m_nClusterCnt; i++ )
		nLabel += nLabelFlag[i];
	if ( nLabel == m_nClusterCnt )
		bSucFlag = true;
	else
		bSucFlag = false;

	delete[] nLabelFlag;
	delete[] U;
	delete[] dbLnSigma;


 //   TGrayImage<int>* pDispMap = new TGrayImage<int>(nWidth, nHeight);
	//for ( int r = 0; r < nHeight; r++ )
	//	for ( int c = 0; c < nWidth; c++ )
	//	{
	//		if ( (nLabel = (*pMap)(r, c)) >= 0 )
	//		{
	//			if (nLabel == 0)              //QIN_DEBUG: Here I only consider at most three classes
	//				(*pDispMap)(r, c) = 128;
	//			if (nLabel == 1)
	//				(*pDispMap)(r, c) = 0;
	//			if (nLabel == 2)
	//				(*pDispMap)(r, c) = 255;
	//			if (nLabel == 3)
	//				(*pDispMap)(r, c) = 64;
	//			if (nLabel == 4)
	//				(*pDispMap)(r, c) = 192;
	//			if (nLabel == 5)
	//				(*pDispMap)(r, c) = 32;
	//			if (nLabel == 6)
	//				(*pDispMap)(r, c) = 224;
	//		}
	//		else  // boundary label
	//			(*pDispMap)(r, c) = nLabel;
	//	}

	//pDispMap->Save("GMM_Pixel_Result.bmp");
	//delete pDispMap;

	delete m_pVecTemp1;
	delete m_pVecTemp2;

	return pMap;
}

/// ************************************************************************************************************************************* ///
/// ********************************************************** KM Clustering ************************************************************///
// 
// KM clustering in univariate feature case.
// 
// Parameter in
//  nMode: mode - no effect when multivariate
//
// Return
//  the clustering result
// 
TGrayImage<int>* TPixelBasedClustering::KMClustering(int nInitMode, bool& bSucFlag, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	if (m_nFeatureCnt == 1)
	{
		if (!(nMode & FAST))
			KMParamEst(nInitMode, dbRatio, nInit1, nInit2, nIterFinal);
		else
		{
			// For fast mode, the parameter is based on histogram. This significantly
			// improves the speed. However, it assumes that the features are univariate
			// and is of range [0 255], and thus is only useful for intensity feature.
			int *pHist = m_ppSrcImg[0]->GetIntensityHisto( -0.4, 255.6, 256, m_pMask );  // histogram (interesting respresentation)
			KMParamEst(nInitMode, pHist, 256, dbRatio, nInit1, nInit2, nIterFinal);
			delete[] pHist;
		}

		//// Show the parameters
		//ShowParam();

		return KMClassify(bSucFlag);
	}
	else if (m_nFeatureCnt > 1)
	{
		KMParamEst(nInitMode, dbRatio, nInit1, nInit2, nIterFinal);
		return KMClassify(bSucFlag);
	}
	else 
		return NULL;
}

/// ************************************************************************************************************************************* ///
/// ********************************************************** OTSU Thresholding ************************************************************///
// 
// OTSU thresholding in univariate feature case.
// 
// Parameter in
//  nMode: mode - no effect when multivariate
//
// Return
//  the clustering result
// 

// Calculate histogram
//Shuhrat otsu
TGrayImage<int>* TPixelBasedClustering::OTSUClustering(bool& bSucFlag, int& o_threshold){

if (m_nFeatureCnt == 1)
	{
	//double dbMax, dbMin;
		//	((TImage*) m_ppSrcImg[0])->GetIntensityRange(dbMax, dbMin, 0);			
			int *histData = m_ppSrcImg[0]->GetRawIntensityHisto(256, m_pMask );  // histogram (interesting respresentation)


	// Total number of pixels
	int total=0;
	for (int i=0;i<256;i++)
		total+=histData[i];

	float sum = 0;
	for (int t=0 ; t<256 ; t++) sum += t * histData[t];

	float sumB = 0;
	int wB = 0;
	int wF = 0;

	float varMax = 0;
	int threshold = 0;

	for (int t=0 ; t<256 ; t++) {
	   wB += histData[t];               // Weight Background
	   if (wB == 0) continue;

	   wF = total - wB;                 // Weight Foreground
	   if (wF == 0) break;

	   sumB += (float) (t * histData[t]);

	   float mB = sumB / wB;            // Mean Background
	   float mF = (sum - sumB) / wF;    // Mean Foreground

	   // Calculate Between Class Variance
	   float varBetween = (float)wB * (float)wF * (mB - mF) * (mB - mF);

	   // Check if new maximum found
	   if (varBetween > varMax) {
		  varMax = varBetween;
		  threshold = t;
   }
}
	cout<<"The threshold in OTSUClustering is"<<threshold;
	o_threshold=threshold;
        delete[] histData; 
		return OTSUClassify(bSucFlag,threshold);
	}
	else 
		return NULL;
}

// Calculate histogram
TGrayImage<int>* TPixelBasedClustering::OTSUClustering(bool& bSucFlag, int& threshold, bool bFast){

if (m_nFeatureCnt == 1)
		return OTSUClassify(bSucFlag,threshold);	
	else 
		return NULL;
}


// 
// Parameter estimation for KM based on histogram
// 
void TPixelBasedClustering::KMParamEst(int nInitMode, int* pHist, int nHistBinCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	// Allocate memory
	// Alloc(m_nClusterCnt, 1);
	TVector **ppMuBest = AllocVectors(m_nClusterCnt, 1);  // for preserving the best

	// Find the range and variance of each feature
	double* VecMin = new double[1];
	double* VecMax = new double[1];
	double* Variance = new double[1];
	VecMin[0] = 0;
	VecMax[0] = nHistBinCnt - 1;
	
	int i, nCnt = 0;
	double dbMu = 0;
	Variance[0] = 0;
	for ( i = 0; i < nHistBinCnt - 1; i++ )
	{
		dbMu += ( pHist[i] * i );
		nCnt += pHist[i];
	}
	dbMu /= nCnt;

	for ( i = 0; i < nHistBinCnt - 1; i++ )
		Variance[0] += ( pHist[i] * ( i - dbMu ) * ( i - dbMu ) );

	Variance[0] /= nCnt;


	// Estimation
	double dbEnergy, dbMinEnergy = -999999999999;
	int j, nIter = 0;  // iteration number
	do {
		PixelInitParam(nInitMode, pHist, nHistBinCnt, Variance);
		dbEnergy = KMIterations( dbRatio, nInit1, pHist, nHistBinCnt ); // 10 times?

		if (dbMinEnergy < 0 ||  // first iteration
		    dbEnergy < dbMinEnergy)  
			{
				dbMinEnergy = dbEnergy;
				for (i = 0; i < m_nClusterCnt; i ++)
					ppMuBest[i]->CopyFrom(m_ppMu[i]);
			}
	} while ( ++nIter < nInit2 );

	// Set the parameters with the obtained best and continue
	for (i = 0; i < m_nClusterCnt; i ++)
		m_ppMu[i]->CopyFrom(ppMuBest[i]);

	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);
	
	// Further iterations to converge
	dbMinEnergy = KMIterations( dbRatio, nIterFinal, pHist, nHistBinCnt );

	//printf("The final minimum energy is %.6f\n", dbMinEnergy);

    // Update Other parameters needed
	int nTotalCnt = 0, nMinIdx;
	int* nSampleCnt = new int[m_nClusterCnt];
	memset(nSampleCnt, 0, sizeof(int) * m_nClusterCnt);
	double dbTmp, dbMin;
	for ( i = 0; i < nHistBinCnt; i++ )
	{
		dbMin = 999999999999;
		for ( j = 0; j < m_nClusterCnt; j++ )
		{
			dbTmp = ( i - (*m_ppMu[j])(0) )*( i - (*m_ppMu[j])(0) );
			if ( dbTmp < dbMin )
			{
				dbMin = dbTmp;
				nMinIdx = j;
			}
		}
		(*m_ppCovInv[nMinIdx])(0, 0) += ( i * i * pHist[i] );
		nSampleCnt[nMinIdx] += pHist[i];
	}

	for (i = 0; i < m_nClusterCnt; i ++)
	{
		if ( nSampleCnt[i] > 0 )
			(*m_ppCovInv[i])(0, 0) = (*m_ppCovInv[i])(0, 0) / nSampleCnt[i] - (*m_ppMu[i])(0) * (*m_ppMu[i])(0);
		else
			nSampleCnt[i]++;

		if ( (*m_ppCovInv[i])(0, 0) <= ZERO_THRESH )
			(*m_ppCovInv[i])(0, 0) += Variance[0];

		m_dbDeterm[i] = sqrt((*m_ppCovInv[i])(0, 0));

		(*m_ppCovInv[i])(0, 0) = 1 / (*m_ppCovInv[i])(0, 0);

		nTotalCnt += nSampleCnt[i];
	}

	for (i = 0; i < m_nClusterCnt; i ++)
		m_dbPri[i] = ((double) nSampleCnt[i]) / ((double) nTotalCnt);  // prior

	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	//SortCluster();

	for (i = 0; i < m_nClusterCnt; i ++)
		delete ppMuBest[i];
	delete[] ppMuBest;
	delete[] nSampleCnt;
	delete[] VecMin;
	delete[] VecMax;
	delete[] Variance;
}

// 
// Parameter estimation for KM based on multiple features
// 
void TPixelBasedClustering::KMParamEst(int nInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	// Allocate memory
	// Alloc(m_nClusterCnt, m_nFeatureCnt);
	TVector **ppMuBest = AllocVectors(m_nClusterCnt, m_nFeatureCnt);  // for preserving the best

	// Find the range and variance of each feature
	double* VecMin = new double[m_nFeatureCnt];
	double* VecMax = new double[m_nFeatureCnt];
	double* Mu = new double[m_nFeatureCnt];
	double* Variance = new double[m_nFeatureCnt];
	int i, j, k;
	for ( i = 0; i < m_nFeatureCnt; i++ )
	{
		((TImage*) m_ppSrcImg[i])->GetIntensityRange(VecMax[i], VecMin[i], m_pMask);
		Mu[i] = m_ppSrcImg[i]->GetIntensityMean(m_pMask);
	}

	int nCnt = 0;
	double dbAvgVar = 0.0;
	for ( i = 0; i < m_nFeatureCnt; i++ )
	{
		Variance[i] = 0;
		nCnt = 0;
		for ( j = 0; j < m_ppSrcImg[i]->Height(); j++ )
			for ( k = 0; k < m_ppSrcImg[i]->Width(); k++ )
			{
				if ( m_pMask && ! (*m_pMask)(j, k) )
					continue;
				Variance[i] += ( ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) * ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) );
				nCnt++;
			}
		Variance[i] /= nCnt;
		dbAvgVar += Variance[i];
	}
	dbAvgVar /= m_nFeatureCnt;

	// Parameter estimation
	double dbEnergy, dbMinEnergy = -999999999999;
	int nIter = 0;  // iteration number and sample count
	do {
		PixelInitParam(nInitMode, VecMin, VecMax, Variance);
		dbEnergy = KMIterations( dbRatio, nInit1 );  //10 times?

		if (dbMinEnergy < 0 ||  // first iteration
		    dbEnergy < dbMinEnergy)  
		{
			dbMinEnergy = dbEnergy;
			for (i = 0; i < m_nClusterCnt; i ++)
				ppMuBest[i]->CopyFrom(m_ppMu[i]);

		}
		//printf("%d", nIter);
	} while ( ++nIter < nInit2 );

	// Set the parameters with the obtained best and continue
	for ( i = 0; i < m_nClusterCnt; i++ )
		m_ppMu[i]->CopyFrom(ppMuBest[i]);

	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);

	// Further iterations to converge
	dbMinEnergy = KMIterations( dbRatio, nIterFinal );

	//printf("The final minimum energy is %.6f\n", dbMinEnergy);

    // Update other parameters needed
	int nTotalCnt = 0, nMinIdx;
	int* nSampleCnt = new int[m_nClusterCnt];
	memset(nSampleCnt, 0, sizeof(int) * m_nClusterCnt);
	double dbTmp, dbMin;
	for (int r = m_ppSrcImg[0]->Height() - 1; r >= 0; r --)
		for (int c = m_ppSrcImg[0]->Width() - 1; c >= 0; c --)
		{
			if (m_pMask && !(*m_pMask)(r, c))
				continue;

			dbMin = 999999999999;
			for ( i = 0; i < m_nClusterCnt; i++ )
			{
				dbTmp = 0;
				for ( j = 0; j < m_nFeatureCnt; j++ )
					dbTmp += ( ( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) )*( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) ) );
				if ( dbTmp < dbMin )
				{
					dbMin = dbTmp;
					nMinIdx = i;
				}
			}
			for ( i = 0; i < m_nFeatureCnt; i++ )
				for ( j = i ; j < m_nFeatureCnt; j++ )
				    (*m_ppCovInv[nMinIdx])(i, j) += ( (*m_ppSrcImg[i])(r, c) * (*m_ppSrcImg[j])(r, c) );
			nSampleCnt[nMinIdx]++;
		}


	for ( k = 0; k < m_nClusterCnt; k++ )
	{
		if ( nSampleCnt[k] > 0 )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ )  // covariance
			{
				(*m_ppCovInv[k])(i, i) = (*m_ppCovInv[k])(i, i) / nSampleCnt[k] -  (*m_ppMu[k])(i) * (*m_ppMu[k])(i);
				for ( j = i + 1; j < m_nFeatureCnt; j++ )
				{
					(*m_ppCovInv[k])(i, j) = (*m_ppCovInv[k])(i, j) / nSampleCnt[k] -  (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
					(*m_ppCovInv[k])(j, i) = (*m_ppCovInv[k])(i, j);
				}
			}
		}
		else
			nSampleCnt[k]++;

		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += dbAvgVar;
		}

		m_dbDeterm[k] = sqrt(m_ppCovInv[k]->Determinant());  // determinant

		m_ppCovInv[k]->Inv(m_ppCovInv[k]);  // inverse of covariance

		nTotalCnt += nSampleCnt[k];
	}

	for ( k = 0; k < m_nClusterCnt; k++ )
		m_dbPri[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);

	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	//SortCluster();

	//// Show the parameters
	//ShowParam();

	for (i = 0; i < m_nClusterCnt; i ++)
		delete ppMuBest[i];
	delete[] ppMuBest;
	delete[] nSampleCnt;
	delete[] VecMin;
	delete[] VecMax;
	delete[] Mu;
	delete[] Variance;
}

// 
// The univariate KM operation in fast mode (i.e. using histogram)
// 
double TPixelBasedClustering::KMIterations(double dbRatio, int nIterations, int* pHist, int nHistBinCnt)
{
	int* nSampleCnt = new int[m_nClusterCnt];  // for the current compuation
	int nIter = 0, i, j, nMinIdx;
	double dbOldEnergy, dbNewEnergy = -999999999999, dbMin, dbTmp;
	TVector **ppMu = AllocVectors(m_nClusterCnt, 1);  // for preserving the best

	do {
		dbOldEnergy = dbNewEnergy;
		dbNewEnergy = 0;
		memset(nSampleCnt, 0, sizeof(int) * m_nClusterCnt);
		SetZeroVector(ppMu , m_nClusterCnt);
		for ( i = 0; i < nHistBinCnt; i++ )
		{
			dbMin = 999999999999;
			for ( j = 0; j < m_nClusterCnt; j++ )
			{
				dbTmp = ( i - (*m_ppMu[j])(0) )*( i - (*m_ppMu[j])(0) );
				if ( dbTmp < dbMin )
				{
					dbMin = dbTmp;
					nMinIdx = j;
				}
			}
		    (*ppMu[nMinIdx])(0) += ( i * pHist[i] );
			nSampleCnt[nMinIdx] += pHist[i];
			dbNewEnergy += ( dbMin * pHist[i] );
		}

		if ( nIter == 0 || ( ( nIter < nIterations ) &&( fabs(dbNewEnergy -dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) )
		{
			// Update the parameters
			for (i = 0; i < m_nClusterCnt; i ++)
				if ( nSampleCnt[i] > 0 )
					(*m_ppMu[i])(0) = (*ppMu[i])(0)/nSampleCnt[i];
		}

		nIter++;

	} while ( ( nIter == 1 ) || ( ( nIter <= nIterations ) && ( fabs(dbNewEnergy -dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) );

	for (i = 0; i < m_nClusterCnt; i ++)
		delete ppMu[i];
	delete[] ppMu;
	delete[] nSampleCnt;

	return dbNewEnergy;
}

// 
// The multivariate KM operation
// 
double TPixelBasedClustering::KMIterations(double dbRatio, int nIterations)
{
	int* nSampleCnt = new int[m_nClusterCnt];  // for the current compuation
	int nIter = 0, i, j, nMinIdx;
	double dbOldEnergy, dbNewEnergy = -999999999999, dbMin, dbTmp;
	TVector **ppMu = AllocVectors(m_nClusterCnt, m_nFeatureCnt);  // for preserving the best
	
	do {
		dbOldEnergy = dbNewEnergy;
		dbNewEnergy = 0;
		memset(nSampleCnt, 0, sizeof(int) * m_nClusterCnt);
		SetZeroVector(ppMu , m_nClusterCnt);
		for (int r = m_ppSrcImg[0]->Height() - 1; r >= 0; r --)
			for (int c = m_ppSrcImg[0]->Width() - 1; c >= 0; c --)
			{
				if (m_pMask && !(*m_pMask)(r, c))
					continue;

				dbMin = 999999999999;
				for ( i = 0; i < m_nClusterCnt; i++ )
				{
					dbTmp = 0;
					for ( j = 0; j < m_nFeatureCnt; j++ )
						dbTmp += ( ( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) ) * ( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) ) );
					if ( dbTmp < dbMin )
					{
						dbMin = dbTmp;
						nMinIdx = i;
					}
				}
				for ( i = 0; i < m_nFeatureCnt; i++ )
					(*ppMu[nMinIdx])(i) += (*m_ppSrcImg[i])(r, c);
				nSampleCnt[nMinIdx]++;
				dbNewEnergy += dbMin;
			}

		if ( nIter == 0 || ( ( nIter < nIterations ) &&( fabs(dbNewEnergy -dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) )
		{
			// Update the parameters
			for (i = 0; i < m_nClusterCnt; i ++)
				for ( j = 0; j < m_nFeatureCnt; j++ )
				{
					if ( nSampleCnt[i] > 0 )
						(*m_ppMu[i])(j) = (*ppMu[i])(j)/nSampleCnt[i];
				}
		}

		nIter++;

	} while ( ( nIter == 1 ) || ( ( nIter <= nIterations ) && ( fabs(dbNewEnergy -dbOldEnergy) > dbRatio * fabs(dbOldEnergy) ) ) );

	for (i = 0; i < m_nClusterCnt; i ++)
		delete ppMu[i];
	delete[] ppMu;
	delete[] nSampleCnt;

	return dbNewEnergy;
}

// 
// Classification based on the obtained KM
//
// Return
//  the classification map
// 
TGrayImage<int>* TPixelBasedClustering::KMClassify(bool& bSucFlag)
{
	int nWidth = m_ppSrcImg[0]->Width();
	int nHeight = m_ppSrcImg[0]->Height();
	TGrayImage<int>* pMap = new TGrayImage<int>(nWidth, nHeight);

	// Classification
	double dbTmp, dbMin;
	int i, j, nLabel;
	int* nLabelFlag = new int[m_nClusterCnt];
	memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));

	for (int r = 0; r < nHeight; r ++)
		for (int c = 0; c < nWidth; c ++)
		{
			if (m_pMask && !(*m_pMask)(r, c))
			{
				(*pMap)(r, c) = -2;
				continue;
			}

			dbMin = 999999999999;
			for ( i = 0; i < m_nClusterCnt; i++ )
			{
				dbTmp = 0;
				for ( j = 0; j < m_nFeatureCnt; j++ )
					dbTmp += ( ( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) )*( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) ) );
				if ( dbTmp < dbMin )
				{
					dbMin = dbTmp;
					nLabel = i;
				}
			}

			(*pMap)(r, c) = nLabel;
			if ( nLabelFlag[nLabel] == 0 )
				nLabelFlag[nLabel] = 1;
		}
	
	nLabel = 0;
	for ( i = 0; i < m_nClusterCnt; i++ )
		nLabel += nLabelFlag[i];
	if ( nLabel == m_nClusterCnt )
		bSucFlag = true;
	else
		bSucFlag = false;

	delete[] nLabelFlag;

 //   TGrayImage<int>* pDispMap = new TGrayImage<int>(nWidth, nHeight);
	//for ( int r = 0; r < nHeight; r++ )
	//	for ( int c = 0; c < nWidth; c++ )
	//	{
	//		if ((nLabel = (*pMap)(r, c)) >= 0)
	//		{
	//			if (nLabel == 0)              //QIN_DEBUG: Here I only consider at most three classes
	//				(*pDispMap)(r, c) = 128;
	//			if (nLabel == 1)
	//				(*pDispMap)(r, c) = 0;
	//			if (nLabel == 2)
	//				(*pDispMap)(r, c) = 255;
	//			if (nLabel == 3)
	//				(*pDispMap)(r, c) = 64;
	//			if (nLabel == 4)
	//				(*pDispMap)(r, c) = 192;
	//			if (nLabel == 5)
	//				(*pDispMap)(r, c) = 32;
	//			if (nLabel == 6)
	//				(*pDispMap)(r, c) = 224;
	//		}
	//		else  // boundary label
	//			(*pDispMap)(r, c) = nLabel;
	//	}

	//pDispMap->Save("KM_Pixel_Result.bmp");
	//delete pDispMap;

	return pMap;
}

TGrayImage<int>* TPixelBasedClustering::OTSUClassify(bool& bSucFlag,int threshold)
{
	int nWidth = m_ppSrcImg[0]->Width();
	int nHeight = m_ppSrcImg[0]->Height();
	TGrayImage<int>* pMap = new TGrayImage<int>(nWidth, nHeight);

	// Classification
		int i, nLabel;
		int* nLabelFlag = new int[m_nClusterCnt];
		memset(nLabelFlag, 0, m_nClusterCnt * sizeof(int));

	for (int r = 0; r < nHeight; r ++){
		for (int c = 0; c < nWidth; c ++)
		{
			if (m_pMask && !(*m_pMask)(r, c))
			{
				(*pMap)(r, c) = -2;
				continue;
			}
			if ((*m_ppSrcImg[0])(r, c)>threshold){
				(*pMap)(r, c) = 2;
			if ( nLabelFlag[0] == 0 )
				nLabelFlag[0] = 1;}
			else {
				(*pMap)(r, c) = 1;
			if ( nLabelFlag[1] == 0 )
				nLabelFlag[1] = 1;
				}
			}
		}
nLabel = 0;
	for ( i = 0; i < m_nClusterCnt; i++ )
		nLabel += nLabelFlag[i];
	if ( nLabel == m_nClusterCnt )
		bSucFlag = true;
	else
		bSucFlag = false;

	delete[] nLabelFlag;
	return pMap;
}

/// ************************************************************************************************************************************* ///

/// ********************************************************** KMGMM Clustering ************************************************************///
// 
// KMGMM clustering
// 
// Parameter in
//  nMode: fast mode - no effect in multivariate mode
//
// Return
//  the clustering result
TGrayImage<int>* TPixelBasedClustering::KMGMMClustering(int nInitMode, bool& bSucFlag, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	if (m_nFeatureCnt == 1)
	{
		if (!(nMode & FAST))
			KMGMMParamEst( nInitMode, dbRatio, nInit1, nInit2, nIterFinal );
		else
		{
			// For fast mode, the parameters are based on histogram. This significantly
			// improves the speed. However, it assumes that the features are univariate
			// and of range [0 255], and thus are only useful for intensity feature.
			int *pHist = m_ppSrcImg[0]->GetIntensityHisto( -0.4, 255.6, 256, m_pMask );  // histogram (interesting respresentation)
			KMGMMParamEst( nInitMode, pHist, 256, dbRatio, nInit1, nInit2, nIterFinal );
			delete[] pHist;
		}

		//// Show the parameters
		//ShowParam();

		return GaussMixClassify( bSucFlag );
	} else if (m_nFeatureCnt > 1)
	{
		KMGMMParamEst( nInitMode, dbRatio, nInit1, nInit2, nIterFinal );
		return GaussMixClassify( bSucFlag );
	}

	else
		return NULL;

}

///// ********************************************************** KMGMM Clustering ************************************************************///
//// 
//// KMGMM clustering in univariate feature case.
//// 
//// Parameter in
////  pFeature:  the univariate input feature image
////  pMask: the mask image
////  nClusterCnt: number of classes
////  nMode: mode
////
//// Return
////  the clustering result
//// 
//TGrayImage<int>* TPixelBasedClustering::KMGMMClustering( int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>* pFeature, TMonoImage* pMask, int nClusterCnt, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal )
//{
//	if (!(nMode & FAST))
//		KMGMMParamEst( nInitMode, dbRatio, nInit1, nInit2, nIterFinal );
//	else
//	{
//		// For fast mode, the parameters are based on histogram. This significantly
//		// improves the speed. However, it assumes that the features are univariate
//		// and of range [0 255], and thus are only useful for intensity feature.
//		int *pHist = pFeature->GetIntensityHisto( -0.4, 255.6, 256, pMask );  // histogram (interesting respresentation)
//		KMGMMParamEst( nInitMode, pHist, 256, dbRatio, nInit1, nInit2, nIterFinal );
//		delete[] pHist;
//	}
//
//	//// Show the parameters
//	//ShowParam();
//
//	return GaussMixClassify( bSucFlag );
//}
//
//// 
//// KMGMM clustering in multivariate feature case.
//// 
//// Parameter in
////  ppFeature: the multivariate input feature image
////  pMask: the mask image
////  nClusterCnt: number of classes
////  nFeatureCnt: feature space dimensionality
////
//// Return
////  the clustering result
//// 
//TGrayImage<int>* TPixelBasedClustering::KMGMMClustering( int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>** ppFeatures, int nFeatureCnt, TMonoImage* pMask, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal )
//{
//	KMGMMParamEst( nInitMode, dbRatio, nInit1, nInit2, nIterFinal );
//	return GaussMixClassify( bSucFlag );
//}

// 
// Parameter estimation for KMGMM based on histogram
// 
void TPixelBasedClustering::KMGMMParamEst(int nInitMode, int* pHist, int nHistBinCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	// Allocate memory
	// Alloc(m_nClusterCnt, 1);
	TVector **ppMuBest = AllocVectors(m_nClusterCnt, 1);  // for preserving the best

	// Find the range and variance of each feature
	double* VecMin = new double[1];
	double* VecMax = new double[1];
	double* Variance = new double[1];
	VecMin[0] = 0;
	VecMax[0] = nHistBinCnt - 1;
	
	int i, nCnt = 0;
	double dbMu = 0;
	Variance[0] = 0;
	for ( i = 0; i < nHistBinCnt - 1; i++ )
	{
		dbMu += ( pHist[i] * i );
		nCnt += pHist[i];
	}
	dbMu /= nCnt;

	for ( i = 0; i < nHistBinCnt - 1; i++ )
		Variance[0] += ( pHist[i] * ( i - dbMu ) * ( i - dbMu ) );

	Variance[0] /= nCnt;

	// Estimation
	double dbEnergy, dbMinEnergy = -999999999999;
	int j, nIter = 0;  // iteration number
	do {
		PixelInitParam(nInitMode, pHist, nHistBinCnt, Variance);
		dbEnergy = KMIterations( dbRatio, nInit1, pHist, nHistBinCnt );  // 10 times?

		if (dbMinEnergy < 0 ||  // first iteration
		    dbEnergy < dbMinEnergy)  
			{
				dbMinEnergy = dbEnergy;
				for ( i = 0; i < m_nClusterCnt; i++ )
					ppMuBest[i]->CopyFrom(m_ppMu[i]);
			}
	} while ( ++nIter < nInit2 );

    // Set the parameters with the obtained best and continue with EM based GMM
	for (i = 0; i < m_nClusterCnt; i ++)
	{
		m_ppMu[i]->CopyFrom(ppMuBest[i]);
		delete ppMuBest[i];
	}	
	delete[] ppMuBest;

    // Update other parameters needed for GMM further
	int nTotalCnt = 0, nMinIdx;
	int* nSampleCnt = new int[m_nClusterCnt];
	memset(nSampleCnt, 0, sizeof(int) * m_nClusterCnt);
	double dbTmp, dbMin;
	for ( i = 0; i < nHistBinCnt; i++ )
	{
		dbMin = 999999999999;
		for ( j = 0; j < m_nClusterCnt; j++ )
		{
			dbTmp = ( i - (*m_ppMu[j])(0) )*( i - (*m_ppMu[j])(0) );
			if ( dbTmp < dbMin )
			{
				dbMin = dbTmp;
				nMinIdx = j;
			}
		}
		(*m_ppCovInv[nMinIdx])(0, 0) += ( i * i * pHist[i] );
		nSampleCnt[nMinIdx] += pHist[i];
	}

	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		if ( nSampleCnt[i] > 0 )
			(*m_ppCovInv[i])(0, 0) = (*m_ppCovInv[i])(0, 0) / nSampleCnt[i] -  (*m_ppMu[i])(0) * (*m_ppMu[i])(0);
		else
			nSampleCnt[i]++;

		if ( (*m_ppCovInv[i])(0, 0) <= ZERO_THRESH )
			(*m_ppCovInv[i])(0, 0) += Variance[0];

		m_dbDeterm[i] = sqrt( (*m_ppCovInv[i])(0, 0) );  // determinant
		(*m_ppCovInv[i])(0, 0) = 1 / (*m_ppCovInv[i])(0, 0);  // inverse of covariance

		nTotalCnt += nSampleCnt[i];
	}
    
	for ( i = 0; i < m_nClusterCnt; i++ )
		m_dbPri[i] = ((double) nSampleCnt[i]) / ((double) nTotalCnt);  // prior

	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);
	
	// Further iterations to converge
	dbEnergy = GaussMixEMIterations( dbRatio, nInit1, pHist, nHistBinCnt );

	//printf("The final minimum energy is %.6f\n", dbMinEnergy);

	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	//SortCluster();

	//// Show the parameters
	//ShowParam();
}

// 
// Parameter estimation for KMGMM based on multiple features
// 
void TPixelBasedClustering::KMGMMParamEst(int nInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal)
{
	// Allocate memory
	// Alloc(m_nClusterCnt, m_nFeatureCnt);
	TVector **ppMuBest = AllocVectors(m_nClusterCnt, m_nFeatureCnt);  // for preserving the best

	// Find the range and variance of each feature
	double* VecMin = new double[m_nFeatureCnt];
	double* VecMax = new double[m_nFeatureCnt];
	double* Mu = new double[m_nFeatureCnt];
	double* Variance = new double[m_nFeatureCnt];
	int i;
	for ( i = 0; i < m_nFeatureCnt; i++ )
	{
		((TImage*) m_ppSrcImg[i])->GetIntensityRange(VecMax[i], VecMin[i], m_pMask);
		Mu[i] = m_ppSrcImg[i]->GetIntensityMean(m_pMask);
	}

	int nCnt = 0;
	double dbAvgVar = 0.0;
	for ( i = 0; i < m_nFeatureCnt; i++ )
	{
		Variance[i] = 0;
		nCnt = 0;
		for ( int j = 0; j < m_ppSrcImg[i]->Height(); j++ )
			for ( int k = 0; k < m_ppSrcImg[i]->Width(); k++ )
			{
				if ( m_pMask && ! (*m_pMask)(j, k) )
					continue;
				Variance[i] += ( ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) * ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) );
				nCnt++;
			}
		Variance[i] /= nCnt;
		dbAvgVar += Variance[i];
	}
	dbAvgVar /= m_nFeatureCnt;

	// Parameter estimation
	double dbEnergy, dbMinEnergy = -999999999999;
	int nIter = 0;  // iteration number and sample count
	do {
		PixelInitParam(nInitMode, VecMin, VecMax, Variance);
		dbEnergy = KMIterations( dbRatio, nInit1 ); //10 times?

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
	for ( i = 0; i < m_nClusterCnt; i++ )
	{
		m_ppMu[i]->CopyFrom(ppMuBest[i]);
		delete ppMuBest[i];
	}	
	delete[] ppMuBest;

    // Update Other parameters needed for GMM further
	int j, k, nTotalCnt = 0, nMinIdx;
	int* nSampleCnt = new int[m_nClusterCnt];
	memset(nSampleCnt, 0, sizeof(int) * m_nClusterCnt);
	double dbTmp, dbMin;
	for (int r = m_ppSrcImg[0]->Height() - 1; r >= 0; r --)
		for (int c = m_ppSrcImg[0]->Width() - 1; c >= 0; c --)
		{
			if (m_pMask && !(*m_pMask)(r, c))
				continue;

			dbMin = 999999999999;
			for ( i = 0; i < m_nClusterCnt; i++ )
			{
				dbTmp = 0;
				for ( j = 0; j < m_nFeatureCnt; j++ )
					dbTmp += ( ( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) )*( (*m_ppSrcImg[j])(r, c) - (*m_ppMu[i])(j) ) );
				if ( dbTmp < dbMin )
				{
					dbMin = dbTmp;
					nMinIdx = i;
				}
			}
			for ( i = 0; i < m_nFeatureCnt; i++ )
				for ( j = i ; j < m_nFeatureCnt; j++)
				    (*m_ppCovInv[nMinIdx])(i, j) += ( (*m_ppSrcImg[i])(r, c) * (*m_ppSrcImg[j])(r, c) );
			
			nSampleCnt[nMinIdx]++;
		}

	for ( k = 0; k < m_nClusterCnt; k++ )
	{
		if ( nSampleCnt[k] > 0 )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ )  // covariance
				for ( j = i; j < m_nFeatureCnt; j++ )
					(*m_ppCovInv[k])(i, j) = (*m_ppCovInv[k])(i, j) / nSampleCnt[k] -  (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
		}
		else
			nSampleCnt[k]++;

		for ( i = 0; i < m_nFeatureCnt; i++ )
			for ( j = i + 1; j < m_nFeatureCnt; j++ )
				(*m_ppCovInv[k])(j, i) = (*m_ppCovInv[k])(i, j);

		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += dbAvgVar;
		}

		m_dbDeterm[k] = sqrt( m_ppCovInv[k]->Determinant() );  // determinant
		m_ppCovInv[k]->Inv(m_ppCovInv[k]);  // inverse of covariance

		nTotalCnt += nSampleCnt[k];
	}

	for ( k = 0; k < m_nClusterCnt; k++ )
		m_dbPri[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);  // prior

	//printf("The initial minimum energy is %.6f\n", dbMinEnergy);
	
	// Further iterations to converge
	dbMinEnergy = GaussMixEMIterations( dbRatio, nIterFinal);

	//printf("The final minimum energy is %.6f\n", dbMinEnergy);

	//// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	//SortCluster();

	//// Show the parameters
	//ShowParam();
}

/// ************************************************************************************************************************************* ///

/// ********************************************************** CMLL Pixel-based Clustering ************************************************************///
TGrayImage<int>* TPixelBasedClustering::CMLLClustering( int nIter, int nBetaMode, double dbBeta1, double dbBeta2,
									                    int nInitMode, int nInitInitMode, bool& bInitSucFlag, double dbRatio, int nInit1, int nInit2, int nIterFinal, 
														double dbSAInitT, bool bOutIntermResult )
{
	// Allocate memory
	// Alloc(m_nClusterCnt, m_nFeatureCnt);

	int nHeight = m_ppSrcImg[0]->Height();
	int nWidth = m_ppSrcImg[0]->Width();

	// Initialization
	TGrayImage<int>* pMap;
	if ( nInitMode == RAND1 )
	{
		// Find the range and variance of each feature
		double* VecMin = new double[m_nFeatureCnt];
		double* VecMax = new double[m_nFeatureCnt];
		double* Mu = new double[m_nFeatureCnt];
		double* Variance = new double[m_nFeatureCnt];
		int i;
		for ( i = 0; i < m_nFeatureCnt; i++ )
		{
			((TImage*) m_ppSrcImg[i])->GetIntensityRange(VecMax[i], VecMin[i], m_pMask);
			Mu[i] = m_ppSrcImg[i]->GetIntensityMean(m_pMask);
		}

		int nCnt = 0;
		for ( i = 0; i < m_nFeatureCnt; i++ )
		{
			Variance[i] = 0;
			nCnt = 0;
			for ( int j = 0; j < m_ppSrcImg[i]->Height(); j++ )
				for ( int k = 0; k < m_ppSrcImg[i]->Width(); k++ )
				{
					if ( m_pMask && ! (*m_pMask)(j, k) )
						continue;
					Variance[i] += ( ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) * ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) );
					nCnt++;
				}
			Variance[i] /= nCnt;
		}

		int j;
		for ( i = 0; i < m_nClusterCnt; i++ )
		{
			m_ppCovInv[i]->SetZero();
			m_dbPri[i] = 1.0 / m_nClusterCnt;  // assume the priors equal
		}
		
		double dbRange, dbMin;
		for ( i = 0; i < m_nFeatureCnt; i++ )
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
		for ( i = 0; i < m_nFeatureCnt; i++ )
			dbDeterm *= ( 1 / (*m_ppCovInv[0])(i, i) );
		dbDeterm = sqrt(dbDeterm);   // Note it is the square root of the determinant
		for ( i = 0; i < m_nClusterCnt; i++ )
			m_dbDeterm[i] = dbDeterm;

		delete[] VecMin;
		delete[] VecMax;
		delete[] Mu;
		delete[] Variance;
        
		pMap = GaussMixClassify(bInitSucFlag);
	}

	if ( nInitMode == RAND2 )
	{
		// Find the average variance over all features
		double dbAvgVar = 0.0;
		double* Mu = new double[m_nFeatureCnt];
		double* Variance = new double[m_nFeatureCnt];
		int i;
		for ( i = 0; i < m_nFeatureCnt; i++ )
			Mu[i] = m_ppSrcImg[i]->GetIntensityMean(m_pMask);

		int nCnt = 0;
		for ( i = 0; i < m_nFeatureCnt; i++ )
		{
			Variance[i] = 0;
			nCnt = 0;
			for ( int j = 0; j < m_ppSrcImg[i]->Height(); j++ )
				for ( int k = 0; k < m_ppSrcImg[i]->Width(); k++ )
				{
					if ( m_pMask && ! (*m_pMask)(j, k) )
						continue;
					Variance[i] += ( ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) * ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) );
					nCnt++;
				}
			Variance[i] /= nCnt;
			dbAvgVar += Variance[i];
		}
		
		dbAvgVar /= m_nFeatureCnt;

		pMap = new TGrayImage<int>( nWidth, nHeight );
		GenerateRandomMap(pMap);
		UpdateParam(pMap, dbAvgVar);
		bInitSucFlag = true;
	}


	if ( nInitMode == KM_PIXEL )
	{
		bool bSucFlag;

		TGrayImage<int>* pRetTemp = 0; // when KMClustering is run, it allocates an image. This needs to be deleted.
		pRetTemp = KMClustering(nInitInitMode, bSucFlag, FAST, dbRatio, nInit1, nInit2, nIterFinal);
		if (pRetTemp)
			delete pRetTemp; // this HAS to be deleted.
		
        pMap = GaussMixClassify(bInitSucFlag);
	}

	if ( nInitMode == GMM_PIXEL )
	{
		pMap = GaussMixClustering(nInitInitMode, bInitSucFlag, FAST, dbRatio, nInit1, nInit2, nIterFinal);
	}

	if ( nInitMode == KMGMM_PIXEL )
	{
		pMap = KMGMMClustering(nInitInitMode, bInitSucFlag, FAST, dbRatio, nInit1, nInit2, nIterFinal);
	}

	// Output intermediate results for each iteration
	char szOut[64];
	szOut[0] = 'o';
	szOut[1] = 'u';
	szOut[2] = 't';
	szOut[6] = '.';
	szOut[7] = 'b';
	szOut[8] = 'm';
	szOut[9] = 'p';
	szOut[10] = 0;
	TGrayImage<int>* pDispMap;
	if (bOutIntermResult)
		pDispMap = new TGrayImage<int>( nWidth, nHeight );

	// Refine by using simulated annealing
	TAnnealingUtil* pAnnealUtil = new TAnnealingUtil(dbSAInitT, 0.98);
	TSampler* pScanUtil = new TSampler();
	double TotalEnergy;

	// Compute total number of site involved
	int i, j, nTotalCnt = 0;
	for ( i = 0; i < nHeight; i++ )
		for ( j = 0; j < nWidth; j++ )
		{
			if ( m_pMask && !(*m_pMask)(i, j) )
				continue;
			nTotalCnt ++;
		}
	
	double dbBeta;
	if ( nBetaMode == FIXED )
	{
		dbBeta = dbBeta1;
		cout << "Fixed Beta is: " << dbBeta << endl;
	}

	double dbFisher;
	for ( i = 0; i < nIter; i++ )
	{
		if ( nBetaMode == MLFIXED )
		{
			dbBeta = UpdateBeta( ComputeEdgePairCnt(pMap, m_pMask) / 6.0 / nTotalCnt );
			dbBeta *= dbBeta1;
			cout << "Adjusted Beta is: " << dbBeta << endl;
		}
		
		if ( nBetaMode == MLFISH )
		{
			dbBeta = UpdateBeta( ComputeEdgePairCnt(pMap, m_pMask) / 6.0 / nTotalCnt );
		    dbFisher = ComputeMinFisher();
			dbBeta *= ( dbBeta1 * dbFisher / dbBeta2 / (1 + dbFisher / dbBeta2) );
			cout << "Fisher: " << dbFisher << ", Adjusted Beta: " << dbBeta << endl;
		}

		TotalEnergy = CMLLAnnealing(pMap, dbBeta, pAnnealUtil, pScanUtil);
		
		cout << "Annealing iteration " << i << " is finished: energy = " << TotalEnergy << endl;

		// compute the next temparature
		pAnnealUtil->NextTemperature();
		//if ( pAnnealUtil->GetTemparature() < 0.1 )
			if ( !UpdateParam(pMap, ZERO_THRESH ) )
				goto END;

		if (bOutIntermResult)
		{
			// Output intermediate results for each iteration
			int k = i, l;
			for ( j = 2; j >= 0; j-- )  // 4 digits
			{
				l = k % 10;
				k = k / 10;
				szOut[3 + j] = '0' + l;
			}


			for ( int r = 0; r < nHeight; r++ )
				for ( int c = 0; c < nWidth; c++ )
				{
					if ( ( l = (*pMap)(r, c) ) >= 0 )  // class label
					{
						if (l == 0)
							(*pDispMap)(r, c) = 128;
						if (l == 1)
							(*pDispMap)(r, c) = 0;
						if (l == 2)
							(*pDispMap)(r, c) = 255;
						if (l == 3)
							(*pDispMap)(r, c) = 64;
						if (l == 4)
							(*pDispMap)(r, c) = 192;
						if (l == 5)
							(*pDispMap)(r, c) = 32;
						if (l == 6)
							(*pDispMap)(r, c) = 224;
					}
					else  // boundary label
						(*pDispMap)(r, c) = l;
				}
//			pDispMap->Save(szOut);
		}

	}

END:
	if (bOutIntermResult)
		delete pDispMap;
	delete pAnnealUtil;
	delete pScanUtil;
	
	return pMap;
}

// 
// Simulated annealing for CMLL
// 
// Parameter in
//  pMap: the previous segmentation map
//  pMask: the mask
//  ppFeatures: the input multivariate image
//  nFeatureCnt: number of features
//  nClusterCnt: number of clusters
//  dbBeta: the MLL penalty weight
//  pAnnealUtil: utility class for simulated annealing
//  pScanUtil: utility class for sampling
//
// Parameter out
//  pMap: the current segmentation map
// 
// Return
//  the total energy
/////////////////////////////////////////////////////////////////////////////
double TPixelBasedClustering::CMLLAnnealing(TGrayImage<int>*& pMap, double dbBeta, 
					 TAnnealingUtil* pAnnealUtil, TSampler* pScanUtil)
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int i;
	double* dbLnSigma = new double[m_nClusterCnt];
	for ( i = 0; i < m_nClusterCnt; i++ )
		dbLnSigma[i] = log(m_dbDeterm[i]);

	double dbTotalEnergy = 0;
	double* dbEnergies = new double[m_nClusterCnt];
	TVector* pVecFeature = new TVector(m_nFeatureCnt, false);  // feature vector
	TVector* pVecTemp1 = new TVector(m_nFeatureCnt, false);  // temporary feature vector
	TVector* pVecTemp2 = new TVector(m_nFeatureCnt, false);  // temporary feature vector
	int r, c, nRow, nCol, k;
	pScanUtil->InitScan(nWidth * nHeight);
	while ( ( k = pScanUtil->GetNextScan() ) >= 0 )
	{
		nRow = k / nWidth;
		nCol = k % nWidth;

		if ( m_pMask && !(*m_pMask)(nRow, nCol) )
			continue;

		for ( i = 0; i < m_nFeatureCnt; i++ )
			(*pVecFeature)(i) = (*m_ppSrcImg[i])(nRow, nCol);

		for ( k = 0; k < m_nClusterCnt; k++ )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ )
				(*pVecTemp1)(i) = (*pVecFeature)(i) - (*m_ppMu[k])(i);
			pVecTemp1->Mul(pVecTemp2, m_ppCovInv[k]);
			dbEnergies[k] = pVecTemp2->VecMul(pVecTemp1) / 2 + dbLnSigma[k];
		}

		for (i = 0; i < 8; i ++)
		{
			r = nRow + neighbor8[i].Row;
			c = nCol + neighbor8[i].Col;
			if ( r < 0 || r >= nHeight || c < 0 || c >= nWidth )
				continue;
			if ( m_pMask && !(*m_pMask)(r, c) )
				continue;

			k = (*pMap)(r, c);
			dbEnergies[k] -= dbBeta;
		}

		(*pMap)(nRow, nCol) = pAnnealUtil->GibbsSampler(dbEnergies, m_nClusterCnt);
		dbTotalEnergy += ( dbEnergies[(*pMap)(nRow, nCol)] + 8 * dbBeta );
	}

	delete[] dbLnSigma ;
	delete[] dbEnergies;
	delete pVecFeature;
	delete pVecTemp1;
	delete pVecTemp2;

	return dbTotalEnergy;
}


/// ************************************************************************************************************************************* ///

/// ********************************************************** VMLL Pixel-based Clustering ************************************************************///
TGrayImage<int>* TPixelBasedClustering::VMLLClustering( int nIter, 
		                                                int nBetaMode, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2,
									                    int nInitMode, int nInitInitMode, bool& bInitSucFlag, double dbRatio, int nInit1, int nInit2, int nIterFinal,
														double dbSAInitT, bool bOutIntermResult )
{
	// Allocate memory
	// Alloc(m_nClusterCnt, m_nFeatureCnt);

	int nHeight = m_ppSrcImg[0]->Height();
	int nWidth = m_ppSrcImg[0]->Width();

	// Initialization
	TGrayImage<int>* pMap;
	if ( nInitMode == RAND1 )
	{
		// Find the range and variance of each feature
		double* VecMin = new double[m_nFeatureCnt];
		double* VecMax = new double[m_nFeatureCnt];
		double* Mu = new double[m_nFeatureCnt];
		double* Variance = new double[m_nFeatureCnt];
		int i;
		for ( i = 0; i < m_nFeatureCnt; i++ )
		{
			((TImage*) m_ppSrcImg[i])->GetIntensityRange(VecMax[i], VecMin[i], m_pMask);
			Mu[i] = m_ppSrcImg[i]->GetIntensityMean(m_pMask);
		}

		int nCnt = 0;
		for ( i = 0; i < m_nFeatureCnt; i++ )
		{
			Variance[i] = 0;
			nCnt = 0;
			for ( int j = 0; j < m_ppSrcImg[i]->Height(); j++ )
				for ( int k = 0; k < m_ppSrcImg[i]->Width(); k++ )
				{
					if ( m_pMask && ! (*m_pMask)(j, k) )
						continue;
					Variance[i] += ( ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) * ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) );
					nCnt++;
				}
			Variance[i] /= nCnt;
		}

		int j;
		for ( i = 0; i < m_nClusterCnt; i++ )
		{
			m_ppCovInv[i]->SetZero();
			m_dbPri[i] = 1.0 / m_nClusterCnt;  // assume the priors equal
		}
		
		double dbRange, dbMin;
		for ( i = 0; i < m_nFeatureCnt; i++ )
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
		for ( i = 0; i < m_nFeatureCnt; i++ )
			dbDeterm *= ( 1 / (*m_ppCovInv[0])(i, i) );
		dbDeterm = sqrt(dbDeterm);   // Note it is the square root of the determinant
		for ( i = 0; i < m_nClusterCnt; i++ )
			m_dbDeterm[i] = dbDeterm;

		delete[] VecMin;
		delete[] VecMax;
		delete[] Mu;
		delete[] Variance;
        
		pMap = GaussMixClassify(bInitSucFlag);
	}

	if ( nInitMode == RAND2 )
	{
		// Find the average variance over all features
		double dbAvgVar = 0.0;
		double* Mu = new double[m_nFeatureCnt];
		double* Variance = new double[m_nFeatureCnt];
		int i;
		for ( i = 0; i < m_nFeatureCnt; i++ )
			Mu[i] = m_ppSrcImg[i]->GetIntensityMean(m_pMask);

		int nCnt = 0;
		for ( i = 0; i < m_nFeatureCnt; i++ )
		{
			Variance[i] = 0;
			nCnt = 0;
			for ( int j = 0; j < m_ppSrcImg[i]->Height(); j++ )
				for ( int k = 0; k < m_ppSrcImg[i]->Width(); k++ )
				{
					if ( m_pMask && ! (*m_pMask)(j, k) )
						continue;
					Variance[i] += ( ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) * ( (*m_ppSrcImg[i])(j, k) - Mu[i] ) );
					nCnt++;
				}
			Variance[i] /= nCnt;
			dbAvgVar += Variance[i];
		}
		
		dbAvgVar /= m_nFeatureCnt;

		pMap = new TGrayImage<int>( nWidth, nHeight );
		GenerateRandomMap(pMap);
		UpdateParam(pMap, dbAvgVar);
		bInitSucFlag = true;
	}

	if ( nInitMode == KM_PIXEL )
	{
		bool bSucFlag;
		TGrayImage<int>* pRetTemp = 0; // when KMClustering is run, it allocates an image. This needs to be deleted.
		pRetTemp = KMClustering(nInitInitMode, bSucFlag, FAST, dbRatio, nInit1, nInit2, nIterFinal);
		if (pRetTemp)
			delete pRetTemp; // this HAS to be deleted.

        pMap = GaussMixClassify(bInitSucFlag);
	}

	if ( nInitMode == GMM_PIXEL )
	{
		pMap = GaussMixClustering(nInitInitMode, bInitSucFlag, FAST, dbRatio, nInit1, nInit2, nIterFinal);
	}

	if ( nInitMode == KMGMM_PIXEL )
	{
		pMap = KMGMMClustering(nInitInitMode, bInitSucFlag, FAST, dbRatio, nInit1, nInit2, nIterFinal);
	}

	// Output intermediate results for each iteration
	char szOut[64];
	szOut[0] = 'o';
	szOut[1] = 'u';
	szOut[2] = 't';
	szOut[6] = '.';
	szOut[7] = 'b';
	szOut[8] = 'm';
	szOut[9] = 'p';
	szOut[10] = 0;
	TGrayImage<int>* pDispMap;
	if (bOutIntermResult)
		pDispMap = new TGrayImage<int>( nWidth, nHeight );

	// Refine by using simulated annealing
	TAnnealingUtil* pAnnealUtil = new TAnnealingUtil(dbSAInitT, 0.98);
	TSampler* pScanUtil = new TSampler();
	double TotalEnergy;

	// Compute total number of site involved
	int i, j, nTotalCnt = 0;
	for ( i = 0; i < nHeight; i++ )
		for ( j = 0; j < nWidth; j++ )
		{
			if ( m_pMask && !(*m_pMask)(i, j) )
				continue;
			nTotalCnt ++;
		}
	
	double dbBeta, dbAlpha;
	if ( nBetaMode == FIXED )
	{
		dbBeta = dbBeta1;
		cout << "Fixed Beta is: " << dbBeta << endl;
	}

	double dbFisher;
	for ( i = 0; i < nIter; i++ )
	{
		dbAlpha = UpdateAlpha( nIter, dbAlpha1, dbAlpha2 );

		if ( nBetaMode == MLFIXED )
		{
			dbBeta = UpdateBeta( ComputeEdgePairCnt(pMap, m_pMask) / 6.0 / nTotalCnt );
			dbBeta *= dbBeta1;
			cout << "Adjusted Beta is: " << dbBeta << endl;
		}
		
		if ( nBetaMode == MLFISH )
		{
			dbBeta = UpdateBeta( ComputeEdgePairCnt(pMap, m_pMask) / 6.0 / nTotalCnt );
		    dbFisher = ComputeMinFisher();
			dbBeta *= ( dbBeta1 * dbFisher / dbBeta2 / (1 + dbFisher / dbBeta2) );
			cout << "Fisher: " << dbFisher << ", Adjusted Beta: " << dbBeta << endl;
		}

		TotalEnergy = VMLLAnnealing(pMap, dbBeta, dbAlpha, pAnnealUtil, pScanUtil);
		
		cout << "Annealing iteration " << i << " is finished: energy = " << TotalEnergy << endl;

		// compute the next temparature
		pAnnealUtil->NextTemperature();
		//if ( pAnnealUtil->GetTemparature() < 0.1 )
			if ( !UpdateParam(pMap, ZERO_THRESH ) )
				goto END;

		if (bOutIntermResult)
		{
			// Output intermediate results for each iteration
			int k = i, l;
			for ( j = 2; j >= 0; j-- )  // 4 digits
			{
				l = k % 10;
				k = k / 10;
				szOut[3 + j] = '0' + l;
			}


			for ( int r = 0; r < nHeight; r++ )
				for ( int c = 0; c < nWidth; c++ )
				{
					if ( ( l = (*pMap)(r, c) ) >= 0 )  // class label
					{
						if (l == 0)
							(*pDispMap)(r, c) = 128;
						if (l == 1)
							(*pDispMap)(r, c) = 0;
						if (l == 2)
							(*pDispMap)(r, c) = 255;
						if (l == 3)
							(*pDispMap)(r, c) = 64;
						if (l == 4)
							(*pDispMap)(r, c) = 192;
						if (l == 5)
							(*pDispMap)(r, c) = 32;
						if (l == 6)
							(*pDispMap)(r, c) = 224;
					}
					else  // boundary label
						(*pDispMap)(r, c) = l;
				}
//			pDispMap->Save(szOut);
		}

	}

END:
	if (bOutIntermResult)
		delete pDispMap;

	delete pAnnealUtil;
	delete pScanUtil;
	
	return pMap;
}

// 
// Simulated annealing for VMLL
// 
// Parameter in
//  pMap: the previous segmentation map
//  pMask: the mask
//  ppFeatures: the input multivariate image
//  nFeatureCnt: number of features
//  nClusterCnt: number of clusters
//  dbAlpha: the feature model weight
//  dbBeta: the MLL penalty weight
//  pAnnealUtil: utility class for simulated annealing
//  pScanUtil: utility class for sampling
//
// Parameter out
//  pMap: the current segmentation map
// 
// Return
//  the total energy
/////////////////////////////////////////////////////////////////////////////
double TPixelBasedClustering::VMLLAnnealing(TGrayImage<int>*& pMap, double dbAlpha, double dbBeta, 
					 TAnnealingUtil* pAnnealUtil, TSampler* pScanUtil)
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int i;
	double* dbLnSigma = new double[m_nClusterCnt];
	for ( i = 0; i < m_nClusterCnt; i++ )
		dbLnSigma[i] = log(m_dbDeterm[i]);

	double dbTotalEnergy = 0;
	double* dbEnergies = new double[m_nClusterCnt];
	TVector* pVecFeature = new TVector(m_nFeatureCnt, false);  // feature vector
	TVector* pVecTemp1 = new TVector(m_nFeatureCnt, false);  // temporary feature vector
	TVector* pVecTemp2 = new TVector(m_nFeatureCnt, false);  // temporary feature vector
	int r, c, nRow, nCol, k;
	pScanUtil->InitScan(nWidth * nHeight);
	while ( ( k = pScanUtil->GetNextScan() ) >= 0 )
	{
		nRow = k / nWidth;
		nCol = k % nWidth;

		if ( m_pMask && !(*m_pMask)(nRow, nCol) )
			continue;

		for ( i = 0; i < m_nFeatureCnt; i++ )
			(*pVecFeature)(i) = (*m_ppSrcImg[i])(nRow, nCol);

		for ( k = 0; k < m_nClusterCnt; k++ )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ )
				(*pVecTemp1)(i) = (*pVecFeature)(i) - (*m_ppMu[k])(i);
			pVecTemp1->Mul(pVecTemp2, m_ppCovInv[k]);
			dbEnergies[k] = dbAlpha * ( pVecTemp2->VecMul(pVecTemp1) / 2 + dbLnSigma[k] );
		}

		for (i = 0; i < 8; i ++)
		{
			r = nRow + neighbor8[i].Row;
			c = nCol + neighbor8[i].Col;
			if ( r < 0 || r >= nHeight || c < 0 || c >= nWidth )
				continue;
			if ( m_pMask && !(*m_pMask)(r, c) )
				continue;

			k = (*pMap)(r, c);
			dbEnergies[k] -= dbBeta;
		}

		(*pMap)(nRow, nCol) = pAnnealUtil->GibbsSampler(dbEnergies, m_nClusterCnt);
		dbTotalEnergy += ( dbEnergies[(*pMap)(nRow, nCol)] + 8 * dbBeta );
	}

	delete[] dbLnSigma ;
	delete[] dbEnergies;
	delete pVecFeature;
	delete pVecTemp1;
	delete pVecTemp2;

	return dbTotalEnergy;
}

//
// Update feature model weight for VMLL
//
// Param in
// nIter:    current iteration number
// dbAlpha1: C1 in the updating rule
// dbAlpha2: Ratio in the updating rule
// Updating rule: C1*Ratio^nIter + 1;
//
// Return
// The weight

double TPixelBasedClustering::UpdateAlpha(int nIter, double dbAlpha1, double dbAlpha2)
{
	double dbAlpha;

	dbAlpha = dbAlpha1 * pow(dbAlpha2, nIter) + 1.0 / m_nFeatureCnt;
//	cout << "The Alpha is " << dbAlpha<< endl;
	return dbAlpha;
}

//
// Compute the number of pairs of neighboring sites belonging to different
// classes
// 
// Param in
//  pMap: the image
//  pMask: the mask
/////////////////////////////////////////////////////////////////////////////
int TPixelBasedClustering::ComputeEdgePairCnt(TGrayImage<int>* pMap, TMonoImage* pMask)
{
	int** data = pMap->GetData();
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();

	int r, c, nLabel, nEdgeCnt = 0;
	for (int nRow = 0; nRow < nHeight; nRow++ )
		for (int nCol = 0; nCol < nWidth; nCol++ )
		{
			if ( pMask && !(*pMask)(nRow, nCol) )
				continue;
			nLabel = data[nRow][nCol];
			for ( int j = 0; j < 8; j++ )
			{
				r = nRow + neighbor8[j].Row;
				c = nCol + neighbor8[j].Col;
				if ( r < 0 || r >= nHeight || c < 0 || c >= nWidth )
					continue;
				if ( pMask && !(*pMask)(r, c) )
					continue;
				if ( data[r][c] != nLabel )
					nEdgeCnt ++;
			}
		}

	return nEdgeCnt;
}

//
// Get weight for edge penalty
// 
// Return
//  The weight
// 
double TPixelBasedClustering::UpdateBeta(double dbEdgeRatio)
{
	double dbBeta;
	if (dbEdgeRatio >= betaLUT[0][0])  // larger than all recorded, select the smallest beta
		dbBeta = betaLUT[0][1];
	else
	{
		int i;
		for (i = 1; i < BETA_CNT; i ++)
		{
			if (dbEdgeRatio > betaLUT[i][0])
				break;
		}

		if (i >= BETA_CNT)  // smaller than all recorded, select the largest beta
			dbBeta = betaLUT[BETA_CNT - 1][1];
		else
			dbBeta = betaLUT[i - 1][1] + 
			         (betaLUT[i][1] - betaLUT[i - 1][1]) /
			         (betaLUT[i][0] - betaLUT[i - 1][0]) *
			         (dbEdgeRatio - betaLUT[i - 1][0]);  // interpolate
	}

//	cout << "The ratio of edge pixels is " << dbEdgeRatio << endl;
//	cout << "The beta is " << dbBeta << endl;
	return dbBeta;
}

// 
// Estimate the parameters of the Gaussian distributed clusters using the obtained
// classification map
// 
// Parameter in
//  pMap: the obtained classification map
//  pMask: the mask
//  ppFeatures: the input multivariate image
// 
bool TPixelBasedClustering::UpdateParam( TGrayImage<int>* pMap, double dbRegCoef)
{
	bool bUpdateFlag = true;

	SetZeroVector(m_ppMu, m_nClusterCnt);
	SetZeroMatrix(m_ppCovInv, m_nClusterCnt);
	memset(m_dbPri, 0, sizeof(double) * m_nClusterCnt);

	int i, j, k;
	for ( int r = 0; r < pMap->Height(); r++ )
		for ( int c = 0; c < pMap->Width(); c++ )
		{
			if ( m_pMask && !(*m_pMask)(r, c) )
				continue;
			if ( ( k = (*pMap)(r, c) ) < 0 )
				continue;

			for ( i = 0; i < m_nFeatureCnt; i++ )
			{
				(*m_ppMu[k])(i) += (*m_ppSrcImg[i])(r, c);
				(*m_ppCovInv[k])(i, i) += ( (*m_ppSrcImg[i])(r, c) * (*m_ppSrcImg[i])(r, c) );
				for ( j = i + 1; j < m_nFeatureCnt; j++ )
					(*m_ppCovInv[k])(i, j) += ( (*m_ppSrcImg[i])(r, c) * (*m_ppSrcImg[j])(r, c) );
			}

			m_dbPri[k]++;
		}

	double dbPriTotal = 0;
	// Compute the parameters
	for ( k = 0; k < m_nClusterCnt; k++ )
	{
		m_nClassPixelCnt[k] = (int) m_dbPri[k]; // store this for later.

		if ( m_dbPri[k] > 0 )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ ) 
				(*m_ppMu[k])(i) /= m_dbPri[k];  // mean
			for ( i = 0; i < m_nFeatureCnt; i++ )  // covariance
			{
				(*m_ppCovInv[k])(i, i) = (*m_ppCovInv[k])(i, i) / m_dbPri[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(i);
				for ( j = i + 1; j < m_nFeatureCnt; j++ )
				{
					(*m_ppCovInv[k])(i, j) = (*m_ppCovInv[k])(i, j) / m_dbPri[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
					(*m_ppCovInv[k])(j, i) = (*m_ppCovInv[k])(i, j);
				}
			}
		}
		else
		{
			bUpdateFlag = false;
			m_dbPri[k]++;
		}

		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
		{
			for ( i = 0; i < m_nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += dbRegCoef;
		}
		m_dbDeterm[k] = sqrt(m_ppCovInv[k]->Determinant());  // determinant
		m_ppCovInv[k]->Inv(m_ppCovInv[k]);  // inverse of covariance

		dbPriTotal += m_dbPri[k];
	}

	for ( k = 0; k < m_nClusterCnt; k++ )
		m_dbPri[k] = m_dbPri[k] / dbPriTotal;  // prior
	
	return bUpdateFlag;
}


// 
// Compute the minimum Fisher criterion.
// 
double TPixelBasedClustering::ComputeMinFisher()
{
	double dbMin = 100000000;  // a big enough number
	double dbVal;
	int i, j;
	if ( m_nFeatureCnt == 1 )
		for (i = 0; i < m_nClusterCnt - 1; i ++)
	    {
		    for (j = i + 1; j < m_nClusterCnt; j ++)
		    {
			    dbVal = (*m_ppMu[i])(0) - (*m_ppMu[j])(0);
			    dbVal *= dbVal;
			    dbVal /= 1 / (*m_ppCovInv[i])(0, 0) * m_dbPri[i] +
			             1 / (*m_ppCovInv[j])(0, 0) * m_dbPri[j];
			    if (dbVal < dbMin)
				    dbMin = dbVal;
		    }
	    }
	else if ( m_nFeatureCnt > 1 ) //QIN_ADD
	{
		TMatrix* WSMat1 = new TMatrix(m_nFeatureCnt, m_nFeatureCnt);
		TMatrix* WSMat2 = new TMatrix(m_nFeatureCnt, m_nFeatureCnt);
		TVector* BSVec = new TVector(m_nFeatureCnt, false);
		TVector* BSVec_T = new TVector(m_nFeatureCnt, false);
		TMatrix *INVWSMat;
		for ( i = 0; i < m_nClusterCnt - 1; i++ )
	    {
		    for ( j = i + 1; j < m_nClusterCnt; j++ )
		    {
				for ( int k = 0; k < m_nFeatureCnt; k++ )
					(*BSVec)(k) = ((*m_ppMu[i])(k) - (*m_ppMu[j])(k));
				
				m_ppCovInv[i]->Inv(WSMat1);
				*WSMat1 *= m_dbPri[i];
				m_ppCovInv[j]->Inv(WSMat2);
				*WSMat2 *= m_dbPri[j];
				(*WSMat1) += (*WSMat2);

				INVWSMat = WSMat1->Inv();

				BSVec->Mul(BSVec_T, INVWSMat);
	            dbVal = BSVec_T->VecMul(BSVec);

				if (dbVal < 0)
				{
					printf("Fisher calculation wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
					//_getchar_nolock();
				}

				if (dbVal < dbMin)
				      dbMin = dbVal;
			}
		}
		delete WSMat1;
		delete WSMat2;
		delete INVWSMat;
		delete BSVec;
		delete BSVec_T;
	}
	else
	{
	    cout << "Error: The number of features should be above 0!" << endl;
		//_getchar_nolock();
		exit(1);
	}

	return dbMin;
}

/// Get segmentation statistics object that lists all the class statistics.
TSegmentationStatistics* TPixelBasedClustering::GetSegStats(TGrayImage<int>* pMap)
{

	if (pMap)
	{
		// do an update based on the map.
		UpdateParam(pMap, ZERO_THRESH);
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
