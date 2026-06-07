// EstimateUtil.cpp : implementation of the TEstimateUtil class
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "EstimateUtil.h"
#include "Matrix.h"
//#include "matrix_gabor.h"   // Use sorting function inside

#define Pi 3.1415926535897
#define MIN_VAR 4.0
#define BETA 1

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

//
// Constructor for EstimateUtil
// 
TEstimateUtil::TEstimateUtil()
{
	m_nClassCnt = 0;
	m_nFeatureCnt = 0;

	m_ppMu = 0;
	m_ppCovInv = 0;
	m_dbDeterm = 0;
	m_dbPri = 0;

	m_pVecTemp1 = 0;
	m_pVecTemp2 = 0;

//	srand((unsigned) time(NULL));
}


//
// Destructor
// 
TEstimateUtil::~TEstimateUtil()
{
	Clear();
}

//
// Allocate memories needed
// 
void TEstimateUtil::Alloc(int nClassCnt, int nFeatureCnt)
{
	if (m_nClassCnt == nClassCnt && m_nFeatureCnt == nFeatureCnt)
		return;
	Clear();

	m_nClassCnt = nClassCnt;
	m_nFeatureCnt = nFeatureCnt;

	m_ppMu = AllocVectors(nClassCnt, nFeatureCnt);
	m_ppCovInv = AllocMatrice(nClassCnt, nFeatureCnt);
	m_dbDeterm = new double[nClassCnt];
	m_dbPri = new double[nClassCnt];

	m_pVecTemp1 = new TVector(nFeatureCnt, false);  // a temperoray vector during the computation
	m_pVecTemp2 = new TVector(nFeatureCnt, false);  // a temperoray vector during the computation
}

// 
// Clear memory
// 
void TEstimateUtil::Clear()
{
	m_nClassCnt = 0;
	m_nFeatureCnt = 0;

	if (m_ppMu)
	{
		for ( int i = 0; i < m_nClassCnt; i++ )
			delete m_ppMu[i];
		delete[] m_ppMu;
		m_ppMu = 0;
	}
	if (m_ppCovInv)
	{
		for ( int i = 0; i < m_nClassCnt; i++ )
			delete m_ppCovInv[i];
		delete[] m_ppCovInv;
		m_ppCovInv = 0;
	}
	if (m_dbDeterm)
	{
		delete[] m_dbDeterm;
		m_dbDeterm = 0;
	}
	if (m_dbPri)
	{
		delete[] m_dbPri;
		m_dbPri = 0;
	}
	if (m_pVecTemp1)
	{
		delete m_pVecTemp1;
		m_pVecTemp1 = 0;
	}
	if (m_pVecTemp2)
	{
		delete m_pVecTemp2;
		m_pVecTemp2 = 0;
	}
}

//
// Allocate vectors (for means) for a given number of classes
// 
TVector** TEstimateUtil::AllocVectors(int nClassCnt, int nFeatureCnt)
{
	TVector** ppRet = 0;
	if ( nClassCnt >= 1 && nFeatureCnt >= 1 )
	{
		ppRet = new TVector*[nClassCnt];
		for ( int i = 0; i < nClassCnt; i++ )
			ppRet[i] = new TVector(nFeatureCnt, false);  // row vector
	}
	return ppRet;
}

//
// Allocate matrice (for covariance) for a given number of classes
// 
TMatrix** TEstimateUtil::AllocMatrice(int nClassCnt, int nFeatureCnt)
{
	TMatrix** ppRet = 0;
	if ( nClassCnt >= 1 && nFeatureCnt >= 1 )
	{
		ppRet = new TMatrix*[nClassCnt];
		for ( int i = 0; i < nClassCnt; i++ )
			ppRet[i] = new TMatrix(nFeatureCnt, nFeatureCnt);
	}
	return ppRet;
}

// 
// Set all elements in the vectors to be zero
// 
void TEstimateUtil::SetZeroVector(TVector** ppVectors, int nClassCnt)
{
	for ( int i = 0; i < nClassCnt; i++ )
		ppVectors[i]->SetZero();
}

// 
// Set all elements in the matrice to be zero
// 
void TEstimateUtil::SetZeroMatrix(TMatrix** ppMatrice, int nClassCnt)
{
	for ( int i = 0; i < nClassCnt; i++ )
		ppMatrice[i]->SetZero();
}

//
// Arrange the clusters in the ascending order of the first mean feature 
// values, so that the resulting map looks more pleasant
// 
void TEstimateUtil::SortCluster()
{
	TVector* pVecSwap;
	TMatrix* pMatSwap;
	double dbPriSwap, dbDetermSwap;
	for ( int i = 0; i < m_nClassCnt; i++ )
		for ( int j = i + 1; j < m_nClassCnt; j++ )
		{
			if ( (*m_ppMu[i])(0) > (*m_ppMu[j])(0) )
			{
				pVecSwap = m_ppMu[i]; m_ppMu[i] = m_ppMu[j]; m_ppMu[j] = pVecSwap;
				pMatSwap = m_ppCovInv[i]; m_ppCovInv[i] = m_ppCovInv[j]; m_ppCovInv[j] = pMatSwap;
				dbPriSwap = m_dbPri[i]; m_dbPri[i] = m_dbPri[j]; m_dbPri[j] = dbPriSwap;
				dbDetermSwap = m_dbDeterm[i]; m_dbDeterm[i] = m_dbDeterm[j]; m_dbDeterm[j] = dbDetermSwap;
			}
		}
}

// 
// Clustering on Gaussian distributed class mixture. Univariate feature case.
// 
// Parameter in
//  pImg:  the univariate input feature image
//  pMask: the mask image
//  nClassCnt: number of classes
//  nMode: mode
//
// Return
//  the clustering result
// 
TGrayImage<int>* TEstimateUtil::GaussMixClustering(TGrayImage<FLOAT>* pImg, TMonoImage* pMask, int nClassCnt, int nMode)
{
	GaussMixParamEst(pImg, pMask, nClassCnt, nMode);
	return GaussMixClassify(&pImg, pMask);
}

// 
// Clustering on Gauss distributed class mixture. Multivariate feature case.
// 
// Parameter in
//  ppImg: the multivariate input feature image
//  pMask: the mask image
//  nClassCnt: number of classes
//  nFeatureCnt: feature space dimensionality
//
// Return
//  the clustering result
// 
TGrayImage<int>* TEstimateUtil::GaussMixClustering(TGrayImage<FLOAT>** ppImg, TMonoImage* pMask, int nClassCnt, int nFeatureCnt)
{
	GaussMixParamEst(ppImg, pMask, nClassCnt, nFeatureCnt);
	return GaussMixClassify(ppImg, pMask);
}

// 
// Parameter estimation for Gaussian mixture using expectation maximization based on single feature.
// 
void TEstimateUtil::GaussMixParamEst(TGrayImage<FLOAT>* pFeature, TMonoImage* pMask, int nClassCnt, int nMode)
{
	if (!(nMode & FAST))
	{
		GaussMixParamEst(&pFeature, pMask, 1, nClassCnt);
		return;
	}

	// For fast mode, the parameter is based on histogram. This significantly
	// improves the speed. However, it assumes that the features are univariate
	// and is of range [0 255], and thus is only useful for intensity feature.
	int *pHist = pFeature->GetIntensityHisto(-0.4, 255.6, 256, pMask);  // histogram (interesting respresentation)
	GaussMixParamEst(pHist, 256, nClassCnt);

	//// Show the parameters
	//ShowGaussianMixParam();

	delete[] pHist;
}

// 
// Parameter estimation for Gaussian mixture using expectation maximization based on histogram
// 
void TEstimateUtil::GaussMixParamEst(int* pHist, int nHistBinCnt, int nClassCnt)
{
	// Allocate memory
	Alloc(nClassCnt, 1);
	TVector **ppMuBest = AllocVectors(nClassCnt, 1);  // for preserving the best
	TMatrix **ppCovInvBest = AllocMatrice(nClassCnt, 1);  // for preserving the best
	double* dbPriBest = new double[nClassCnt];  // for preserving the best
	double* dbDetermBest = new double[nClassCnt];  // for preserving the best

	// Find the range of each feature
	TVector VecMin(1, false);
	TVector VecMax(1, false);
	VecMin(0) = 0;
	VecMax(0) = nHistBinCnt - 1;

	// Estimation
	double dbEnergy, dbMinEnergy = -1;
	int i, nSampleCnt = 0;  // iteration number and sample count
	do {
		InitGaussMixParam(VecMin, VecMax);  // initialize randomly the parameters
		dbEnergy = GaussMixEMIterations(0.00001, 20, pHist, nHistBinCnt);

		if (dbMinEnergy < 0 ||  // first iteration
		    dbEnergy < dbMinEnergy)  
			{
				dbMinEnergy = dbEnergy;
				for (i = 0; i < nClassCnt; i ++)
				{
					ppMuBest[i]->CopyFrom(m_ppMu[i]);
					ppCovInvBest[i]->CopyFrom(m_ppCovInv[i]);
				}
				memcpy(dbPriBest, m_dbPri, nClassCnt * sizeof(double));
				memcpy(dbDetermBest, m_dbDeterm, nClassCnt * sizeof(double));
			}
	} while (nSampleCnt ++ < 60);

	// Set the parameters with the obtained best and continue with the EM
	for (i = 0; i < nClassCnt; i ++)
	{
		m_ppMu[i]->CopyFrom(ppMuBest[i]);
		m_ppCovInv[i]->CopyFrom(ppCovInvBest[i]);
		cout << (*m_ppMu[i])(0) << ",";
	}
	cout << endl;
	memcpy(m_dbPri, dbPriBest, nClassCnt * sizeof(double));
	memcpy(m_dbDeterm, dbDetermBest, nClassCnt * sizeof(double));

	// Further iterations to converge
	dbEnergy = GaussMixEMIterations(0.000001, 80, pHist, nHistBinCnt);

	// Arrange the clusters in the ascending order of the first mean feature 
	// values, so that the resulting map looks more pleasant
	SortCluster();

	for (i = 0; i < nClassCnt; i ++)
	{
		delete ppMuBest[i];
		delete ppCovInvBest[i];
	}
	delete[] ppMuBest;
	delete[] ppCovInvBest;
	delete[] dbPriBest;
	delete[] dbDetermBest;
}

// 
// Parameter estimation for Gaussian mixture using expectation maximization based on multiple features
// 
void TEstimateUtil::GaussMixParamEst(TGrayImage<FLOAT>** ppFeatures, TMonoImage* pMask,
                                     int nClassCnt, int nFeatureCnt)
{
	// Allocate memory
	Alloc(nClassCnt, nFeatureCnt);
	TVector **ppMuBest = AllocVectors(nClassCnt, nFeatureCnt);  // for preserving the best
	TMatrix **ppCovInvBest = AllocMatrice(nClassCnt, nFeatureCnt);  // for preserving the best
	double* dbPriBest = new double[nClassCnt];  // for preserving the best
	double* dbDetermBest = new double[nClassCnt];  // for preserving the best

	// Find the range of each feature
	TVector VecMin(nFeatureCnt, false);
	TVector VecMax(nFeatureCnt, false);
	int i;
	for (i = 0; i < nFeatureCnt; i ++)
		((TImage*) ppFeatures[i])->GetIntensityRange(VecMax(i), VecMin(i), pMask);

	// Parameter estimation
	double dbEnergy, dbMinEnergy = -1;
	int nSampleCnt = 0;  // iteration number and sample count
	do {
		InitGaussMixParam(VecMin, VecMax);  // initialize randomly the parameters
		dbEnergy = GaussMixEMIterations(0.01, 5, ppFeatures, pMask);

		if (dbMinEnergy < 0 ||  // first iteration
		    dbEnergy < dbMinEnergy)  
		{
			dbMinEnergy = dbEnergy;
			for (i = 0; i < nClassCnt; i ++)
			{
				ppMuBest[i]->CopyFrom(m_ppMu[i]);
				ppCovInvBest[i]->CopyFrom(m_ppCovInv[i]);
			}
			memcpy(dbPriBest, m_dbPri, nClassCnt * sizeof(double));
			memcpy(dbDetermBest, m_dbDeterm, nClassCnt * sizeof(double));
		}
		printf("%d", nSampleCnt);
	} while (++nSampleCnt < 3);

	// Set the parameters with the obtained best and continue with the EM
	for (i = 0; i < nClassCnt; i ++)
	{
		m_ppMu[i]->CopyFrom(ppMuBest[i]);
		m_ppCovInv[i]->CopyFrom(ppCovInvBest[i]);
	}
	memcpy(m_dbPri, dbPriBest, nClassCnt * sizeof(double));
	memcpy(m_dbDeterm, dbDetermBest, nClassCnt * sizeof(double));

	// Further iterations to converge
	dbEnergy = GaussMixEMIterations(0.0001, 20, ppFeatures, pMask);

	// Arrange the clusters in the ascending order of the first mean feature 
	// values, so that the resulting map looks more pleasant
	SortCluster();

	//// Show the parameters
	//ShowGaussianMixParam();

	for (i = 0; i < nClassCnt; i ++)
	{
		delete ppMuBest[i];
		delete ppCovInvBest[i];
	}
	delete[] ppMuBest;
	delete[] ppCovInvBest;
	delete[] dbPriBest;
	delete[] dbDetermBest;
}


// 
// The EM for Gaussian mixture estimation in fast mode (i.e. using histogram)
// 
double TEstimateUtil::GaussMixEMIterations(double dbRatio, int nIterations,
                                           int* pHist, int nHistBinCnt)
{
	TVector **ppMuComp = AllocVectors(m_nClassCnt, 1);
	TMatrix **ppCovComp = AllocMatrice(m_nClassCnt, 1);
	double* dbPri = new double[m_nClassCnt];  // for the current compuation
	int iter = 0, i;
	double dbOldEnergy, dbNewEnergy = -1, dbPriTotal;
	TVector VecFeature(1, false);  // feature vector
	TVector VecProb(m_nClassCnt, false);  // probability
	do {
		dbOldEnergy = dbNewEnergy;
		dbNewEnergy = 0;
		SetZeroVector(ppMuComp, m_nClassCnt);
		SetZeroMatrix(ppCovComp, m_nClassCnt);
		memset(dbPri, 0, sizeof(double) * m_nClassCnt);
		for (i = 0; i < nHistBinCnt; i ++)
		{
			// E step
			VecFeature(0) = i;
			dbNewEnergy += pHist[i] * GaussMixExpectation(VecProb, VecFeature);

			// M step
			VecProb *= pHist[i];
			GaussMixMaximization(ppMuComp, ppCovComp, dbPri, VecProb, VecFeature);
		}

		// Compute the total of the priors
		dbPriTotal = 0;
		for (i = 0; i < m_nClassCnt; i ++)
			dbPriTotal += dbPri[i];

		// Update the parameters
		for (i = 0; i < m_nClassCnt; i ++)
		{
			(*m_ppMu[i])(0) = (*ppMuComp[i])(0) / dbPri[i];  // mean
			(*ppCovComp[i])(0, 0) = (*ppCovComp[i])(0, 0) / dbPri[i] -   // covariance
				                    (*m_ppMu[i])(0) * (*m_ppMu[i])(0);
			if ( (*ppCovComp[i])(0, 0) < MIN_VAR )
				(*ppCovComp[i])(0, 0) = MIN_VAR;
			ppCovComp[i]->Inv(m_ppCovInv[i]);  // inverse of covariance
			m_dbPri[i] = dbPri[i] / dbPriTotal;  // prior
			m_dbDeterm[i] = sqrt(ppCovComp[i]->Determinant());  // determinant
		}

		if ( dbOldEnergy < 0 )  // First iteration
			continue;
	} while ( ( iter ++ < nIterations ) && ( fabs(dbOldEnergy - dbNewEnergy) > dbRatio * dbNewEnergy ) );

	for (i = 0; i < m_nClassCnt; i ++)
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
// The EM for Gaussian mixture estimation
// 
double TEstimateUtil::GaussMixEMIterations(double dbRatio, int nIterations,
                                           TGrayImage<FLOAT>** ppFeatures,
                                           TMonoImage* pMask)
{
	TVector **ppMuComp = AllocVectors(m_nClassCnt, m_nFeatureCnt);
	TMatrix **ppCovComp = AllocMatrice(m_nClassCnt, m_nFeatureCnt);
	double* dbPri = new double[m_nClassCnt];  // for the current compuation
	int iter = 0, i, k;
	double dbOldEnergy, dbNewEnergy = -1, dbPriTotal;
	TVector VecFeature(m_nFeatureCnt, false);  // feature vector
	TVector VecProb(m_nClassCnt, false);  // probability
	do {
		dbOldEnergy = dbNewEnergy;
		dbNewEnergy = 0;
		SetZeroVector(ppMuComp, m_nClassCnt);
		SetZeroMatrix(ppCovComp, m_nClassCnt);
		memset(dbPri, 0, sizeof(double) * m_nClassCnt);
		for (int r = ppFeatures[0]->Height() - 1; r >= 0; r --)
			for (int c = ppFeatures[0]->Width() - 1; c >= 0; c --)
			{
				if (pMask && !(*pMask)(r, c))
					continue;

				// E step
				for (i = 0; i < m_nFeatureCnt; i ++)
					VecFeature(i) = (*ppFeatures[i])(r, c);
				dbNewEnergy += GaussMixExpectation(VecProb, VecFeature);

				// M step
				GaussMixMaximization(ppMuComp, ppCovComp, dbPri, VecProb, VecFeature);
			}

		// Compute the total of the priors
		dbPriTotal = 0;
		for (k = 0; k < m_nClassCnt; k ++)
			dbPriTotal += dbPri[k];
		// Update the parameters
		for (k = 0; k < m_nClassCnt; k ++)
		{
			for (i = 0; i < m_nFeatureCnt; i ++)  // mean
				(*m_ppMu[k])(i) = (*ppMuComp[k])(i) / dbPri[k];
			for (i = 0; i < m_nFeatureCnt; i ++)  // covariance
			{
				for (int j = i; j < m_nFeatureCnt; j ++)
					(*ppCovComp[k])(i, j) = (*ppCovComp[k])(i, j) / dbPri[k] - 
					                        (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
			}
			for (i = 0; i < m_nFeatureCnt; i ++)
			{
				if ((*ppCovComp[k])(i, i) < MIN_VAR)
					(*ppCovComp[k])(i, i) = MIN_VAR;
				for (int j = i + 1; j < m_nFeatureCnt; j ++)
					(*ppCovComp[k])(j, i) = (*ppCovComp[k])(i, j);
			}
			ppCovComp[k]->Inv(m_ppCovInv[k]);  // inverse of covariance
			m_dbPri[k] = dbPri[k] / dbPriTotal;  // prior
			m_dbDeterm[k] = sqrt(ppCovComp[k]->Determinant());  // determinant
		}

		if (dbOldEnergy < 0)
			continue;

	} while ( ( iter ++ < nIterations ) && ( fabs(dbOldEnergy - dbNewEnergy) > dbRatio * dbNewEnergy ) );

	for (i = 0; i < m_nClassCnt; i ++)
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
// Initialize the Gaussian mixutre parameters randomly
//
void TEstimateUtil::InitGaussMixParam(TVector& VecMin, TVector& VecMax)
{
	int i, k;
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_ppCovInv[k]->SetZero();
		m_dbPri[k] = 1.0 / m_nClassCnt;  // assume the priors equal
	}

	double dbRange, dbMin;
	for ( i = 0; i < m_nFeatureCnt; i++ )
	{
		dbMin = VecMin(i);
		dbRange = VecMax(i) - dbMin;
		for ( k = 0; k < m_nClassCnt; k++ )
		{
			// Choose means randomly
			(*m_ppMu[k])(i) = dbMin + rand() * dbRange / RAND_MAX;
			// Evenly divide the range to m_nClassCnt parts, and use half as the
			// standard deviation. Different features are assumed indepent first
			(*m_ppCovInv[k])(i, i) = 4 /
				((dbRange / m_nClassCnt) * (dbRange / m_nClassCnt));
		}
	}

	double dbDeterm = 1;
	for ( i = 0; i < m_nFeatureCnt; i++ )
		dbDeterm *= ( 1 / (*m_ppCovInv[0])(i, i) );
	dbDeterm = sqrt(dbDeterm);   // Note it is the square root of the determinant
	for ( k = 0; k < m_nClassCnt; k++ )
		m_dbDeterm[k] = dbDeterm;
}

// 
// The expectation step for Gaussian mixture estimation
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
double TEstimateUtil::GaussMixExpectation(TVector& VecProb, TVector& VecFeature)
{
	double dbTotal = 0;
	int i;
	for (i = 0; i < m_nClassCnt; i ++)
	{
		for (int j = VecFeature.GetElementCnt() - 1; j >= 0; j --)
			(*m_pVecTemp1)(j) = VecFeature(j) - (*m_ppMu[i])(j);
		m_pVecTemp1->Mul(m_pVecTemp2, m_ppCovInv[i]);
		VecProb(i) = exp(- m_pVecTemp2->VecMul(m_pVecTemp1) / 2) / m_dbDeterm[i] *
		             m_dbPri[i];
		dbTotal += VecProb(i);
	}

	//// Normalize
	//if ( dbTotal <= ZERO_THRESH )
	//	dbTotal = ZERO_THRESH;
	for ( i = 0; i < m_nClassCnt; i++ )
		VecProb(i) = VecProb(i) / dbTotal;

	return -log(dbTotal);
}

// 
// The maximization step for Gaussian mixture estimation
// 
// Param in
//  VecProb: the probability
//  VecFeature: the feature
//
// Param out
//  ppMu, ppCov, dbPri: the new sum of mean, covariance and prior
// 
void TEstimateUtil::GaussMixMaximization(TVector** ppMu, TMatrix** ppCov, double *dbPri, 
										 TVector& VecProb, TVector& VecFeature)
{
	double dbTmp;
	for ( int k = m_nClassCnt - 1; k >= 0; k-- )
	{
		for ( int i = 0; i < m_nFeatureCnt; i++ )
		{
			dbTmp = VecFeature(i) * VecProb(k);
			(*ppMu[k])(i) += dbTmp;
			for ( int j = i; j < m_nFeatureCnt; j++ )
				(*ppCov[k])(i, j) += dbTmp * VecFeature(j);
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
TGrayImage<int>* TEstimateUtil::GaussMixClassify(TGrayImage<FLOAT>** ppFeatures, TMonoImage* pMask)
{
	int nWidth = ppFeatures[0]->Width();
	int nHeight = ppFeatures[0]->Height();
	TGrayImage<int>* pMap = new TGrayImage<int>(nWidth, nHeight);

	// Classification
	TVector VecFeature(m_nFeatureCnt, false);  // feature vector
	TVector VecProb(m_nClassCnt, false);  // probability
	double dbEnergy = 0, dbMaxProb;
	int nLabel;
	for (int r = 0; r < nHeight; r ++)
	{
		for (int c = 0; c < nWidth; c ++)
		{
			if (pMask && !(*pMask)(r, c))
			{
				(*pMap)(r, c) = -2;
				continue;
			}

			int i;
			for (i = 0; i < m_nFeatureCnt; i ++)
				VecFeature(i) = (*ppFeatures[i])(r, c);
			for (i = 0; i < m_nClassCnt; i ++)
			{
				for (int j = m_nFeatureCnt - 1; j >= 0; j --)
					(*m_pVecTemp1)(j) = VecFeature(j) - (*m_ppMu[i])(j);
				m_pVecTemp1->Mul(m_pVecTemp2, m_ppCovInv[i]);
				VecProb(i) = exp(- m_pVecTemp2->VecMul(m_pVecTemp1) / 2) / m_dbDeterm[i];
			}

			dbMaxProb = 0;
			for (i = 0; i < m_nClassCnt; i ++)
			{
				if (VecProb(i) > dbMaxProb)
				{
					dbMaxProb = VecProb(i);
					nLabel = i;
				}
			}

			(*pMap)(r, c) = nLabel;
		}
	}

	return pMap;
}

// 
// Estimate the parameters of the Gaussain distributed clusters using obtained
// segmentation map
// 
// Parameter in
//  ppFeatures: the feature vector image
//  pMask: the mask image
//  pMap: the map
//  nClassCnt: number of classes
//  nFeatureCnt: number of features
// 
BOOL TEstimateUtil::GaussMixParamUpdate(TGrayImage<FLOAT>** ppFeatures,
                                        TMonoImage* pMask,
                                        TGrayImage<int>* pMap,
                                        int nClassCnt, int nFeatureCnt)
{
	Alloc(nClassCnt, nFeatureCnt);
	SetZeroVector(m_ppMu, nClassCnt);
	SetZeroMatrix(m_ppCovInv, nClassCnt);
	memset(m_dbPri, 0, sizeof(double) * nClassCnt);

	int i, j, k;
	for (int r = ppFeatures[0]->Height() - 1; r >= 0; r --)
		for (int c = ppFeatures[0]->Width() - 1; c >= 0; c --)
		{
			if (pMask && !(*pMask)(r, c))
				continue;
			if ((k = (*pMap)(r, c)) < 0)
				continue;

			for (i = 0; i < nFeatureCnt; i ++)
				(*m_ppMu[k])(i) += (*ppFeatures[i])(r, c);
			for (i = 0; i < nFeatureCnt; i ++)
				for (j = i; j < nFeatureCnt; j ++)
				{
					(*m_ppCovInv[k])(i, j) += (*ppFeatures[i])(r, c) * (*ppFeatures[j])(r, c);
				}

			m_dbPri[k]++;
		}

	// Compute the total of the priors
	double dbPriTotal = 0;
	for (k = 0; k < nClassCnt; k ++)
		dbPriTotal += m_dbPri[k];
	// Compute the parameters
	for (k = 0; k < nClassCnt; k ++)
	{
		if (m_dbPri[k] > 0)
		{
			(*m_ppMu[k]) *= 1 / m_dbPri[k];  // mean
			for (i = 0; i < nFeatureCnt; i ++)  // covariance
				for (j = i; j < nFeatureCnt; j ++)
				{
					(*m_ppCovInv[k])(i, j) = (*m_ppCovInv[k])(i, j) / m_dbPri[k] - 
											 (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
				}
			for (i = 0; i < nFeatureCnt; i ++)
				for (int j = i + 1; j < nFeatureCnt; j ++)
				{
					(*m_ppCovInv[k])(j, i) = (*m_ppCovInv[k])(i, j);
				}

			m_dbDeterm[k] = sqrt(m_ppCovInv[k]->Determinant());  // determinant
			m_ppCovInv[k]->Inv(m_ppCovInv[k]);  // inverse of covariance
			m_dbPri[k] = m_dbPri[k] / dbPriTotal;  // prior
		}
		else
			return false;
	}

	// Show the parameters
	// ShowGaussianMixParam();
	return true;
}

// 
// Display the Gaussian mixture parameters
// 
void TEstimateUtil::ShowGaussianMixParam()
{
	int i, k;
	for (k = 0; k < m_nClassCnt; k ++)
	{
		cout << "Mu: ";
		for (i = 0; i < m_nFeatureCnt; i ++)
			cout << (*m_ppMu[k])(i) << ", ";
		cout << "Cov: ";
		for (i = 0; i < m_nFeatureCnt; i ++)
			for (int j = 0; j < m_nFeatureCnt; j ++)
			{
				cout << (*m_ppCovInv[k])(i, j) << ", ";
			}
		cout << "Pri: " << m_dbPri[k] << endl;
	}
}

// 
// Generate a random map
//
// Param in
//  nClusterCnt: the number of clusters
// 
// Param out
//  pMap: generated map
// 
void TEstimateUtil::GenerateRandomMap(TGrayImage<int>* pMap, TMonoImage* pMask, 
                                      int nClusterCnt)
{
	for (int i = pMap->Height() - 1; i >= 0; i --)
		for (int j = pMap->Width() - 1; j >= 0; j --)
	{
		if (pMask && !(*pMask)(i, j))
			continue;
		(*pMap)(i, j) = rand() * nClusterCnt / (RAND_MAX + 1.0);
	}
}

// 
// Get the maximum lable of the map
// 
int TEstimateUtil::GetMaxLabel(TGrayImage<int>* pMap, TMonoImage* pMask)
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int nMax = -1;
	for (int r = 0; r < nHeight; r ++)
		for (int c = 0; c < nWidth; c ++)
	{
		if (pMask && !(*pMask)(r, c))
			continue;

		if ((*pMap)(r, c) > nMax)
			nMax = (*pMap)(r, c);
	}

	return nMax;
}


// 
// Get the number of positive labels of the map
// 
int TEstimateUtil::GetNumPositiveLabels(TGrayImage<int>* pMap, TMonoImage* pMask)
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();

	// max label in image.
	int nMax = GetMaxLabel(pMap, pMask) + 1;  // because the label begins with 0

	// indicates whether each gray level exists.
	char* cCnts = new char[nMax];
	memset(cCnts, 0, nMax*sizeof(char));

	for (int r = 0; r < nHeight; r ++)
		for (int c = 0; c < nWidth; c ++)
	{
		if (pMask && !(*pMask)(r, c))
			continue;

		if ((*pMap)(r, c) > 0)
			cCnts[(*pMap)(r, c)] = 1;
	}

	int cnt = 0;

	for (int i = 0; i < nMax; i++)
	{
		if (cCnts[i] == 1)
			cnt++;
	}

	delete [] cCnts;

	return cnt;
}

//
// Generate the consecutive mapping
// 
// Param in
//  pMap: the segmentation map
//  pMask: the mask
// 
// Param out
//  nMappingItemCnt: the number of mapping items
// 
// Return
//  the mapping lookup table
//
int* TEstimateUtil::GenerateConsecutiveMapping(int& nMappingItemCnt,
                                               TGrayImage<int>* pMap, TMonoImage* pMask)
{
	// Find the maximum label of the map
	nMappingItemCnt = GetMaxLabel(pMap, pMask) + 1;  // because the label begins with 0

	if(nMappingItemCnt<1)
		return 0;

	// Generate the mapping from non-consecutive labels to consecutive labels
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int* nMapping = new int[nMappingItemCnt];
	memset(nMapping, 0, nMappingItemCnt * sizeof(int));
	for ( int r = 0; r < nHeight; r++ )
		for ( int c = 0; c < nWidth; c++ )
		{
			if (pMask && !(*pMask)(r, c))
				continue;

			// Mark the labels (excluding the boundary label) that show up in the map
			if ( (*pMap)(r, c) >= 0 )
				nMapping[(*pMap)(r, c)] = 1;
		}

	int i;
	for ( i = 1; i < nMappingItemCnt; i++ )
		nMapping[i] += nMapping[i - 1];
	for ( i = 0; i < nMappingItemCnt; i++ )
		nMapping[i]--;  // labels should begin with 0

	return nMapping;
}

int* TEstimateUtil::GenerateMergingLUT(int& nMappingItemCnt, TGrayImage<int>* pMapOld, TGrayImage<int>* pMapNew, TMonoImage* pMask)
{
    // Find the maximum label of the map
    nMappingItemCnt = GetMaxLabel(pMapOld, pMask) + 1;  // because the label begins with 0

    if(nMappingItemCnt<1)
        return 0;

    // Generate the mapping from non-consecutive labels to consecutive labels
    int nWidth = pMapOld->Width();
    int nHeight = pMapOld->Height();
    int* nMapping = new int[nMappingItemCnt];
    memset(nMapping, 0, nMappingItemCnt * sizeof(int));
    for ( int r = 0; r < nHeight; r++ )
        for ( int c = 0; c < nWidth; c++ )
        {
            if (pMask && !(*pMask)(r, c))
                continue;

            // Old labels become new labels
            if ( (*pMapOld)(r, c) >= 0 )
                nMapping[(*pMapOld)(r, c)] = (*pMapNew)(r,c);
        }
    return nMapping;
}

//
// Re-labelize the map so that the labels are consecutive
// 
// Return
//  the number of segments in the map
// 
int TEstimateUtil::RelabelMap(TGrayImage<int>* pMap, TMonoImage* pMask)
{
	int nMaxCnt;
	int* nLabels = GenerateConsecutiveMapping(nMaxCnt, pMap, pMask);
	if(nLabels==0)
		return 0;

	// Update the labels in the map
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	for ( int r = 0; r < nHeight; r++ )
		for ( int c = 0; c < nWidth; c++ )
		{
			if (pMask && !(*pMask)(r, c))
				continue;

			if ((*pMap)(r, c) >= 0)
				// Not boundary label
				(*pMap)(r, c) = nLabels[(*pMap)(r, c)];
		}

	int nRet = nLabels[nMaxCnt - 1] + 1;
	//Debug::WriteLine("nRet = " + nRet);
	delete[] nLabels;
	return nRet;
}

TGrayImage<int>* TEstimateUtil :: DispMap(TGrayImage<int>* pMap, TMonoImage* pMask, bool bAddBlackBox)
{
	int nHeight = pMap->Height();
	int nWidth = pMap->Width();

	TGrayImage<int>* pDispMap = new TGrayImage<int>(nWidth, nHeight);

	int r, c, l;
	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;

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

	if ( bAddBlackBox )
	{
		TGrayImage<int>* pDispMap1 = new TGrayImage<int>(nWidth + 2, nHeight + 2);
		for ( r = 1; r < nHeight + 2; r++ )
			for ( c = 1; c < nWidth + 2; c++ )
			{
				if ( r == 0 || c == 0 || r == nHeight + 1 || c == nWidth + 1 )
					(*pDispMap1)(r, c) = 0;
				else
					(*pDispMap1)(r, c) = (*pDispMap)(r-1, c-1);
			}

		delete pDispMap;

		pDispMap = pDispMap1;
	}

	return pDispMap;
}

// Sort the map according to the average feature value in each region
void TEstimateUtil::SortMap(TGrayImage<int>*& pMap, TMonoImage* pMask, TGrayImage<FLOAT>** ppFeatImg, int nFeatureCnt)
{
	int nRCnt = RelabelMap(pMap, pMask);

	int nHeight = pMap->Height();
	int nWidth = pMap->Width();

	int* pRPixelCnt = new int[nRCnt];
	double* pVal = new double[nRCnt];
	memset(pRPixelCnt, 0, nRCnt*sizeof(int));
	memset(pVal, 0, nRCnt*sizeof(double));

	int r, c, l;
	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;

			if ( (*pMap)(r, c) < 0 )
				continue;

			pVal[(*pMap)(r, c)] += (*ppFeatImg[0])(r, c);
			pRPixelCnt[(*pMap)(r, c)]++;
		}
	
	for ( l = 0; l < nRCnt; l++ )
		pVal[l] /= pRPixelCnt[l];

	double* pValSort = new double[nRCnt];
	sort_(pValSort, pRPixelCnt, pVal, nRCnt);
	delete[] pValSort;
	delete[] pVal;

	int* pIdxMapping = new int[nRCnt];
	for ( r = 0; r < nRCnt; r++ )
	{
		for ( c = 0; c < nRCnt; c++ )
			if ( pRPixelCnt[c] == r )
				break;
		pIdxMapping[r] = c;
	}
	delete[] pRPixelCnt;

	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;

			if ( (*pMap)(r, c) < 0 )
				continue;

			(*pMap)(r, c) = pIdxMapping[(*pMap)(r, c)];
		}

	delete[] pIdxMapping;
}

// Add black box to a result segmentation map
TGrayImage<int>* TEstimateUtil::AddBlackBox(TGrayImage<int>* pMap)
{
	int nHeight = pMap->Height();
	int nWidth = pMap->Width();

	int r, c;
	TGrayImage<int>* pMap1 = new TGrayImage<int>(nWidth + 2, nHeight + 2);
	for ( r = 1; r < nHeight + 2; r++ )
		for ( c = 1; c < nWidth + 2; c++ )
		{
			if ( r == 0 || c == 0 || r == nHeight + 1 || c == nWidth + 1 )
				(*pMap1)(r, c) = 0;
			else
				(*pMap1)(r, c) = (*pMap)(r-1, c-1);
		}

	return pMap1;
}

// 
// Erase the boundary
// 
// Param in
//  pMap: the obtained segmenation map
//  ppFeatures: the feature images
//  pMask: the mask
//  nFeatureCnt: number of features
// 
// Param out
//  pMap: the map with the boundary erased
// 
void TEstimateUtil::EraseBoundary(TGrayImage<int>* pMap,
                                  TGrayImage<FLOAT>** ppFeatImg, int nFeatureCnt,
                                  TMonoImage* pMask)
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();

	// Relabel to make them consecutive
	int nMaxCnt;
	int* nForwardMapping = GenerateConsecutiveMapping(nMaxCnt, pMap, pMask);
	int r, c;
	for (r = pMap->Height() - 1; r >= 0; r --)
		for (c = pMap->Width() - 1; c >= 0; c --)
	{
		if (pMask && !(*pMask)(r, c))
			continue;

		if ((*pMap)(r, c) >= 0)
			// Not boundary label
			(*pMap)(r, c) = nForwardMapping[(*pMap)(r, c)];
	}
	int nClassCnt = nForwardMapping[nMaxCnt - 1] + 1;

	// Generate the reverse mapping
	int* nReverseMapping = new int[nClassCnt];
	int k, l;
	for (k = 0; k < nClassCnt; k ++)
		for (l = 0; l < nMaxCnt; l ++)
	{
		if (nForwardMapping[l] == k)
		{
			nReverseMapping[k] = l;
			break;
		}
	}

	// Erase boundary sites that are within the interior of the same class
	int nRow, nCol;
	int nLabel;
	for (r = 0; r < nHeight; r ++)
		for (c = 0; c < nWidth; c ++)
	{
		if (pMask && !(*pMask)(r, c))
			continue;
		if ((*pMap)(r, c) < 0)
		{
			// Check if there are two different non-boundary labels in the
			// neighborhood
			nLabel = -1;
			int j;
			for (j = 0; j < 8; j ++)
			{
				nRow = r + neighbor8[j].Row;
				nCol = c + neighbor8[j].Col;
				if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
					continue;
				if (pMask && !(*pMask)(nRow, nCol))  // not valid
					continue;
				if ((l = (*pMap)(nRow, nCol)) < 0)  // also boundary
					continue;
				if (l != nLabel && nLabel >= 0)
					break;
				nLabel = l;
			}

			if (j >= 8 && nLabel >= 0)  // inside the same class
				(*pMap)(r, c) = nLabel;
		}
	}

	// Erase the border
	TVector VecFeature(nFeatureCnt, false);  // feature vector
	TVector VecProb(nClassCnt, false);  // probability
	double dbMinE;
	for (k = 0; k < 1; k ++)  // several iterations to refine
	{
		// Compute the cluster parameters
		GaussMixParamUpdate(ppFeatImg, pMask, pMap, nClassCnt, nFeatureCnt);
		for (int i = 0; i < nClassCnt; i ++)
			m_dbPri[i] = 1.0 / nClassCnt;  // set the priors for all classes to be equal, since we use spatial context as prior

		for (r = 0; r < nHeight; r ++)
			for (c = 0; c < nWidth; c ++)
		{
			if (pMask && !(*pMask)(r, c))
				continue;
			if (IsBoundary(pMap, pMask, r, c))
			{
				// Compute the probabilities
				int i;
				for (i = 0; i < nFeatureCnt; i ++)
					VecFeature(i) = (*ppFeatImg[i])(r, c);
				GaussMixExpectation(VecProb, VecFeature);

				// Transform to energy
				for (i = 0; i < nClassCnt; i ++)
					VecProb(i) = -log(VecProb(i));

				// Account spatial context constraint
				for (int j = 0; j < 8; j ++)
				{
					nRow = r + neighbor8[j].Row;
					nCol = c + neighbor8[j].Col;
					if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
						continue;
					if (pMask && !(*pMask)(nRow, nCol))  // not valid
						continue;
					if ((l = (*pMap)(nRow, nCol)) < 0)  // also boundary
						continue;
					VecProb(l) -= BETA;
				}

				// Determine the label
				nLabel = -1;
				dbMinE = 100000000;  // a large enough number
				for (i = 0; i < nClassCnt; i ++)
				{
					if (VecProb(i) < dbMinE)
					{
						dbMinE = VecProb(i);
						nLabel = i;
					}
				}

				(*pMap)(r, c) = nLabel;
			}
		}
	}

	for (r = pMap->Height() - 1; r >= 0; r --)
		for (c = pMap->Width() - 1; c >= 0; c --)
	{
		if (pMask && !(*pMask)(r, c))
			continue;
		(*pMap)(r, c) = nReverseMapping[(*pMap)(r, c)];
	}

	delete nForwardMapping;
	delete nReverseMapping;
}

//
// Check if the current site is at boundary
// 
BOOL TEstimateUtil::IsBoundary(TGrayImage<int>* pMap, TMonoImage* pMask, int nRow, int nCol)
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int r, c;
	int l = (*pMap)(nRow, nCol);
	for (int i = 0; i < 8; i ++)
	{
		r = nRow + neighbor8[i].Row;
		c = nCol + neighbor8[i].Col;
		if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
			continue;
		if (pMask && !(*pMask)(r, c))
			continue;

		if ((*pMap)(r, c) != l)
			return true;
	}

	return false;
}

// 
// Estimate the noise variance (i.e. noise variance assuming constant multi-dimensional image features corrupted by noise) //QIN_MOD
//
double TEstimateUtil::EstVar(TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt, TMonoImage* pMask) //QIN_MOD
{
	int nWin = (int) sqrt( (double) (nFeatureCnt * 10) ); // window size (odd number) for computing the local variance
	if ( nWin % 2 == 0 )
		nWin++;
	if ( nWin < 7 )
		nWin = 7;

	int nSrcWidth = ppSrcImg[0]->Width(); //assuming each component image has the same size
	int nSrcHeight = ppSrcImg[0]->Height();

	if (nWin > nSrcWidth)
		nWin = nSrcWidth;
	else if (nWin > nSrcHeight)
		nWin = nSrcHeight;

	int nSize = nWin * nWin;

	int nDstWidth = nSrcWidth - (nWin - 1);
	int nDstHeight = nSrcHeight - (nWin - 1);
	TGrayImage<FLOAT>* pVar = new TGrayImage<FLOAT>(nDstWidth, nDstHeight);
	
	if (nFeatureCnt == 1)
	{
	  	// Compute the local mean and local square mean
    	// First, Average in horizontal direction
    	TGrayImage<FLOAT>* pAvgH = new TGrayImage<FLOAT>(nDstWidth, nSrcHeight);
    	TGrayImage<FLOAT>* pVarH = new TGrayImage<FLOAT>(nDstWidth, nSrcHeight);
	    float fVal;
	    int r, c;
	    for (r = 0; r < nSrcHeight; r ++)
	    {
	    	(*pAvgH)(r, 0) = 0;
	    	(*pVarH)(r, 0) = 0;
	    	for (c = 0; c < nWin; c ++)
	    	{
	    		fVal = (*ppSrcImg[0])(r, c);
		    	(*pAvgH)(r, 0) += fVal;
		    	(*pVarH)(r, 0) += fVal * fVal;
	    	}

	    	for (c = 1; c < nDstWidth; c ++)
	    	{
		    	fVal = (*ppSrcImg[0])(r, c - 1);
		    	(*pAvgH)(r, c) = (*pAvgH)(r, c - 1) - fVal;
		    	(*pVarH)(r, c) = (*pVarH)(r, c - 1) - fVal * fVal;

		    	fVal = (*ppSrcImg[0])(r, c + nWin - 1);
	    		(*pAvgH)(r, c) += fVal;
	    		(*pVarH)(r, c) += fVal * fVal;
	    	}
    	}
    	// Average in vertical direction
    	TGrayImage<FLOAT>* pAvg = new TGrayImage<FLOAT>(nDstWidth, nDstHeight);
    	for (c = 0; c < nDstWidth; c ++)
	    {
	    	(*pAvg)(0, c) = 0;
	    	(*pVar)(0, c) = 0;
	    	for (r = 0; r < nWin; r ++)
	    	{
	    		(*pAvg)(0, c) += (*pAvgH)(r, c);
	    		(*pVar)(0, c) += (*pVarH)(r, c);
	    	}

	    	for (r = 1; r < nDstHeight; r ++)
	    	{
	    		(*pAvg)(r, c) = (*pAvg)(r - 1, c) - (*pAvgH)(r - 1, c);
	    		(*pVar)(r, c) = (*pVar)(r - 1, c) - (*pVarH)(r - 1, c);

	    		(*pAvg)(r, c) += (*pAvgH)(r + nWin - 1, c);
	    		(*pVar)(r, c) += (*pVarH)(r + nWin - 1, c);
	    	}
	    }
	    delete pAvgH;
	    delete pVarH;

	    for (r = 0; r < nDstHeight; r ++)
	    	for (c = 0; c < nDstWidth; c ++)
	    {
	    	fVal = (*pAvg)(r, c) / nSize;
	    	if ( abs(fVal) <= ZERO_THRESH )
		    	(*pVar)(r, c) = 0;
		    else
		    {
		    	fVal *= fVal;
		    	(*pVar)(r, c) = (*pVar)(r, c) / nSize - fVal;  // variance for additive noise
    //			(*pVar)(r, c) = ((*pVar)(r, c) / nSize - fVal) / fVal;  // variance for multiplicative noise
	    	}
	    }
	    delete pAvg;
	}
	//else if (nFeatureCnt > 1)
	//{
 //       // Compute the local mean and the determinant of the local covariance matrix
	//    int r, c, k, m, s, t;
	//    TMatrix* CovMat = new TMatrix(nFeatureCnt, nFeatureCnt);
	//    double* dbMean = new double[nFeatureCnt];
	//    for (r = 0; r < nDstHeight; r++)
	//	{
	//		printf("%d", r);
	//	    for (c = 0; c < nDstWidth; c++)
	//	    {
	//		    CovMat -> SetZero();
	//		    memset(dbMean, 0, sizeof(double)*nFeatureCnt);
	//		    for (k = r; k < r + nWin; k++)
	//		    {
	//			    for (m = c; m < c + nWin; m++)
	//				    for (s = 0; s < nFeatureCnt; s++)
	//				    {
	//				    	dbMean[s] += (*ppSrcImg[s])(k, m);
	//						(*CovMat)(s, s) += (*ppSrcImg[s])(k, m)*(*ppSrcImg[s])(k, m);
	//				    	for (t = s + 1; t < nFeatureCnt; t++)
	//						{
	//						    (*CovMat)(s, t) += (*ppSrcImg[s])(k, m)*(*ppSrcImg[t])(k, m);
 //                               (*CovMat)(t, s) = (*CovMat)(s, t);
	//						}
	//				    }
	//		    }

	//	    	for (s = 0; s < nFeatureCnt; s++)
	//		    {
	//		    	(*CovMat)(s, s) = (*CovMat)(s, s)/nSize - (dbMean[s]/nSize)*(dbMean[s]/nSize);
	//		    	for (t = s + 1; t < nFeatureCnt; t++)
	//		    	{
	//		    		(*CovMat)(s, t) = (*CovMat)(s, t)/nSize - (dbMean[s]/nSize)*(dbMean[t]/nSize);
	//		    	    (*CovMat)(t, s) = (*CovMat)(s, t);
	//		    	}
	//		    }

	//		    (*pVar)(r, c) = (float) CovMat -> Determinant();
	//	    }
	//	}

	//    delete CovMat;
	//    delete[] dbMean;
	//}

	else if (nFeatureCnt > 1)
	{
		//double dbMax = -999999999999;
        // Compute the local mean and the determinant of the local covariance matrix
	    int r, c, k, m, s, t;
	    TMatrix* CovMat = new TMatrix(nFeatureCnt, nFeatureCnt);
		double** ppCovMatData;
	    double* pMean = new double[nFeatureCnt];
		float tmp1, tmp2;
	    for (r = 0; r < nDstHeight; r++)
		{
		    for (c = 0; c < nDstWidth; c++)
		    {
				if ( c == 110 )
					c = 110;
			    CovMat -> SetZero();
				ppCovMatData = CovMat->m_ppData;
			    memset(pMean, 0, sizeof(double)*nFeatureCnt);
			    for (k = r; k < r + nWin; k++)
			    {
				    for (m = c; m < c + nWin; m++)
					    for (s = 0; s < nFeatureCnt; s++)
					    {
							tmp1 = (*ppSrcImg[s])(k, m);
					    	pMean[s] += tmp1;
							ppCovMatData[s][s] += tmp1*tmp1;
					    	for (t = s + 1; t < nFeatureCnt; t++)
							{
								tmp2 = (*ppSrcImg[t])(k, m);
							    ppCovMatData[s][t] += tmp1*tmp2;
                                ppCovMatData[t][s] = ppCovMatData[s][t];
							}
					    }
			    }

		    	for (s = 0; s < nFeatureCnt; s++)
			    {
			    	ppCovMatData[s][s] = ppCovMatData[s][s]/nSize - (pMean[s]/nSize)*(pMean[s]/nSize);
			    	for (t = s + 1; t < nFeatureCnt; t++)
			    	{
			    		ppCovMatData[s][t] = ppCovMatData[s][t]/nSize - (pMean[s]/nSize)*(pMean[t]/nSize);
			    	    ppCovMatData[t][s] = ppCovMatData[s][t];
			    	}
			    }

			    (*pVar)(r, c) = (float) CovMat -> Determinant();

				//if ( (*pVar)(r, c) > dbMax )
				//	dbMax = (*pVar)(r, c);
			}
			//printf("%d", r);
		}

	    delete CovMat;
	    delete[] pMean;
	}
	else
	{ 
		cout << "Error: The number of features should be above 0!" << endl;
		exit(1);
	}

	double dbNoiseVar;
	if (pMask)
	{
		TMonoImage* pNewMask = new TMonoImage(pMask);
		pNewMask->Depad(nWin / 2, nWin / 2, nWin / 2, nWin / 2);
		dbNoiseVar = FindMode(pVar, pNewMask);
		delete pNewMask;
	}
	else
		dbNoiseVar = FindMode(pVar, 0);
	delete pVar;
	return dbNoiseVar;
	//return 0;
}

// 
// Find the mode of the distribution
// 
double TEstimateUtil::FindMode(TGrayImage<FLOAT>* pImg, TMonoImage* pMask)
{
	double dbMin, dbMax;
	((TImage*) pImg)->GetIntensityRange(dbMax, dbMin, pMask);
	int* pHist = pImg->GetIntensityHisto(dbMin, dbMax, 256, pMask);
	int nMaxCnt = 0, nIndex = 0;
	for (int i = 0; i < 256; i ++)
	{
		if (nMaxCnt < pHist[i])
		{
			nMaxCnt = pHist[i];
			nIndex = i;
		}
	}
	delete[] pHist;
	return dbMin + (dbMax - dbMin) / 256 * (nIndex + 0.5);
}

void sort_(double *Y, int *I, double *A, int length)
{
    int i, j;
    double max, *tmp;

    tmp = (double *) calloc(length, sizeof(double));

    for (i=0;i<length;i++)
        tmp[i] = A[i];

    max = tmp[0];
    for (i=1;i<length;i++) {
        if (tmp[i] > max)
            max = tmp[i];
    }

    max = fabs(10*max);

    for (i=0;i<length;i++) {
        Y[i] = tmp[0];
        I[i] = 0;
        for (j=1;j<length;j++) {
            if (tmp[j] < Y[i]) {
                Y[i] = tmp[j];
                I[i] = j;
            }
        }

        tmp[I[i]] = max;
    }

    free(tmp);
}