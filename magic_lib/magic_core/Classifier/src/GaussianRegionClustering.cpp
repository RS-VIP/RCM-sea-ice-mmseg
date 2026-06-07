// GaussianRegionClustering.cpp : implementation of the TGaussianRegionClustering class
//
// Clustering the regions by using the MRF on GIEP graph using Gaussian Statistics
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "GaussianRegionClustering.h"
#include "StrengthPenaltyBoundary.h"
#include "EstimateUtil.h"
#include "Matrix.h"  // For n-D minimum Fisher criterion
#include "RegionBasedClustering.h"
#include "PixelBasedClustering.h"
#include "GIEPGraph.h"

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
TGaussianRegionClustering::TGaussianRegionClustering(int nType, int nInitMode, int nInitInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal,
									 int nMode, double dbInitT, int nIterations, int nClassCnt, bool bKeepBoundaries)
	: TRegionClustering(nType, nInitMode, nInitInitMode, dbRatio, nInit1, 
				nInit2, nIterFinal, nMode, dbInitT, 
				nIterations, nClassCnt)
{
	m_ppMu = 0;
	m_ppCovInv = 0;

	m_dbPrior = 0;
	m_dbCovDet = 0;

	m_dbLnSigma = 0;
	m_dbLnPrior = 0;
	m_bKeepBoundaries = bKeepBoundaries;
}

//
// Destructor
// 
TGaussianRegionClustering::~TGaussianRegionClustering()
{
	if (m_ppMu)
	{
		for ( int i = 0; i < m_nClassCnt; i++ )
			delete m_ppMu[i];
		delete[] m_ppMu;
	}
	if (m_ppCovInv)
	{
		for ( int i = 0; i < m_nClassCnt; i++ )
			delete m_ppCovInv[i];
		delete[] m_ppCovInv;
	}
	if (m_dbPrior) delete[] m_dbPrior; 
	if (m_dbCovDet) delete[] m_dbCovDet;
	if (m_dbLnSigma) delete[] m_dbLnSigma;  //QIN_MOD
	if (m_dbLnPrior) delete[] m_dbLnPrior;  //QIN_MOD
}

// set the cluster means
bool TGaussianRegionClustering::SetMu(TVector** ppClusterMeans)
{
    int nFeatureCnt = m_ppMu[0]->GetElementCnt();
    for(int i=0; i<nFeatureCnt; i++)
        m_ppMu[i] = ppClusterMeans[i];
    return true;
}

// Set the inverse covariance
bool TGaussianRegionClustering::SetCov(TMatrix** ppClusterCovInv)
{
    int nFeatureCnt = m_ppMu[0]->GetElementCnt();
    for(int i=0; i<nFeatureCnt; i++)
        m_ppCovInv[i] = ppClusterCovInv[i];
    return true;
}

// 
// Get the energy corresponding to the unary property
// 
// Param in
//  pp: the region being classified
//  k: the current class investigated
// 
double TGaussianRegionClustering::GetUnaryEnergy(TRegion* pp, int k)
{
	TGaussianRegion* p = (TGaussianRegion*) pp;
	int nFeatureCnt = m_ppMu[k]->GetElementCnt();
	double *dbMu = m_ppMu[k]->m_ppData[0];
	double **dbCoMu = m_ppCovInv[k]->m_ppData;
	double U = 0;
	for (int i = 0; i < nFeatureCnt; i ++)
	{
		U += ((*p->m_pCoMean)(i, i) - 2 * (*p->m_pMean)(i) * dbMu[i] +
		      dbMu[i] * dbMu[i]) * dbCoMu[i][i];
		for (int j = i + 1; j < nFeatureCnt; j ++)
			U += ((*p->m_pCoMean)(i, j) - (*p->m_pMean)(i) * dbMu[j] -
			      (*p->m_pMean)(j) * dbMu[i] + dbMu[i] * dbMu[j]) *
			     dbCoMu[i][j] * 2;
	}
	U = U / 2 + m_dbLnSigma[k];  //Note that the determinant of COV is squared in the program
	if ( m_nType & TRegionClustering::VAFM )
		return ( U * p->m_nPixelCnt / nFeatureCnt * ((TGIEPGraph*) m_pGraph)->GetCurAlpha() ); //QIN_nFeatureCnt
	else
		return ( U * p->m_nPixelCnt / nFeatureCnt ); //QIN_nFeatureCnt
}

// 
// Get the energy corresponding to the binary relation
// 
// Param in
//  pRegion: the region being classified
//  pBoundary: the boundary connecting the two regions of interest
//  k: the current class investigated
// 
double TGaussianRegionClustering::GetBinaryEnergy(TRegion* pRegion, TBoundary* pBoundary, int k)
{
	TRegion* pNeighbor = (TRegion*) pRegion->GetNeighborVertex(pBoundary);
	if (pNeighbor->m_nClassLabel != k)
	{
		if ( m_nType & TRegionClustering::GIEP )
			return ((TStrengthPenaltyBoundary*) pBoundary)->m_dbPenaltyEnergy; //QIN_nFeatureCnt
		else
			return ( pBoundary->m_nLength * ((TGIEPGraph*) m_pGraph)->GetCurBeta() ); //QIN_nFeatureCnt
	}
	else
		return 0;
}

// 
// Initialize the classifier
// TODO: there's a lot of repeated code here (the if clauses). Clean it up.
bool TGaussianRegionClustering::Init()
{
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	if ( m_nClassCnt <= 0 || nFeatureCnt <= 0 )
		return false;

	m_ppMu = new TVector*[m_nClassCnt];
	m_ppCovInv = new TMatrix*[m_nClassCnt];
	for ( int i = 0; i < m_nClassCnt; i++ )
	{
		m_ppMu[i] = new TVector(nFeatureCnt, false);
		m_ppMu[i]->SetZero();
		m_ppCovInv[i] = new TMatrix(nFeatureCnt, nFeatureCnt);
		m_ppCovInv[i]->SetZero();
	}

	m_dbPrior = new double[m_nClassCnt];
	m_dbCovDet = new double[m_nClassCnt];

	m_dbLnSigma = new double[m_nClassCnt];
	m_dbLnPrior = new double[m_nClassCnt];

	
	TMarkovNetClassifier::Init();

	if ( m_nInitMode == RAND1 )
	{
		int i;
		for ( i = 0; i < m_nClassCnt; i++ )
			m_dbPrior[i] = 1.0 / m_nClassCnt;  // assume the priors equal

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

		int j, k, nCnt = 0;
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
		}
		
		double dbRange, dbMin;
		for ( i = 0; i < nFeatureCnt; i++ )
		{
			dbMin = VecMin[i];
			dbRange = VecMax[i] - dbMin;
			for ( j = 0; j < m_nClassCnt; j++ )
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

		double dbCovDet = 1;
		for ( i = 0; i < nFeatureCnt; i++ )
			dbCovDet *= ( 1 / (*m_ppCovInv[0])(i, i) );
		dbCovDet = sqrt(dbCovDet);   // Note it is the square root of the determinant
		for ( i = 0; i < m_nClassCnt; i++ )
			m_dbCovDet[i] = dbCovDet;

		delete[] VecMin;
		delete[] VecMax;
		delete[] Mu;
		delete[] Variance;

		for ( i = 0; i < m_nClassCnt; i++ )
		{
			m_dbLnPrior[i] = log(m_dbPrior[i]);
		    m_dbLnSigma[i] = log(m_dbCovDet[i]);
		}

		m_bInitSucFlag = FinalClassify();
	}

	if ( m_nInitMode == RAND2 )
	{
		// Find the range and variance of each feature
		int nFeatureCnt = m_pGraph->GetFeatureCnt();
		TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
		TMonoImage* pMask = m_pGraph->GetMask();

		double* VecMin = new double[nFeatureCnt];
		double* VecMax = new double[nFeatureCnt];
		double* Mu = new double[nFeatureCnt];
		double* Variance = new double[nFeatureCnt];
		int i;
		for ( i = 0; i < nFeatureCnt; i++ )
		{
			((TImage*) ppFeatures[i])->GetIntensityRange(VecMax[i], VecMin[i], pMask);
			Mu[i] = ppFeatures[i]->GetIntensityMean(pMask);
		}

		int j, k, nCnt = 0;
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
		}
		TPixelBasedClustering* pInitPixelClustering = new TPixelBasedClustering(pMask, ppFeatures, m_nClassCnt, nFeatureCnt);
		pInitPixelClustering->PixelInitParam(m_nInitMode, VecMin, VecMax, Variance);
		TVector** ppMu = pInitPixelClustering->GetMean();
		TMatrix** ppCovInv = pInitPixelClustering->GetCovInv();
		double* dbPrior = pInitPixelClustering->GetPrior();
		double* dbCovDet = pInitPixelClustering->GetCovDet();

		for ( i = 0; i < m_nClassCnt; i++ )
		{
			m_ppMu[i]->CopyFrom(ppMu[i]);
			m_ppCovInv[i]->CopyFrom(ppCovInv[i]);
			m_dbPrior[i] = dbPrior[i];
			m_dbCovDet[i] = dbCovDet[i];
			m_dbLnPrior[i] = log(m_dbPrior[i]);
			m_dbLnSigma[i] = log(m_dbCovDet[i]);
		}

		delete[] VecMin;
		delete[] VecMax;
		delete[] Mu;
		delete[] Variance;
		delete pInitPixelClustering;

		m_bInitSucFlag = FinalClassify();   //bug: 2009.04.09
	}



	if ( m_nInitMode == KM_PIXEL )
	{		
		TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
		TMonoImage* pMask = m_pGraph->GetMask();
		TPixelBasedClustering* pInitPixelClustering = new TPixelBasedClustering(pMask, ppFeatures, m_nClassCnt, nFeatureCnt);

		TGrayImage<int>* pRet = 0;

		pRet = pInitPixelClustering->KMClustering(m_nInitInitMode, m_bInitSucFlag, TPixelBasedClustering::FAST, m_dbInitRatio, m_nInit1, m_nInit2, m_nInitIterFinal);

		TVector** ppMu = pInitPixelClustering->GetMean();
		TMatrix** ppCovInv = pInitPixelClustering->GetCovInv();
		double* dbPrior = pInitPixelClustering->GetPrior();
		double* dbCovDet = pInitPixelClustering->GetCovDet();

		for ( int i = 0; i < m_nClassCnt; i++ )
		{
			m_ppMu[i]->CopyFrom(ppMu[i]);
			m_ppCovInv[i]->CopyFrom(ppCovInv[i]);
			m_dbPrior[i] = dbPrior[i];
			m_dbCovDet[i] = dbCovDet[i];
			m_dbLnPrior[i] = log(m_dbPrior[i]);
			m_dbLnSigma[i] = log(m_dbCovDet[i]);
		}

		if (pRet)
			delete pRet; // this HAS to be deleted.

		delete pInitPixelClustering;
            
		m_bInitSucFlag = FinalClassify();
	}

	if ( m_nInitMode == GMM_PIXEL )
	{		
		TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
		TMonoImage* pMask = m_pGraph->GetMask();

		TPixelBasedClustering* pInitPixelClustering = new TPixelBasedClustering(pMask, ppFeatures, m_nClassCnt, nFeatureCnt);

		TGrayImage<int>* pRet = 0;

		pRet = pInitPixelClustering->GaussMixClustering(m_nInitInitMode, m_bInitSucFlag, TPixelBasedClustering::FAST, m_dbInitRatio, m_nInit1, m_nInit2, m_nInitIterFinal);

		TVector** ppMu = pInitPixelClustering->GetMean();
		TMatrix** ppCovInv = pInitPixelClustering->GetCovInv();
		double* dbPrior = pInitPixelClustering->GetPrior();
		double* dbCovDet = pInitPixelClustering->GetCovDet();

		for ( int i = 0; i < m_nClassCnt; i++ )
		{
			m_ppMu[i]->CopyFrom(ppMu[i]);
			m_ppCovInv[i]->CopyFrom(ppCovInv[i]);
			m_dbPrior[i] = dbPrior[i];
			m_dbCovDet[i] = dbCovDet[i];
			m_dbLnPrior[i] = log(m_dbPrior[i]);
			m_dbLnSigma[i] = log(m_dbCovDet[i]);
		}

		if (pRet)
			delete pRet; // this HAS to be deleted.

		delete pInitPixelClustering;
            
		m_bInitSucFlag = FinalClassify();

	}

	if ( m_nInitMode == KMGMM_PIXEL )
	{
		TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
		TMonoImage* pMask = m_pGraph->GetMask();

		TPixelBasedClustering* pInitPixelClustering = new TPixelBasedClustering(pMask, ppFeatures, m_nClassCnt, nFeatureCnt);

		TGrayImage<int>* pRet = 0;

		pRet = pInitPixelClustering->KMGMMClustering(m_nInitInitMode, m_bInitSucFlag, TPixelBasedClustering::FAST, m_dbInitRatio, m_nInit1, m_nInit2, m_nInitIterFinal);
		
		TVector** ppMu = pInitPixelClustering->GetMean();
		TMatrix** ppCovInv = pInitPixelClustering->GetCovInv();
		double* dbPrior = pInitPixelClustering->GetPrior();
		double* dbCovDet = pInitPixelClustering->GetCovDet();

		for ( int i = 0; i < m_nClassCnt; i++ )
		{
			m_ppMu[i]->CopyFrom(ppMu[i]);
			m_ppCovInv[i]->CopyFrom(ppCovInv[i]);
			m_dbPrior[i] = dbPrior[i];
			m_dbCovDet[i] = dbCovDet[i];
			m_dbLnPrior[i] = log(m_dbPrior[i]);
			m_dbLnSigma[i] = log(m_dbCovDet[i]);
		}

		if (pRet)
			delete pRet; // this HAS to be deleted.

		delete pInitPixelClustering;
            
		m_bInitSucFlag = FinalClassify();

	}

	if ( m_nInitMode == RAND2_REGION )
	{
		// Find the range and variance of each feature
		int nFeatureCnt = m_pGraph->GetFeatureCnt();
		TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
		TMonoImage* pMask = m_pGraph->GetMask();

		double* Mu = new double[nFeatureCnt];
		double* Variance = new double[nFeatureCnt];
		int i;
		for ( i = 0; i < nFeatureCnt; i++ )
			Mu[i] = ppFeatures[i]->GetIntensityMean(pMask);

		int j, k, nCnt = 0;
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
		}			
		dbAvgVar /= nFeatureCnt;
		
		GenerateRandomLabel();
		m_bInitSucFlag = UpdateParam(dbAvgVar);

		delete[] Mu;
		delete[] Variance;
	}

	if ( m_nInitMode == KM_REGION )
	{
		TRegionBasedClustering* pInitRegionClustering = new TRegionBasedClustering(m_pGraph, m_nClassCnt);

		TGrayImage<int>* pRet = 0;

		pRet = pInitRegionClustering->RegionKMClustering(m_nInitInitMode, m_bInitSucFlag, m_dbInitRatio, m_nInit1, m_nInit2, m_nInitIterFinal);

	    if ( pRet )
		{
			int i;
			TVector** ppMu = pInitRegionClustering->GetMean();
			TMatrix** ppCovInv = pInitRegionClustering->GetCovInv();
			double* dbPrior = pInitRegionClustering->GetPrior();
			double* dbCovDet = pInitRegionClustering->GetCovDet();

			for ( i = 0; i < m_nClassCnt; i++ )
			{
				m_ppMu[i]->CopyFrom(ppMu[i]);
				m_ppCovInv[i]->CopyFrom(ppCovInv[i]);
				m_dbPrior[i] = dbPrior[i];
				m_dbCovDet[i] = dbCovDet[i];
				m_dbLnPrior[i] = log(m_dbPrior[i]);
				m_dbLnSigma[i] = log(m_dbCovDet[i]);
			}

			
			delete pRet; // this HAS to be deleted.
			delete pInitRegionClustering;			
 
			m_bInitSucFlag = FinalClassify();
			return true;
		}
		else
			return false;
	}

	// unrecognized initialization type.
	return false;
}

// 
// Compute the minimum Fisher criterion.
// 
double TGaussianRegionClustering::ComputeMinFisher()
{
	int nFeatureCnt = m_ppMu[0]->GetElementCnt();
	double dbMin = 100000000;  // a big enough number
	double dbVal;
	int i, j;
	if (nFeatureCnt == 1){
		for (i = 0; i < m_nClassCnt - 1; i ++){
		    for (j = i + 1; j < m_nClassCnt; j ++){
			    dbVal = (*m_ppMu[i])(0) - (*m_ppMu[j])(0);
			    dbVal *= dbVal;
			    dbVal /= (1 / (*m_ppCovInv[i])(0, 0) * (m_dbPrior[i]/(m_dbPrior[i]+m_dbPrior[j])) +            //Qin_2009_04_06
			             1 / (*m_ppCovInv[j])(0, 0) * (m_dbPrior[j]/(m_dbPrior[i]+m_dbPrior[j])));
			    //dbVal /= (1 / (*m_ppCovInv[i])(0, 0) * m_dbPrior[i] +            //Qin_2009_04_06
			    //         1 / (*m_ppCovInv[j])(0, 0) * m_dbPrior[j]);	
				//dbVal *= ((m_dbPrior[i]*m_dbPrior[j])/((m_dbPrior[i]+m_dbPrior[j])*(m_dbPrior[i]+m_dbPrior[j])));
				if (dbVal < dbMin)
				    dbMin = dbVal;
		    }
	    }
	}else if (nFeatureCnt > 1) //QIN_ADD
	{
		TMatrix* WSMat1 = new TMatrix(nFeatureCnt, nFeatureCnt);
		TMatrix* WSMat2 = new TMatrix(nFeatureCnt, nFeatureCnt);
		TVector* BSVec = new TVector(nFeatureCnt, false);
		TVector* BSVec_T = new TVector(nFeatureCnt, false);
		TMatrix *INVWSMat;
		for ( i = 0; i < m_nClassCnt - 1; i++ )
	    {
		    for ( j = i + 1; j < m_nClassCnt; j++ )
		    {
				for ( int k = 0; k < nFeatureCnt; k++ )
					(*BSVec)(k) = ((*m_ppMu[i])(k) - (*m_ppMu[j])(k));
				
				m_ppCovInv[i]->Inv(WSMat1);
				*WSMat1 *= (m_dbPrior[i]/(m_dbPrior[i]+m_dbPrior[j]));
				m_ppCovInv[j]->Inv(WSMat2);
				*WSMat2 *= (m_dbPrior[j]/(m_dbPrior[i]+m_dbPrior[j]));
				(*WSMat1) += (*WSMat2);

				//m_ppCovInv[i]->Inv(WSMat1);
				//*WSMat1 *= m_dbPrior[i];
				//m_ppCovInv[j]->Inv(WSMat2);
				//*WSMat2 *= m_dbPrior[j];
				//(*WSMat1) += (*WSMat2);

				INVWSMat = WSMat1->Inv();

				BSVec->Mul(BSVec_T, INVWSMat);
	            dbVal = BSVec_T->VecMul(BSVec);
				//dbVal *= ((m_dbPrior[i]*m_dbPrior[j])/((m_dbPrior[i]+m_dbPrior[j])*(m_dbPrior[i]+m_dbPrior[j])));

				if (dbVal < 0)
				{
					printf("Fisher calculation wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
					//_getchar_nolock();
					exit(1);
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
		exit(1);
	}

	return dbMin;
}

/// Update the class parameters/statistics by considering region internal statistics.
/**
 * This version of the function takes dbRegCoef, which is
 * regularization coefficient to use when a class is empty. The
 * coefficient is used to give that class valid statistics so that it
 * can hopefully become un-empty during subsequent relabeling.
 */
bool TGaussianRegionClustering::UpdateParam(double dbRegCoef)
{
	bool bSucFlag = true;

	int *nSampleCnt = new int[m_nClassCnt];
	int i, j, k;
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_ppMu[k]->SetZero();
		m_ppCovInv[k]->SetZero();
		nSampleCnt[k] = 0;
	}

	// Mean and Variance
	int nFeatureCnt = m_ppMu[0]->GetElementCnt();
	double *dbMu, **dbCoMu;
	TGaussianRegion* pCurrent;
	TVertex* pVertex = m_pGraph->GetVertices();
	while (pVertex)
	{
		pCurrent = (TGaussianRegion*) pVertex;
		k = pCurrent->m_nClassLabel;

		dbMu = m_ppMu[k]->m_ppData[0];
		dbCoMu = m_ppCovInv[k]->m_ppData;
		for ( i = 0; i < nFeatureCnt; i++ )
		{
			dbMu[i] += (*pCurrent->m_pMean)(i) * pCurrent->m_nPixelCnt;
			for ( j = i; j < nFeatureCnt; j++ )
				dbCoMu[i][j] += (*pCurrent->m_pCoMean)(i, j) * pCurrent->m_nPixelCnt;
		}
		nSampleCnt[k] += pCurrent->m_nPixelCnt;

		pVertex = pVertex->m_pNext;
	}

	int nTotalCnt = 0;
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_nClassPixelCnt[k] = nSampleCnt[k]; // store this for later use.

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
			bSucFlag = false;
			nSampleCnt[k]++;
		}

		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += dbRegCoef;

		m_dbCovDet[k] = sqrt(m_ppCovInv[k]->Determinant());
		m_dbLnSigma[k] = log(m_dbCovDet[k]);
	    m_ppCovInv[k]->Inv(m_ppCovInv[k]);

		nTotalCnt += nSampleCnt[k];
	}

	// Priors
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_dbPrior[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);
	    m_dbLnPrior[k] = log(m_dbPrior[k]); //QIN_ADD
	}

	delete[] nSampleCnt; //QIN_MOD

	return bSucFlag;
}


//
// Finally update the parameters for each class considering the inner boundary case
//
bool TGaussianRegionClustering::FinalUpdateParam()
{
	bool bSucFlag = true;

	/*TGrayImage<int>* tpMap = m_pGraph->GetMap();
	TMonoImage* pMask = m_pGraph->GetMask();
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();

	// Get maximum number of remaining segments
	TEstimateUtil* pEstUtil = new TEstimateUtil(0, 0);
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

	int nRow, nCol, l2, i, j, k;
	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;
			if ( (*pMap)(r, c) < 0 )
			{
				// Check if there are two different non-boundary labels in the neighborhood
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
		}*/

	TGrayImage<int>* pMap = GetLabeledGraphMap();

	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
	TMonoImage* pMask = m_pGraph->GetMask();

	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int r, c, i, j, k;

	for ( i = 0; i < m_nClassCnt; i++ )
	{
		m_ppMu[i]->SetZero();
		m_ppCovInv[i]->SetZero();
	}
	memset(m_dbPrior, 0, sizeof(double) * m_nClassCnt);

	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;
			if ( ( k = (*pMap)(r, c) ) < 0 )
				continue;

			for ( i = 0; i < nFeatureCnt; i++ )
			{
				(*m_ppMu[k])(i) += (*ppFeatures[i])(r, c);
				(*m_ppCovInv[k])(i, i) += (*ppFeatures[i])(r, c) * (*ppFeatures[i])(r, c);
				for ( j = i + 1; j < nFeatureCnt; j++ )
					(*m_ppCovInv[k])(i, j) += (*ppFeatures[i])(r, c) * (*ppFeatures[j])(r, c);
			}

			m_dbPrior[k]++;
		}

	double dbPriTotal = 0.0;
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_nClassPixelCnt[k] = (int) m_dbPrior[k]; // store this for later use.

		if ( m_dbPrior[k] > 0 )
		{
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppMu[k])(i) /= m_dbPrior[k];  // mean

			for ( i = 0; i < nFeatureCnt; i++ )  // covariance
			{
				(*m_ppCovInv[k])(i, i) = (*m_ppCovInv[k])(i, i) / m_dbPrior[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(i);
				for ( j = i + 1; j < nFeatureCnt; j++ )
				{
					(*m_ppCovInv[k])(i, j) = (*m_ppCovInv[k])(i, j) / m_dbPrior[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
					(*m_ppCovInv[k])(j, i) = (*m_ppCovInv[k])(i, j);
				}
			}
		}
		else
		{
			m_dbPrior[k]++;
			bSucFlag = false;
		};

		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += REG_COEF;	
		m_dbCovDet[k] = sqrt(m_ppCovInv[k]->Determinant());  // determinant
		m_dbLnSigma[k] = log(m_dbCovDet[k]);
		m_ppCovInv[k]->Inv(m_ppCovInv[k]);  // inverse of covariance
		dbPriTotal += m_dbPrior[k];
	}
	    
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_dbPrior[k] = m_dbPrior[k] / dbPriTotal;  // prior
		m_dbLnPrior[k] = log(m_dbPrior[k]); //QIN_ADD
	}

	delete pMap;
	return bSucFlag;
}

//
// Update the parameters for each class according to RAG
//
bool TGaussianRegionClustering::UpdateParam()
{
	bool bSucFlag = true;

	int *nSampleCnt = new int[m_nClassCnt];
	int i, j, k;
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_ppMu[k]->SetZero();
		m_ppCovInv[k]->SetZero();
		nSampleCnt[k] = 0;
	}

	// Mean and Variance
	int nFeatureCnt = m_ppMu[0]->GetElementCnt();
	double *dbMu, **dbCoMu;
	TGaussianRegion* pCurrent;
	TVertex* pVertex = m_pGraph->GetVertices();
	while (pVertex)
	{
		pCurrent = (TGaussianRegion*) pVertex;
		k = pCurrent->m_nClassLabel;

		dbMu = m_ppMu[k]->m_ppData[0];
		dbCoMu = m_ppCovInv[k]->m_ppData;
		for ( i = 0; i < nFeatureCnt; i++ )
		{
			dbMu[i] += (*pCurrent->m_pMean)(i) * pCurrent->m_nPixelCnt;
			for ( j = i; j < nFeatureCnt; j++ )
				dbCoMu[i][j] += (*pCurrent->m_pCoMean)(i, j) * pCurrent->m_nPixelCnt;
		}
		nSampleCnt[k] += pCurrent->m_nPixelCnt;

		pVertex = pVertex->m_pNext;
	}

	int nTotalCnt = 0;
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_nClassPixelCnt[k] = nSampleCnt[k]; // store this for later use.

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
			bSucFlag = false;
			nSampleCnt[k]++;
		}

		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += REG_COEF;

		m_dbCovDet[k] = sqrt(m_ppCovInv[k]->Determinant());
		m_dbLnSigma[k] = log(m_dbCovDet[k]);
	    m_ppCovInv[k]->Inv(m_ppCovInv[k]);

		nTotalCnt += nSampleCnt[k];
	}

	// Priors
	for ( k = 0; k < m_nClassCnt; k++ )
	{
		m_dbPrior[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);
	    m_dbLnPrior[k] = log(m_dbPrior[k]); //QIN_ADD
	}

	delete[] nSampleCnt; //QIN_MOD

	return bSucFlag;
}

//
// Update the parameters for each class according to the classification map
//
bool TGaussianRegionClustering::UpdateParam( TGrayImage<int>* pMap, int* nClassLabel )
{
	bool bSucFlag = true;

	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	TMonoImage* pMask = m_pGraph->GetMask();

	int i;
	for ( i = 0; i < m_nClassCnt; i++ )
	{
		m_ppMu[i]->SetZero();
		m_ppCovInv[i]->SetZero();
	}
	memset(m_dbPrior, 0, sizeof(double) * m_nClassCnt);

	int nHeight = pMap->Height();
	int nWidth = pMap->Width();

	int r, c, j, k;
	for ( r = 0; r < nHeight; r++ )
		for ( c = 0; c < nWidth; c++ )
		{
			if ( pMask && !(*pMask)(r, c) )
				continue;
			if ( ( k = (*pMap)(r, c) ) < 0 )
				continue;

			for ( i = 0; i < nFeatureCnt; i++ )
			{
				(*m_ppMu[k])(i) += (*ppFeatures[i])(r, c);
			    (*m_ppCovInv[k])(i, i) += (*ppFeatures[i])(r, c) * (*ppFeatures[i])(r, c);
				for ( j = i + 1; j < nFeatureCnt; j++ )
					(*m_ppCovInv[k])(i, j) += (*ppFeatures[i])(r, c) * (*ppFeatures[j])(r, c);
			}

			m_dbPrior[k]++;
		}

	double dbPriTotal = 0.0;
	for ( k = 0; k < m_nClassCnt; k++ )
	{

		m_nClassPixelCnt[k] = (int) m_dbPrior[k]; // store this for later use.

		if ( m_dbPrior[k] > 0 )
		{
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppMu[k])(i) /= m_dbPrior[k];  // mean

			for ( i = 0; i < nFeatureCnt; i++ )  // covariance
			{
				(*m_ppCovInv[k])(i, i) = (*m_ppCovInv[k])(i, i) / m_dbPrior[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(i);
				for ( j = i + 1; j < nFeatureCnt; j++ )
				{
					(*m_ppCovInv[k])(i, j) = (*m_ppCovInv[k])(i, j) / m_dbPrior[k] - (*m_ppMu[k])(i) * (*m_ppMu[k])(j);
					(*m_ppCovInv[k])(j, i) = (*m_ppCovInv[k])(i, j);
				}
			}
			nClassLabel[k] = 1;
		}
		else
		{
			bSucFlag = false;
			m_dbPrior[k]++;
			nClassLabel[k] = 0;
		}


		if ( m_ppCovInv[k]->Eigen() <= ZERO_THRESH )
			for ( i = 0; i < nFeatureCnt; i++ )
				(*m_ppCovInv[k])(i, i) += REG_COEF;	
		m_dbCovDet[k] = sqrt(m_ppCovInv[k]->Determinant());  // determinant
		m_ppCovInv[k]->Inv(m_ppCovInv[k]);  // inverse of covariance

		dbPriTotal += m_dbPrior[k];
	}
    
	for ( k = 0; k < m_nClassCnt; k++ )
		m_dbPrior[k] = m_dbPrior[k] / dbPriTotal;  // prior

	return bSucFlag;
}


// 
// Get the classification map
// by mapping each pixel of the image to its class category 
// 
TImage* TGaussianRegionClustering::GetClassMap()
{
	/*TGrayImage<int>* tpMap = m_pGraph->GetMap();
	TMonoImage* pMask = m_pGraph->GetMask();
	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();

	// Get maximum number of remaining segments
	TEstimateUtil* pEstUtil = new TEstimateUtil(0, 0);
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

	delete[] nLabels;*/
	
	TGrayImage<int>* pMap = GetLabeledGraphMap();

	TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
	TMonoImage* pMask = m_pGraph->GetMask();

	int nFeatureCnt = m_pGraph->GetFeatureCnt();
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	int l, r, c, i, nRow, nCol;

/*	int nRow, nCol, l2, i;
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
		}*/

	// Erase the border
	TVector* pVecFeature = new TVector(nFeatureCnt, false);  // feature vector
	double* dbVecProb = new double[m_nClassCnt];
	TVector* pVecTemp1 = new TVector(nFeatureCnt, false);  // a temperoray vector
	TVector* pVecTemp2 = new TVector(nFeatureCnt, false);  // a temperoray vector

	double dbMin;
	int nRefineCnt = 1, j, k, nMinIdx;

	int* nClassLabel = new int[m_nClassCnt]; // indicates existence of class in image.
	int* nPixelClassLabel = new int[m_nClassCnt]; // indicates existence of class among neighbours of pixel.
	int* nClassLabelUsed; // pointer to which of the above will be used when classifying boundary pixels.
	bool bIndepBoundary; // whether a pixel is an independent boundary pixel with no other pixels of any class beside it.
							// In these cases, the pixel must be classified on its own from amongst all the classes in image.

	TVector** ppMu = new TVector*[m_nClassCnt];
	TMatrix** ppCovInv = new TMatrix*[m_nClassCnt];
	for ( i = 0; i < m_nClassCnt; i++ )
	{
		ppMu[i] = new TVector(nFeatureCnt, false);
		ppCovInv[i] = new TMatrix(nFeatureCnt, nFeatureCnt);
	}

	double* dbPrior = new double[m_nClassCnt];
	double* dbCovDet = new double[m_nClassCnt];

	for ( int m = 0; m < nRefineCnt; m++ )  // several iterations to refine
	{
		for ( i = 0; i < m_nClassCnt; i++ )
		{
			ppMu[i]->SetZero();
			ppCovInv[i]->SetZero();
		}
		memset(dbPrior, 0, sizeof(double) * m_nClassCnt);

		for ( r = 0; r < nHeight; r++ )
			for ( c = 0; c < nWidth; c++ )
			{
				if ( pMask && !(*pMask)(r, c) )
					continue;
				if ( ( k = (*pMap)(r, c) ) < 0 )
					continue;

				for ( i = 0; i < nFeatureCnt; i++ )
				{
					(*ppMu[k])(i) += (*ppFeatures[i])(r, c);
					(*ppCovInv[k])(i, i) += (*ppFeatures[i])(r, c) * (*ppFeatures[i])(r, c);
					for ( j = i + 1; j < nFeatureCnt; j++ )
						(*ppCovInv[k])(i, j) += (*ppFeatures[i])(r, c) * (*ppFeatures[j])(r, c);
				}

				dbPrior[k]++;
			}

		double dbPriTotal = 0.0;
		for ( k = 0; k < m_nClassCnt; k++ )
		{
			if ( dbPrior[k] > 0 )
			{
				for ( i = 0; i < nFeatureCnt; i++ )
					(*ppMu[k])(i) /= dbPrior[k];  // mean

				for ( i = 0; i < nFeatureCnt; i++ )  // covariance
				{
					(*ppCovInv[k])(i, i) = (*ppCovInv[k])(i, i) / dbPrior[k] - (*ppMu[k])(i) * (*ppMu[k])(i);
					for ( j = i + 1; j < nFeatureCnt; j++ )
					{
						(*ppCovInv[k])(i, j) = (*ppCovInv[k])(i, j) / dbPrior[k] - (*ppMu[k])(i) * (*ppMu[k])(j);
						(*ppCovInv[k])(j, i) = (*ppCovInv[k])(i, j);
					}
				}
				nClassLabel[k] = 1;
			}
			else
			{
				nClassLabel[k] = 0;
				dbPrior[k]++;
			}

			if ( ppCovInv[k]->Eigen() <= ZERO_THRESH )
				for ( i = 0; i < nFeatureCnt; i++ )
					(*ppCovInv[k])(i, i) += REG_COEF;	
			dbCovDet[k] = sqrt(ppCovInv[k]->Determinant());  // determinant
			ppCovInv[k]->Inv(ppCovInv[k]);  // inverse of covariance

			dbPriTotal += dbPrior[k];
		}
	    
		for ( k = 0; k < m_nClassCnt; k++ )
			dbPrior[k] = dbPrior[k] / dbPriTotal;  // prior

		for ( r = 0; r < nHeight; r++ )
			for ( c = 0; c < nWidth; c++ )
			{
				if ( pMask && !(*pMask)(r, c) )
					continue;

				// Classify the boundary pixels using a pixel based
				// MRF model. First, check to see if it's a boundary pixel
				// and then if it is, determine what classes the neighbouring
				// pixels are. Choose from amongst these neighbouring classes
				// to classify the boundary pixel.

				memset(nPixelClassLabel, 0, sizeof(int) * m_nClassCnt); // reset the neighbouring class indicator list.

				if (IsBoundary(pMap, pMask, r, c, nPixelClassLabel, bIndepBoundary) && !m_bKeepBoundaries)
				{
					if (bIndepBoundary) // pixel has no neighbours with valid class, so classify it using all classes as candidates.
					{
						nClassLabelUsed = nClassLabel;
					}
					else // pixel has neighbours with valid class, so classify it using only those neighbouring classes.
					{
						nClassLabelUsed = nPixelClassLabel;
					}

					// Compute the probabilities
					for ( i = 0; i < nFeatureCnt; i++ )
						(*pVecFeature)(i) = (*ppFeatures[i])(r, c);

					for ( i = 0; i < m_nClassCnt; i++ )
					{
						if ( nClassLabelUsed [i] )
						{
							for ( j = 0; j < nFeatureCnt; j++ )
								(*pVecTemp1)(j) = (*pVecFeature)(j) - (*ppMu[i])(j);
							pVecTemp1->Mul(pVecTemp2, ppCovInv[i]);
							dbVecProb[i] = ( ( pVecTemp2->VecMul(pVecTemp1) / 2 + log(dbCovDet[i]) ) / nFeatureCnt );  //QIN_nFeatureCnt
						}
						else
							dbVecProb[i] = 999999999999;
					}

					if ( m_nType & TRegionClustering::VAFM )
					{
						for ( i = 0; i < m_nClassCnt; i++ )
							dbVecProb[i] *= ((TGIEPGraph*) m_pGraph)->GetCurAlpha();
					}

					// Account spatial context constraint
					for ( j = 0; j < 8; j++ )
					{
						nRow = r + neighbor8[j].Row;
						nCol = c + neighbor8[j].Col;

						if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
							continue;
						if (pMask && !(*pMask)(nRow, nCol))  // not valid
							continue;
						if ( ( l = (*pMap)(nRow, nCol) ) < 0 )  // also boundary
							continue;

						// Method to account for spatial context. First, assume that all dbVecProb at this pixel had m_nClassCnt*beta added.
						// If all neighbours are the same class, subtract m_nClassCnt*Beta from that class (i.e. no penalty)
						// If all neighbours are different, each class loses 1*Beta (i.e. (m_nClassCnt - 1)*beta penalty).
						// Note that this pixel can only be one of the classes represented by the neighbours since it's a boundary point.
						// Note that the dbVecProb never actually had the m_nClassCnt*beta penalties added explicitly, but because
						// we're looking for the minimum energy, it still works out.
						// i.e.		dbVecProb[i] + (m_nClassCnt - n_i) * beta > dbVecProb[k] + (m_nClassCnt - n_k) * beta
						//	implies
						//		dbVecProb[i] - n_i * beta > dbVecProb[k] - n_k * beta
						//	where n_i and n_k are the number of betas to subtract as calculated by this function for
						//	classes i and k.
						//
						// TODO: For REALLY noisy
						// images like cartoon_log, the
						// boundary pixel spatial
						// context isn't enough to
						// create smooth results at the
						// boundaries. Each boundary
						// pixel in a neighbourhood
						// might be assigned a
						// different class. Even when
						// the classes to be assigned
						// are limited to those that
						// occur in the immediate
						// neighbourhood, the results
						// still aren't smooth enough.
						// In these cases, a pure
						// spatial context (with no
						// feature model) model seems
						// to be much better. Perhaps
						// boundary pixels should only
						// use spatial context.

						dbVecProb[l] -= ((TGIEPGraph*) m_pGraph)->GetCurBeta(); // QIN_nFeatureCnt
						
					}

					// Determine the label
					dbMin = 999999999999;  // a large enough number
					for ( i = 0; i < m_nClassCnt; i++ )
					{
						if ( dbVecProb[i]< dbMin )
						{
							dbMin = dbVecProb[i];
							nMinIdx = i;
						}
					}

					(*pMap)(r, c) = nMinIdx;
		

				}
			}
	}

	for ( i = 0; i < m_nClassCnt; i++ )
	{
		delete ppMu[i];
		delete ppCovInv[i];
	}
	delete[] ppMu;
	delete[] ppCovInv;
	delete[] dbPrior;
	delete[] dbCovDet;

	delete[] nClassLabel;
	delete[] nPixelClassLabel;

	delete pVecFeature;
	delete[] dbVecProb;
	delete pVecTemp1;
	delete pVecTemp2;

	return pMap;
}

double TGaussianRegionClustering::GetFinalMLLEnergy( int& nRetCnt, TGrayImage<int>* pMap, TRegionClustering* pc, double dbBeta )
{
	TGaussianRegionClustering* pClassifier = (TGaussianRegionClustering*) pc;
		
	int nClusterCnt = pClassifier->m_nClassCnt;
	int nFeatureCnt = pClassifier->m_pGraph->GetFeatureCnt();
	TGrayImage<FLOAT>** ppFeatures = pClassifier->m_pGraph->GetImageSource();
	TMonoImage* pMask = pClassifier->m_pGraph->GetMask();

	TGrayImage<int>* pTmp = new TGrayImage<int>(pMap);
	TEstimateUtil* pEstimateUtil = new TEstimateUtil();
	nRetCnt = pEstimateUtil->RelabelMap(pTmp, pMask);
	delete pTmp;
	delete pEstimateUtil;


	int nWidth = pMap->Width();
	int nHeight = pMap->Height();

	double* dbLnSigma = new double[nClusterCnt];
	int i;
	for ( i = 0; i < nClusterCnt; i++ )
		dbLnSigma[i] = log(pClassifier->m_dbCovDet[i]);

	double dbTotalEnergy = 0;
	TVector* pVecFeature = new TVector(nFeatureCnt, false);  // feature vector
	TVector* pVecTemp1 = new TVector(nFeatureCnt, false);  // temporary feature vector
	TVector* pVecTemp2 = new TVector(nFeatureCnt, false);  // temporary feature vector
	int nRow, nCol, r, c, k;
	for ( nRow = 0; nRow < nHeight; nRow++ )
		for ( nCol = 0; nCol < nWidth; nCol++ )
		{

			if ( pMask && !(*pMask)(nRow, nCol) )
				continue;

			k = (*pMap)(nRow, nCol);

			for ( i = 0; i < nFeatureCnt; i++ )
				(*pVecFeature)(i) = (*ppFeatures[i])(nRow, nCol);

			for ( i = 0; i < nFeatureCnt; i++ )
				(*pVecTemp1)(i) = (*pVecFeature)(i) - (*pClassifier->m_ppMu[k])(i);
			pVecTemp1->Mul(pVecTemp2, pClassifier->m_ppCovInv[k]);
			dbTotalEnergy += ( ( pVecTemp2->VecMul(pVecTemp1) / 2 + dbLnSigma[k] ) / nFeatureCnt ); //QIN_nFeatureCnt
	
			for ( i = 0; i < 8; i++ )
			{
				r = nRow + neighbor8[i].Row;
				c = nCol + neighbor8[i].Col;
				if ( r < 0 || r >= nHeight || c < 0 || c >= nWidth )
					continue;
				if ( pMask && !(*pMask)(r, c) )
					continue;

				if ( (*pMap)(r, c) != k )
					dbTotalEnergy += dbBeta; //QIN_nFeatureCnt
			}

		}

	delete[] dbLnSigma;
	delete pVecFeature;
	delete pVecTemp1;
	delete pVecTemp2;

	return dbTotalEnergy;
}

/// Get segmentation statistics object that lists all the class statistics.
TSegmentationStatistics* TGaussianRegionClustering::GetSegStats(TGrayImage<int>* pMap)
{

	TGaussianSegmentationStatistics tscc(m_nClassCnt);

	int* nClassLabel = new int[m_nClassCnt];

	if (pMap)
	{
		// do an update based on the map.
		UpdateParam(pMap, nClassLabel);
	}

	TMatrix* temp = 0;

	for (int k = 0; k < m_nClassCnt; k++)
	{
		tscc.m_vMu[k] = *m_ppMu[k];
		temp = m_ppCovInv[k]->Inv();
		tscc.m_vCov[k] = *temp;
		delete temp;
		tscc.m_vPixelCounts[k] = m_nClassPixelCnt[k];
	}
	
	delete [] nClassLabel;

	//return new TGaussianSegmentationStatistics(tscc);
	return tscc.Clone();
}
