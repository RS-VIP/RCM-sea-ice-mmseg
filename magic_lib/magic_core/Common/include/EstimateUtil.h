/** \file
 *
 * \brief TEstimateUtil class
 *
 */

#ifndef TEstimateUtil_Class
#define TEstimateUtil_Class

//#include "ProcObject.h"
#include "GrayImage.h"
//#include "TrueColorImage.h"
#include "MonoImage.h"
#include "Matrix.h"
#include "Vector.h"

// sort function yoinked from matrix_gabor.h
void sort_(double *Y, int *I, double *A, int length);

/// Misc. estimation utilities, Gaussian mixture routines, some class map operations
/**
 * Estimate parameters for Gaussian mixtures
 * Several useful class map ops: consecutive relabeling, relabeling for display, maximum label, etc.
 */
class TEstimateUtil
{
public:
	TEstimateUtil();

// Attributes
public:
	/// fast mode that uses histogram to speed up the clustering
	enum model_mode {FAST       = 0x01};  

protected:
	// Storage for Gaussian parameters
	/// mean
	TVector** m_ppMu;
	/// covariance
	TMatrix** m_ppCovInv; 
	/// determinant of covariance
	double* m_dbDeterm; 
	/// prior
	double* m_dbPri;

	/// Number of classes
	int m_nClassCnt;
	/// Number of features
	int m_nFeatureCnt;

private:
	// The following are temporarily used for saving time in allocating and freeing memeory
	TVector* m_pVecTemp1;  // a temporary vector
	TVector* m_pVecTemp2;  // a temporary vector

// Operations
public:
	/// Gauss mixture Clustering
	TGrayImage<int>* GaussMixClustering(TGrayImage<FLOAT>* pFeature, TMonoImage* pMask, 
	                                    int nClassCnt, int nMode);
	/// Gauss mixture Clustering (multivariate)	
	TGrayImage<int>* GaussMixClustering(TGrayImage<FLOAT>** ppFeatures, TMonoImage* pMask,
	                                    int nClassCnt, int nFeatureCnt);
	/// Compute parameters for Gaussian mixture
	void GaussMixParamEst(TGrayImage<FLOAT>* pFeature, TMonoImage* pMask,
	                      int nClassCnt, int nMode);
	/// Compute parameters for Gaussian mixture, histogram edition.
	void GaussMixParamEst(int nHist[], int nHistBinCnt, int nClassCnt);
	/// Compute parameters for Gaussian mixture, multivariate.
	void GaussMixParamEst(TGrayImage<FLOAT>** ppFeatures, TMonoImage* pMask,
	                      int nClassCnt, int nFeatureCnt);
	/// Update the Gaussian mixture parameters
	BOOL GaussMixParamUpdate(TGrayImage<FLOAT>** ppFeatures, TMonoImage* pMask,
	                         TGrayImage<int>* pMap, int nClassCnt, int nFeatureCnt);
	/// Classify based on obtained Guass mixture parameters
	TGrayImage<int>* GaussMixClassify(TGrayImage<FLOAT>** ppFeatures, TMonoImage* pMask);
	/// Display the Gaussian mixture parameters
	void ShowGaussianMixParam();

	/// Generate a random label assignment
	void GenerateRandomMap(TGrayImage<int> *pMap, TMonoImage* pMask, int nClassCnt);
	/// Re-labelize the map to have consecutive labels.
	int RelabelMap(TGrayImage<int>* pMap, TMonoImage* pMask);

    int* RelabelMapAndReturnIndices(TGrayImage<int>* pMap, TMonoImage* pMask);


	/// Sort the map according to the average feature value in each region
	void SortMap(TGrayImage<int>*& pMap, TMonoImage* pMask, TGrayImage<FLOAT>** ppFeatImg, int nFeatureCnt);
	/// Display the classification map
	TGrayImage<int>* DispMap(TGrayImage<int>* pMap, TMonoImage* pMask, bool bAddBlackBox);
	/// Get the maximum label of the map
	int GetMaxLabel(TGrayImage<int>* pMap, TMonoImage* pMask);

	/// Get the number of positive labels of the map
	int GetNumPositiveLabels(TGrayImage<int>* pMap, TMonoImage* pMask);

	/// Erase the border
	void EraseBoundary(TGrayImage<int>* pMap, TGrayImage<FLOAT>** ppFeatImg, int nFeatureCnt, TMonoImage* pMask);
	/// Arrange the clusters in the ascending order of the first mean feature values, so that the resulting map looks more pleasant
	void SortCluster();
	/// Add black box to a result segmentation map
	TGrayImage<int>* AddBlackBox(TGrayImage<int>* pMap);

	// For public access
	TVector** GetMean() { return m_ppMu; };
	TMatrix** GetCovInv() { return m_ppCovInv; };
	double* GetPrior() { return m_dbPri; };
	double* GetCovDet() { return m_dbDeterm; };

	/// Estimate the local variance of the image (i.e. the variance of noise assuming constant gray level corrupted by noise in univariate case)
	double EstVar(TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt, TMonoImage* pMask); //QIN_MOD
	/// Find the mode of the distribution
	double FindMode(TGrayImage<FLOAT>* pImg, TMonoImage* pMask);

	// Memory manipulation
	TVector** AllocVectors(int nClassCnt, int nFeatureCnt);
	TMatrix** AllocMatrice(int nClassCnt, int nFeatureCnt);
	void SetZeroVector(TVector** ppVectors, int nClassCnt);
	void SetZeroMatrix(TMatrix** ppMatrice, int nClassCnt);

protected:
	/// Memeory manipulation
	void Alloc(int nClassCnt, int nFeatureCnt);
	void Clear();

	/// Initialize the Gaussian mixutre parameters randomly
	void InitGaussMixParam(TVector& VecMin, TVector& VecMax);
	
	/// The EM for Gaussian mixture parameter estimation in fast mode (i.e. using histogram)
	double GaussMixEMIterations(double dbRatio, int nIterations, int* pHist, int nHistBinCnt);

	/// The EM for Gaussian mixture parameter estimation
	double GaussMixEMIterations(double dbRatio, int nIterations, TGrayImage<FLOAT>** ppFeatures, TMonoImage* pMask);
	
	/// The expectation step for Gaussian mixture parameter estimation
	double GaussMixExpectation(TVector& VecProb, TVector& VecFeature);
	
	/// The maximization step for Gaussian mixture estimation
	void GaussMixMaximization(TVector** ppMu, TMatrix** ppCov, double* dbPri, TVector& VecProb, TVector& VecFeature);


	/// Check if the current site is at boundary
	BOOL IsBoundary(TGrayImage<int>* pMap, TMonoImage* pMask, int nRow, int nCol);

// Implementation
public:

    /// Generate consecutive mapping
    int* GenerateConsecutiveMapping(int& nMappingItemCnt,TGrayImage<int>* pMap, TMonoImage* pMask);

    /// Generate the LUT that maps nodes from before a merging step to nodes after a merging step
    int* GenerateMergingLUT(int& nMappingItemCnt, TGrayImage<int>* pMapOld, TGrayImage<int>* pMapNew, TMonoImage* pMask);

    virtual ~TEstimateUtil();

};

#endif

/////////////////////////////////////////////////////////////////////////////
