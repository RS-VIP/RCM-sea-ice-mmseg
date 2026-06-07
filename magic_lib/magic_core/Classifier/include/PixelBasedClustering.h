/** \file
 *
 * \brief TPixelBasedClustering class.
 *
 */

#ifndef TPixelBasedClustering_Class
#define TPixelBasedClustering_Class

#include "GaussianClustering.h"
#include "GrayImage.h"
#include "MonoImage.h"
#include "AnnealingUtil.h"
#include "Sampler.h"

/// Pixel based clustering methods (both simple and those that have spatial context) using Gaussian statistics.
/**
 * This class contains pixel based clustering methods. There are simple ones,
 * like K-means, GMM, KM-GMM that don't consider spatial context and there are
 * spatial context ones like CMLL and VMLL.
 *
 * The simple clustering methods are mainly used for initialization but can
 * also be called to create segmentation results that do not consider spatial
 * context. Note that while MIRGS uses the result from this class as its
 * initial solution, the clustering methods in this class must also be
 * initialized somehow (e.g. for K-means, by random assignment of labels
 * first).
 *
 * TODO: Consider splitting this into subclasses.
 *
 */
class TPixelBasedClustering : public TGaussianClustering
{
public:
	TPixelBasedClustering(TMonoImage* pMask, 
		TGrayImage<FLOAT>** ppSrcImg, 
		int nClusterCnt,
		int nFeatureCnt);

// Attributes
public:
	/// Fast mode in which histogram is used to speed up the clustering process (only for univariate data)
	enum enModelMode { FAST = 0x01, NORMAL = 0x02};
	/// How should Beta for MLL methods be updated per iteration.
	enum enBetaMode { FIXED = 0x00, MLFIXED = 0x01, MLFISH = 0x02 };
	/// Initialization mode.
	enum enInitMode { RAND1, RAND2, KM_PIXEL, GMM_PIXEL, KMGMM_PIXEL };

protected:
	/// Feature image(s)
	TGrayImage<FLOAT>** m_ppSrcImg;
	/// Mask
	TMonoImage* m_pMask;

// Operations
public:

	/// Initialization
	void PixelInitParam( int nInitMode, double* VecMin, double* VecMax, double* Variance );
	/// Initialization (for fast mode)
	void PixelInitParam( int nInitMode, int* pHist, int nHistBinCnt, double* Variance ); // For fast mode

	/// GMM Clustering
	TGrayImage<int>* GaussMixClustering(int nInitMode, bool& bSucFlag, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal);

	/// GMM clustering (univariate)
	// TGrayImage<int>* GaussMixClustering(int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>* pFeature, TMonoImage* pMask, int nClusterCnt, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal); 
	/// GMM clustering (multivariate)
	// TGrayImage<int>* GaussMixClustering(int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>** ppFeatures, int nFeatureCnt, TMonoImage* pMask, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal);

	/// KM clustering
	TGrayImage<int>* KMClustering(int nInitMode, bool& bSucFlag, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal);
	
	/// KM clustering (univariate)
	// TGrayImage<int>* KMClustering(int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>* pFeature, TMonoImage* pMask, int nClusterCnt, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal); // univariate case
	/// KM Clustering (multivariate)
	// TGrayImage<int>* KMClustering(int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>** ppFeatures, int nFeatureCnt, TMonoImage* pMask, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal); //multivariate case

	/// OTSU clustering with thresholding //Shuhrat otsu
	TGrayImage<int>* OTSUClustering(bool& bSucFlag, int& threshold, bool bFast);
		
	/// OTSU clustering (univariate)
	// TGrayImage<int>* OTSUClustering(); // univariate case
	/// OTSU Clustering (multivariate)
	// TGrayImage<int>* OTSUClustering(); //multivariate case //Shuhrat otsu

	/// OTSU clustering automaticaly //Shuhrat otsu
	TGrayImage<int>* OTSUClustering(bool& bSucFlag, int& threshold);
		
	/// OTSU clustering (univariate)
	// TGrayImage<int>* OTSUClustering(); // univariate case
	/// OTSU Clustering (multivariate)
	// TGrayImage<int>* OTSUClustering(); //multivariate case //Shuhrat otsu

	/// KMGMM clustering
	TGrayImage<int>* KMGMMClustering(int nInitMode, bool& bSucFlag, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal);
	/// KMGMM clustering (univariate)
	// TGrayImage<int>* KMGMMClustering(int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>* pFeature, TMonoImage* pMask, int nClusterCnt, int nMode, double dbRatio, int nInit1, int nInit2, int nIterFinal); // univariate case
	/// KMGMM clustering (multivariate)
	// TGrayImage<int>* KMGMMClustering(int nInitMode, bool& bSucFlag, TGrayImage<FLOAT>** ppFeatures, int nFeatureCnt, TMonoImage* pMask, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal); //multivariate case

	///CMLL clustering
	TGrayImage<int>* CMLLClustering(int nIter, int nBetaMode, double dbBeta1, double dbBeta2, int nInitMode, int nInitInitMode, bool& bInitSucFlag, double dbRatio, int nInit1, int nInit2, int nIterFinal, double dbSAInitT, bool bOutIntermResult);

	///VMLL clustering
	TGrayImage<int>* VMLLClustering(int nIter, int nBetaMode, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2, int nInitMode, int nInitInitMode, bool& bInitSucFlag, double dbRatio, int nInit1, int nInit2, int nIterFinal, double dbSAInitT, bool bOutIntermResult);

	/// Get segmentation statistics object that lists all the class statistics.
	virtual TSegmentationStatistics* GetSegStats(TGrayImage<int>* pMap);

protected:
	/// Memeory manipulation
	// bool Alloc(int nClusterCnt, int nFeatureCnt);

	// The following are for MLL based clustering algorithms
	void GenerateRandomMap(TGrayImage<int>* pMap);
	void GenerateRandomHist(int* nHistLabel, int nHistBinCnt);
	int ComputeEdgePairCnt(TGrayImage<int>* pMap, TMonoImage* pMask);
	/// Update Beta parameter
	double UpdateBeta(double dbEdgeRatio);
	/// Update class parameters based on map
	bool UpdateParam(TGrayImage<int>* pMap, double dbRegCoef);
	/// Compute minimum fisher criterion between all pairs of classes.
	double ComputeMinFisher();
	
	// **************************************** GMM related sub-routines **************************************** //
	/// Estimate cluster parameters by using GMM
	void GaussMixParamEst(int nInitMode, int* pHist, int nHistBinCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal);
	void GaussMixParamEst(int nInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal);

    	/// GMM clustering based on the current cluster parameters
	TGrayImage<int>* GaussMixClassify(bool& bSucFlag);
	
	/// The EM for GMM parameter estimation in a fast mode (i.e. using histogram)
	double GaussMixEMIterations(double dbRatio, int nIterations, int* pHist, int nHistBinCnt);

	/// The EM for GMM parameter estimation
	double GaussMixEMIterations(double dbRatio, int nIterations);
	
	/// The E step for GMM
	double GaussMixExpectation(TVector& VecProb, TVector& VecFeature);
	
	/// The M step for GMM
	void GaussMixMaximization(TVector** ppMu, TMatrix** ppCov, double* dbPri, TVector& VecProb, TVector& VecFeature);

	/// **************************************** KM related sub-routines **************************************** ///
	/// Estimate cluster parameters by using KM
	void KMParamEst(int nInitMode, int* pHist, int nHistBinCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal);
	void KMParamEst(int nInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal);

	/// The univariate KM operation in a fast mode (i.e. using histogram)
	double KMIterations(double dbRatio, int nIterations, int* pHist, int nHistBinCnt);

	/// The multivariate KM operation
	double KMIterations(double dbRatio, int nIterations);

 	/// KM clustering based on the current cluster parameters
	TGrayImage<int>* KMClassify(bool& bSucFlag);

		/// KM clustering based on the current cluster parameters
	TGrayImage<int>* OTSUClassify(bool& bSucFlag,int threshold);

	// **************************************** KMGMM related sub-routines **************************************** //
	/// Estimate cluster parameters by using KMGMM
	void KMGMMParamEst(int nInitMode, int* pHist, int nHistBinCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal);

	/// Estimate cluster parameters by using KMGMM	
	void KMGMMParamEst(int nInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal);
	
	// **************************************** CMLL related sub-routines **************************************** //
	/// CMLL Annealing
	double CMLLAnnealing( TGrayImage<int>*& pMap, double dbBeta, TAnnealingUtil* pAnnealUtil, TSampler* pScanUtil);

	// **************************************** VMLL related sub-routines **************************************** //
	/// VMLL Annealing
	double VMLLAnnealing( TGrayImage<int>*& pMap, double dbAlpha, double dbBeta, TAnnealingUtil* pAnnealUtil, TSampler* pScanUtil);
	
	/// Update alpha parameter.
	double UpdateAlpha(int nIter, double dbAlpha1, double dbAlpha2);



// Implementation
public:
	virtual ~TPixelBasedClustering();

};

#endif

/////////////////////////////////////////////////////////////////////////////
