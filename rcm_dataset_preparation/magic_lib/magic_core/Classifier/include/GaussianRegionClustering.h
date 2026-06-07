/** \file
 *
 * \brief TGaussianRegionClustering class
 *
 */

#ifndef TGaussianRegionClustering_Class
#define TGaussianRegionClustering_Class

#include "RegionClustering.h"
#include "EstimateUtil.h"
#include "GrayImage.h"

/// Clustering of Gaussian regions using MRF on GIEP graph for spatial context.
/**
 * This is a subclass that performs the MIRGS clustering with Gaussian
 * statistics. This is the normal MIRGS clustering algorithm.
 */
class TGaussianRegionClustering : public TRegionClustering
{
public:
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
	TGaussianRegionClustering(int nType, int nInitMode, int nInitInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal,
					  int nMode, double dbInitT, int nIterations, int nClusterCnt, bool bKeepBoundaries);

protected:
	
	/// means of the clusters
	TVector** m_ppMu;  
	/// covariance inversion of the clusters
	TMatrix** m_ppCovInv;
	/// prior probs of the clusters
	double* m_dbPrior;
	/// determinant of covariance
	double* m_dbCovDet; 

	/// Logs of standard deviation and priors, for speed consideration
	double* m_dbLnSigma;
	/// Logs of standard deviation and priors, for speed consideration
	double* m_dbLnPrior;

	//Set this flag to true to keep the boundary pixels and skip final classification
	bool m_bKeepBoundaries;

// Operations
public:

    bool SetMu(TVector** ppClusterMeans);
    bool SetCov(TMatrix** ppClusterCovInv);

	/// Initialize by generating an initial clustering (via various methods).
	virtual bool Init();

	/// Update the class parameters/statistics by considering region internal statistics.	
	virtual bool UpdateParam();

	/// Update the class parameters/statistics by considering both region
	/// internal statistics AND those of the boundary pixels.
	/**
	 * Obviously this is only used after the boundary pixels have been
	 * classified in the final clustering stage.
	 */
	virtual bool FinalUpdateParam();

	/// compute the minimum Fisher criterion
	virtual double ComputeMinFisher();  

	/// Gets the classification map based on current labeling of regions.
	/**
	 * It first assigns the pixels of each region to the region's label.
	 * Then, it assigns the boundary pixels to a class by classifying each
	 * pixel individually with the polarimetric IRGS cost function.
	 */
	virtual TImage* GetClassMap();

	TVector** GetMean() { return m_ppMu; };
	TMatrix** GetCovInv() { return m_ppCovInv; };
	double* GetPrior() { return m_dbPrior; };
	double* GetCovDet() { return m_dbCovDet; };

	/// Final value of the MIRGS energy function.
	virtual double GetFinalMLLEnergy( int& nRetCnt, TGrayImage<int>* pMap, TRegionClustering* pClassifier, double dbBeta );

	/// Get segmentation statistics object that lists all the class statistics.
	virtual TSegmentationStatistics* GetSegStats(TGrayImage<int>* pMap);



	/// Unary energy (i.e. feature model energy) for assigning label k to pRegion.
	/**
	 * Param in
	 * pp: the region being classified
	 * k: the current class investigated
	 */	
	virtual double GetUnaryEnergy(TRegion* pRegion, int nClassIndex);

	/// Binary energy (i.e. spatial context with one other region connected to pRegion) for assigning label k to pRegion
	/**
	 * This function is called for as many times as pRegion has neighbours,
	 * each time pBoundary is the edge between pRegion and the neighbour.
	 *
	 * Param in
	 * pp: the region being classified
	 * pBoundary: the boundary connecting the two regions of interest
	 * k: the current class investigated
	 */	
	virtual double GetBinaryEnergy(TRegion* pRegion, TBoundary* pBoundary, int nClassIndex);

protected:
	/// Given a classification map pMap, update the class statistics based on labels in that map.
	/**
	 * nClassLabel is an array that on return indicates which classes are
	 * not empty.
	 */	
	virtual bool UpdateParam( TGrayImage<int>* pMap, int* nClassLabel );

	/// Update the class parameters/statistics by considering region internal statistics.
	/**
	 * This version of the function takes dbRegCoef, which is
	 * regularization coefficient to use when a class is empty. The
	 * coefficient is used to give that class valid statistics so that it
	 * can hopefully become un-empty during subsequent relabeling.
	 */
	virtual bool UpdateParam( double dbRegCoef );

// Implementation
public:
	virtual ~TGaussianRegionClustering();

};

#endif

/////////////////////////////////////////////////////////////////////////////
