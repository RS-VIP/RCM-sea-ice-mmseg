/** \file
 *
 * \brief TRegionClustering class
 *
 */


#ifndef TRegionClustering_Class
#define TRegionClustering_Class

#include "MarkovNetClassifier.h"
#include "EstimateUtil.h"
#include "GrayImage.h"

/// Base class for clustering of regions using MRF on GIEP graph for spatial context.
/**
 * This is the base class that performs the MIRGS clustering.
 */
class TRegionClustering : public TMarkovNetClassifier
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
	TRegionClustering(int nType, int nInitMode, int nInitInitMode, double dbRatio, int nInit1, int nInit2, int nIterFinal,
					  int nMode, double dbInitT, int nIterations, int nClusterCnt);

// Attributes
public:

	/// Initialization methods.
	enum InitModeType { RAND1, RAND2, KM_PIXEL, GMM_PIXEL, KMGMM_PIXEL, RAND2_REGION, KM_REGION, GMM_REGION, KMGMM_REGION };
	/// Type of clustering cost function.
	enum ClusteringCostModeType { GIEP = 0x01, VAFM = 0x02 };

	/// Was initialization succesful?
	bool m_bInitSucFlag;

protected:

	/// Initialization mode.
	int m_nInitMode;

	/// Init mode of inititialization algorithm.
	int m_nInitInitMode;

	/// Parameter for the initialization algorithm, controls the stopping criterion.
	double m_dbInitRatio;


	//The initialization clustering algorithm will be run nInit2 times,
	//each time starting from a new random starting point. Each time, the
	//clustering algorithm will be run for nInit1 iterations. Out of the
	//nInit2 times, the result with the lowest energy is chosen.

	/// Parameter for the initialization algorithm, the number of iterations of the initialization clustering algorithm.
	int m_nInit1;

	/// Parameter for the initialization algorithm, how many times the clustering algorithm is executed
	int m_nInit2;

	// After the lowest energy initilization is chosen, the
	// clustering algorithm is further run for this many interations. This
	// last result is used as the starting point for IRGS.
	int m_nInitIterFinal;

	/// Type of clustering cost function.
	int m_nType;

// Operations
public:

	/// Initialize by generating an initial clustering (via various methods).	
	virtual bool Init() = 0;

	/// final classification
	virtual bool FinalClassify();

	/// Update the class parameters/statistics by considering region internal statistics.
	virtual bool UpdateParam() = 0;

	/// Update the class parameters/statistics by considering both region
	//internal statistics AND those of the boundary pixels.
	/**
	 * Obviously this is only used after the boundary pixels have been
	 * classified in the final clustering stage.
	 */
	virtual bool FinalUpdateParam() = 0;

	/// compute the minimum Fisher criterion (or substitute thereof)
	// TODO: name it more descriptively since for other clustering distance
	// equation / statistics, it may not be a fisher criterion.
	virtual double ComputeMinFisher() = 0;  
	virtual TImage* GetClassMap() = 0;

	int GetInitMethod() { return m_nInitMode; };
	int GetType() { return m_nType; };

	/// Final value of the MIRGS energy function.	
	virtual double GetFinalMLLEnergy( int& nRetCnt, TGrayImage<int>* pMap, TRegionClustering* pClassifier, double dbBeta ) = 0;

protected:

	/// Unary energy (i.e. feature model energy) for assigning label nClassIndex to pRegion.
	/**
	 * Param in
	 * pp: the region being classified
	 * nClassindex: the current class investigated
	 */	
	virtual double GetUnaryEnergy(TRegion* pRegion, int nClassIndex) = 0;

	/// Binary energy (i.e. spatial context with one other region connected to pRegion) for assigning label nClassIndex to pRegion
	/**
	 * This function is called for as many times as pRegion has neighbours,
	 * each time pBoundary is the edge between pRegion and the neighbour.
	 *
	 * Param in
	 * pp: the region being classified
	 * pBoundary: the boundary connecting the two regions of interest
	 * nClassIndex: the current class investigated
	 */
	virtual double GetBinaryEnergy(TRegion* pRegion, TBoundary* pBoundary, int nClassIndex) = 0;

	/// Randomly assign labels to regions for init purposes.
	void GenerateRandomLabel();

	/// Given a classification map pMap, update the class statistics based on labels in that map.
	/**
	 * nClassLabel is an array that on return indicates which classes are
	 * not empty.
	 */
	virtual bool UpdateParam( TGrayImage<int>* pMap, int* nClassLabel ) = 0;

	// Old IsBoundary method. Kept in case needed again.
	// virtual bool IsBoundary(TGrayImage<int>* pMap, int nRow, int nCol);

	/// Checks a pixel location is a boundary point. 
	/** Also, populates nClassLabel array
	 * to indicate which classes are adjacent to pixel. nClassLabel must be
	 * already allocated. bIndepBoundary is true when the pixel is a
	 * boundary but has no neighbours with a valid class.
	 */
	virtual bool IsBoundary(TGrayImage<int>* pMap, TMonoImage* pMask, int nRow, int nCol, int*& nClassLabel, bool& bIndepBoundary);

	/// Gets map where each region has its proper label applied but boundary pixels are still not labeled.
	virtual TGrayImage<int>* GetLabeledGraphMap();



// Implementation
public:
	virtual ~TRegionClustering();

};

#endif

/////////////////////////////////////////////////////////////////////////////
