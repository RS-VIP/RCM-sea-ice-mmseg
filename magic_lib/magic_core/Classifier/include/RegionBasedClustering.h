/** \file
 *
 * \brief TRegionBasedClustering class.
 *
 * In this file, all the GMM and KMGMM methods are commented because they are
 * theoretically inappropriate for region-based processing (Alex Qin)
 *
 */

#ifndef TRegionBasedClustering_Class
#define TRegionBasedClustering_Class

#include "GaussianClustering.h"
#include "RegionAdjacencyGraph.h"

/// Region-based clustering methods with no spatial context using Gaussian statistics.
/** 
 * This class implements region-based clustering methods that do not use
 * spatial context using Gaussian statistics. This is mainly used for
 * initialization of MIRGS but can also be called to create segmentation results
 * that do not consider spatial context. Note that while MIRGS uses the result
 * from this class as its initial solution, the clustering methods in this
 * class must also be initialized somehow (e.g. for K-means, by random
 * assignment of labels first).
 *
 * In this class, all the GMM and KMGMM methods are commented because they are
 * theoretically inappropriate for region-based processing (Alex Qin)
 *
 */
class TRegionBasedClustering : public TGaussianClustering
{
public:
	/// Constructs RAG from the supplied image.
	TRegionBasedClustering(TGrayImage<int>* pMap, TMonoImage* pMask,
		                   TGrayImage<FLOAT>** ppSrcImg, int nClusterCnt, int nFeatureCnt);
public:
	/// Takes already constructed RAG
	TRegionBasedClustering(TRegionAdjacencyGraph* pGraph, int nClusterCnt);
// Attributes
public:
	/// Initialization mode for the clustering methods in this class.
	enum enInitMode { RAND1, RAND2 };

protected:
	
	TRegionAdjacencyGraph* m_pGraph; //region adjacency graph

	bool m_bGraphNew; // a flag indicating whether the graph is created in the class or passed in as a parameter

// Operations
public:
	/// region adjacency graph not for spatial context but to provide access to all regions.
	void SetGraph(TRegionAdjacencyGraph* pGraph) { m_pGraph = pGraph; };

	// *** Region-based GMM is commented out because it is theoretically inappropriate for region-based processing (Alex Qin)
	// Region-based GMM clustering
	// TGrayImage<int>* RegionGMMClustering(int nInitMode, bool& bSucFlag, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal);

	/// Region-based KM clustering
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
	TGrayImage<int>* RegionKMClustering(int nInitMode, bool& bSucFlag, double dbRatio, int nInit1, int nInit2, int nIterFinal);

	// *** Region-based GMM is commented out because it is theoretically inappropriate for region-based processing (Alex Qin)
	// Region-based KM + GMM clustering
	// TGrayImage<int>* RegionKMGMMClustering(int nInitMode, bool& bSucFlag, int nClusterCnt, double dbRatio, int nInit1, int nInit2, int nIterFinal);

	/// Gets a pointer to the region adjacency graph.
	TRegionAdjacencyGraph* GetGraph() { return m_pGraph; };

	/// Get segmentation statistics object that lists all the class statistics.
	virtual TSegmentationStatistics* GetSegStats(TGrayImage<int>* pMap);

protected:
	/// Memeory manipulation
	// bool Alloc( int nClusterCnt );

	/// Initialize cluster parameters
	void RegionInitParam(int nInitMode, double* VecMin, double* VecMax, double* Variance);
	/// Randomly assign the class label for each region
	void GenerateRandomLabel();
	/// Update cluster parameters according to the current region class-label configuration
	bool UpdateParam(double dbRegCoef);
    /// Final update cluster parameters according to region and boundary
	bool FinalUpdateParam(TGrayImage<int>* pMap, double dbRegCoef);

	// **************************************** Region-based GMM sub-routines **************************************** //
	// The iterations for region-based GMM by Expectation Maximization
	// double RegionGMMIteration(double dbRatio, int nIterations);

	// Final class label assignment
	// void RegionGMMFinalClassify(bool& bSucFlag);

	// Obtain the final classification map
    	// TGrayImage<int>* RegionGMMGetClassMap();
    
	// Assign boundary pixels to one of the two neighboring regions
    	// int RegionGMMBoundaryClassify(TVector* VecFeature, int nLabel1, int nLabel2);

	// Assign boundary pixels to one of the two neighboring regions
	// double RegionGMMBoundaryClassify(TVector* VecFeature, int nLabel);
	/// ********************************************************************************************************* ///
	
	// **************************************** Region-based KM sub-routines **************************************** //
	//
    	/// The iterations for region-based KM
	double RegionKMIteration(double dbRatio, int nIterations);

	/// Final class label assignment
	void RegionKMFinalClassify(bool& bSucFlag);

	/// Obtain the final class map
	TGrayImage<int>* RegionKMGetClassMap();

	// Assign boundary pixels to one of the two neighboring regions
	// int RegionKMBoundaryClassify(TVector* VecFeature, int nLabel1, int nLabel2);

	/// Assign boundary pixels to one of the two neighboring regions
	double RegionKMBoundaryClassify(TVector* VecFeature, int nLabel);

// Implementation
public:
	virtual ~TRegionBasedClustering();

};

#endif

/////////////////////////////////////////////////////////////////////////////
