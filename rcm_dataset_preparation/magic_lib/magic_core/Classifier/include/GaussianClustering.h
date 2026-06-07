/** \file
 *
 * \brief TGaussianClustering class.
 *
 */

#ifndef TGaussianClustering_Class
#define TGaussianClustering_Class

//#include "ProcObject.h"
#include "SimpleClustering.h"
#include "Matrix.h"
#include "Vector.h"

/// Base class for some (but not all) clustering methods with Gaussian statistics.
/**
 * See subclasses for which specific clustering methods, there's quite a mix of
 * them.
 *
 */
class TGaussianClustering : public TSimpleClustering
{
public:
	TGaussianClustering(int nClusterCnt);

protected:

	/// number of features (will be set by derived classes).
	int m_nFeatureCnt;	

	/// means of the clusters
	TVector** m_ppMu;
	/// covariance matrice inversion of the clusters	
	TMatrix** m_ppCovInv;
	/// determinant of covariance matrices
	double* m_dbDeterm;
	/// prior probs of the clusters
	double* m_dbPri;


// Operations
public:

	TVector** GetMean() { return m_ppMu; };
	TMatrix** GetCovInv() { return m_ppCovInv; };
	double* GetPrior() { return m_dbPri; };
	double* GetCovDet() { return m_dbDeterm; };

	/// Get segmentation statistics object that lists all the class statistics.
	/**
	 * The statistics will not be updated before being returned
	 */
	virtual TSegmentationStatistics* GetSegStats();

	/// Get segmentation statistics object that lists all the class statistics.
	/**
	 * The statistics will first be updated by the classes in pMap.
	 */
	virtual TSegmentationStatistics* GetSegStats(TGrayImage<int>* pMap) { return 0; }

protected:

	/// Allocate memory.
	bool Alloc();

	/// Memeory manipulation
	TVector** AllocVectors( int nClusterCnt, int nFeatureCnt );
	/// Memeory manipulation	
	TMatrix** AllocMatrice( int nClusterCnt, int nFeatureCnt );

	/// Zero all provided vectors.
	void SetZeroVector( TVector** ppVectors, int nClusterCnt );
	/// Zero all provided matrices.	
	void SetZeroMatrix( TMatrix** ppMatrice, int nClusterCnt );
	/// Clear everything.
	void Clear();

	/// Arrange the clusters in the ascending order of the first mean feature value, so that the resulting map looks more pleasant
	void SortCluster();

	/// Print the cluster parameters to Console
	void ShowParam();

// Implementation
public:
	virtual ~TGaussianClustering();

};

#endif

/////////////////////////////////////////////////////////////////////////////
