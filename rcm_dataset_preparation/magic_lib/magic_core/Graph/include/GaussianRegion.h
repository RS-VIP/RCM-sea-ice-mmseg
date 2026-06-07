/** \file
 *
 * \author Peter Yu
 *
 * \brief TGaussianRegion class.
 *
 *
 */

#ifndef TGaussianRegion_Class
#define TGaussianRegion_Class

#include "Graph.h"
#include "Vector.h"
#include "Region.h"

/// A region in the RAG with Gaussian statistics
/** 
 *  This class represents one region in the Region Adjacency Graph which
 *  contains a mean feature vector and a matrix of co-values (from which the
 *  sample covariance is computed)
 */
class TGaussianRegion : public TRegion
{
public:
	TGaussianRegion(TGraph* pGraph, int nLabel, int nFeatureCnt);

// Attributes
public:
	/// mean vector
	TVector* m_pMean;  
	/// mean matrix of co-values from which sample covariance can be recovered.
	TMatrix* m_pCoMean;

// Operations
public:

    /// Accumulate pixel statistics during graph construction
    virtual void AccumulateStatistics(TGrayImage<FLOAT>** ppSrcImg, int r, int c);

    /// Finalize region statistics once graph construction is complete
    virtual void FinalizeStatistics();

	/// Compute the unary feature model weight of the vertex
	virtual double GetVertexWeight();

	/// Merge with a specified region
	/**
	 * This will erase the border between the two regions in the watershed
	 * map, compute the new region center, new region statistics and
	 * relabels the combined region in the watershed map with the right
	 * value.
	 */
	virtual void MergeWith(TVertex* pNeighbor);

// Implementation
public:
	virtual ~TGaussianRegion();

};

#endif

/////////////////////////////////////////////////////////////////////////////
