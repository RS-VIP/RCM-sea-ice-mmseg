/** \file
 *
 * \brief TRegion class
 *
 */

#ifndef TRegion_Class
#define TRegion_Class


#include "Graph.h"
#include "Boundary.h"
#include "GrayImage.h"
#include "MonoImage.h"

/// Represents a region in the RAG.
class TRegion : public TVertex
{
public:
	TRegion(TGraph* pGraph, int nLabel);

// Attributes
public:
	/// class label of the corresponding region
	int m_nClassLabel;

	/// count of the pixels in current region
	int m_nPixelCnt;

	/// mean of the X position
	double m_dbXMean;
	/// mean of the Y position	
	double m_dbYMean;

	/// mean of the edge strength in the region
	double m_dbEdgeMean;

	/// Bounding rectangle of the region segment
	int m_nTop, m_nBottom, m_nLeft, m_nRight;

	/// A flag indicating if it is a region at border (no matter grid or mask)
	BOOL m_bIsBorder;

// Operations
public:
	/// Compute the unary feature model weight of the vertex
	virtual double GetVertexWeight() = 0;

	/// Merge with a specified region
	virtual void MergeWith(TVertex* pNeighbor) = 0;

	/// Length of boundary of region
	int GetBoundaryLength();

    /// Accumulate pixel statistics during graph construction
    virtual void AccumulateStatistics(TGrayImage<FLOAT>** ppSrcImg, int r, int c);

    /// Finalize region statistics once graph construction is complete
    virtual void FinalizeStatistics();

protected:

	/// preserve the site to other boundary chain
	/**
	 *
	 * Preserve the site to a boundary chain associated with the specified
	 * region segment since it is still boundary site after merging
	 *
	 * Return
	 *  the associated boundary found
	 *  NULL: the associated boundary not found for current region segment
	 *
	 */
	TBoundary* PreserveNonDeletableSite(TRegion* pCurrent, int nLabel,
	                                    POINTex Pt);

	/// get a neighboring segment label on the specifed site. 
	/**
	 * If more than two different segments intersect, just choose the first found
	 *
	 * Param in
	 * Pt: position of the site
	 * 
	 * Param out
	 *  nLabel: the obtained label
	 *  
	 * Return
	 *  true: the specified site is a boundary site (with at least one
	 *   different neighbor; 
	 *  false: not boundary site
	 * 
	 */
	BOOL GetNeighboringRegionLabels(int& nLabel, POINTex Pt); 

// Implementation
public:
	virtual ~TRegion();

};

#endif

/////////////////////////////////////////////////////////////////////////////
