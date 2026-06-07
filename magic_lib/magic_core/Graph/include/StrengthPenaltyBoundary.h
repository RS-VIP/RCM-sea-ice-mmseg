/** \file
 *
 * \brief TStrengthPenaltyBoundary class
 *
 */

#ifndef TStrengthPenaltyBoundary_Class
#define TStrengthPenaltyBoundary_Class

#include "Boundary.h"

/// Base class for region boundary with strength penalty.
class TStrengthPenaltyBoundary : public TBoundary
{
public:
	TStrengthPenaltyBoundary(TGraph* pGraph, int nLabel1, int nLabel2,
	                         TChainCode* pChain, BOOL bOnGrid);

// Attributes
public:
	/// edge penalty energy
	double m_dbPenaltyEnergy;
	/// region energy difference between before merging and after merging
	double m_dbRegionEnergyDiff;

// Operations
public:
	/// Obtain the weight (energy difference between after merging and before merging) associated with current edge
	virtual double GetEdgeWeight();
	/// add a site to this boundary.
	virtual void AddSite(POINTex Pt, int nDir);
	/// Merge with specified edge
	virtual void MergeWith(TEdge* pEdge);

	/// compute feature model region energy difference between before merging and after merging
	virtual void ComputeRegionDiffEnergy() = 0;
	
        /// compute edge penalty energy	
	void ComputePenaltyEnergy();

protected:
	/// compute the edge penalty energy at a given site
	virtual double GetPenaltyEnergyAtSite(POINTex Pt, int nDir) = 0; 

// Implementation
public:
	virtual ~TStrengthPenaltyBoundary();

};

 
#endif

/////////////////////////////////////////////////////////////////////////////
