/** \file
 *
 * \brief TGaussianStrengthPenaltyBoundary class
 *
 */

#ifndef TGaussianStrengthPenaltyBoundary_Class
#define TGaussianStrengthPenaltyBoundary_Class

#include "StrengthPenaltyBoundary.h"
#include "GaussianRegion.h"


/// Region boundary with strength penalty for Gaussian statistics.
class TGaussianStrengthPenaltyBoundary : public TStrengthPenaltyBoundary
{
public:
	TGaussianStrengthPenaltyBoundary(TGraph* pGraph, int nLabel1, int nLabel2,
	                         TChainCode* pChain, BOOL bOnGrid);
// Operations
public:

	/// compute region energy difference between before merging and after merging
	virtual void ComputeRegionDiffEnergy();

	//Compute gaussian difference between two regions.  
	static double GaussianRegionDifferenceEnergy(TGaussianRegion *region1, TGaussianRegion *region2);

protected:
	/// compute the edge penalty energy at a given site
	virtual double GetPenaltyEnergyAtSite(POINTex Pt, int nDir);

// Implementation
public:
	virtual ~TGaussianStrengthPenaltyBoundary();

};

 
#endif

/////////////////////////////////////////////////////////////////////////////
