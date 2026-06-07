// StrengthPenaltyBoundary.cpp : implementation of the TStrengthPenaltyBoundary class
//
// Class for region boundary with strength penalty
/////////////////////////////////////////////////////////////////////////////

//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
//#include <fstream>
#include "StrengthPenaltyBoundary.h"
#include "GIEPGraph.h"

// 
// Constructor
//
TStrengthPenaltyBoundary::TStrengthPenaltyBoundary(TGraph* pGraph,
                                                   int nLabel1, int nLabel2,
                                                   TChainCode* pChain,
                                                   BOOL bOnGrid)
	: TBoundary(pGraph, nLabel1, nLabel2, pChain, bOnGrid)
{
	m_dbPenaltyEnergy = -1;
	m_dbRegionEnergyDiff = -1;
}

// 
// Destructor
// 
TStrengthPenaltyBoundary::~TStrengthPenaltyBoundary()
{
	TChainCode* pCur = m_pChain;
	while (pCur)
	{
		m_pChain = m_pChain->m_pNext;
		delete pCur;
		pCur = m_pChain;
	}
}

// 
// Merge with a specified StrengthPenaltyBoundary
// 
void TStrengthPenaltyBoundary::MergeWith(TEdge* pEdge)
{
	//TBoundary::MergeWith(pEdge);

	// Update penalty energy and the difference between after
	// merging and before merging
	m_dbPenaltyEnergy += ((TStrengthPenaltyBoundary*) pEdge)->m_dbPenaltyEnergy;
	TBoundary::MergeWith(pEdge);
}

// 
// Add a site
// 
// Param in
//  Pt: the location of the site
//  nDir: the boundary direction of the site
// 
void TStrengthPenaltyBoundary::AddSite(POINTex Pt, int nDir)
{
	TBoundary::AddSite(Pt, nDir);

	m_dbPenaltyEnergy += GetPenaltyEnergyAtSite(Pt, nDir);
}

// 
// Compute edge penalty energy
// 
void TStrengthPenaltyBoundary::ComputePenaltyEnergy()
{
	if (!m_pGraph)
		return;

	double dbEnergy = 0;
	TChainCode *pChain = m_pChain;
	POINTex Pt;
	int nDir;

	while (pChain)
	{
		pChain->GetFirstSite(Pt, nDir);
		do {

		    dbEnergy += GetPenaltyEnergyAtSite(Pt, nDir);

		} while (pChain->GetNextSite(Pt, nDir));

		pChain = pChain->m_pNext;
	}

	m_dbPenaltyEnergy = dbEnergy * ((TGIEPGraph*) m_pGraph)->GetCurBeta();
}

// 
// Obtain the weight (energy difference between after merging and
// before merging) associated with current edge
// 
double TStrengthPenaltyBoundary::GetEdgeWeight()
{
	if (!m_bActive)
		return INIT_MAXDIST;  // A large enough number (acutually any positive number is ok)

	//if (m_dbRegionEnergyDiff < 0)  // not initialized
	if (m_dbRegionEnergyDiff == -1)  // not initialized
		ComputeRegionDiffEnergy();
	//if (m_dbPenaltyEnergy < 0)  // not initialized
	if (m_dbPenaltyEnergy == -1)  // not initialized
		ComputePenaltyEnergy();

	return ( m_dbRegionEnergyDiff - m_dbPenaltyEnergy ); //QIN_nFeatureCnt
}
