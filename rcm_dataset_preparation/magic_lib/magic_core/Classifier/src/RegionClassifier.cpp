// RegionClassifier.cpp : implementation of the TRegionClassifier class
//
// Base class for classifying the regions
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "RegionClassifier.h"


//
// Constructor
// 
TRegionClassifier::TRegionClassifier()
{
	m_pGraph = NULL;
}

//
// Destructor
// 
TRegionClassifier::~TRegionClassifier()
{
}

// 
// Compute the overall energy corresponding to the current configuration
// 
double TRegionClassifier::OverallEnergy()
{
	return 0;
}

// 
// Intermediate classification
//
bool TRegionClassifier::IntermediateClassify()
{
	// Break the links between regions belonging different classes, so
	// that they cannot be merged afterwards
	TEdge* pCur = m_pGraph->GetEdges();
	while (pCur)
	{
		if (((TRegion*) pCur->m_pNeighbor1)->m_nClassLabel !=
			((TRegion*) pCur->m_pNeighbor2)->m_nClassLabel)
			((TBoundary*) pCur)->m_bActive = false;  // not active
		else
			((TBoundary*) pCur)->m_bActive = true;  // active
		pCur = pCur->m_pNext;
	}
	return true; 
}

// 
// Final classification
// 
bool TRegionClassifier::FinalClassify()
{
	return IntermediateClassify();
}
