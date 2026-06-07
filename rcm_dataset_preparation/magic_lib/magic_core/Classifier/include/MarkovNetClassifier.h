/** \file
 *
 * \brief TMarkovNetClassifier class
 *
 */




#ifndef TMarkovNetClassifier_Class
#define TMarkovNetClassifier_Class

#include "AnnealingUtil.h"
#include "Sampler.h"
#include "RegionClassifier.h"
#include "Region.h"
#include "Boundary.h"

///Markov network classifier; similar to Markov Random Field but has irregular grid
class TMarkovNetClassifier : public TRegionClassifier
{
public:
	TMarkovNetClassifier( int nAnnealingMode, double dbInitT, int nIterCnt, int nClassCnt);

// Attributes
public:
	/// Annealing mode
	enum AnnealingModeType { SINGLE_STEP = 0x01,  // classification performed in single step
							NO_COOLING = 0x02 };  // no cooling, stay at a single temperature

protected:
	/// simulated annealing utility
	TAnnealingUtil* m_pAnnealUtil;  

	/// scan utility (sampler)
	TSampler* m_pScanUtil;

	/// initial temperature
	double m_dbInitT;

	/// number of iterations
	int m_nIterCnt;

	/// mode.
	int m_nAnnealingMode; 

	/// number of classes
	int m_nClassCnt;

	/// count of number of pixels in each class.
	int* m_nClassPixelCnt;

// Operations
public:
	/// initialize the classifier
	virtual bool Init();  

	/// intermidiate classification
	virtual bool IntermediateClassify();

	/// Overall energy
	virtual double OverallEnergy();

	/// Access annealing util
	TAnnealingUtil* GetAnnealUtil() { return m_pAnnealUtil; };

protected:
	/// Unary energy (i.e. feature model energy) for assigning label nClassIndex to pRegion.
	virtual double GetUnaryEnergy( TRegion* pRegion, int nClassIndex ) = 0;

	/// Binary energy (i.e. spatial context with one other region connected to pRegion) for assigning label nClassIndex to pRegion
	virtual double GetBinaryEnergy( TRegion* pRegion,
	                                TBoundary* pBoundary, int nClassIndex ) = 0;

	/// Perform annealing
	void Annealing();

// Implementation
public:
	virtual ~TMarkovNetClassifier();

};

#endif

/////////////////////////////////////////////////////////////////////////////
