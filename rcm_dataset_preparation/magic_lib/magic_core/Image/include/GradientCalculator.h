/** \file
 *
 * \brief TGradientCalculator class.
 *
 */

#ifndef TGradientCalculator_Class
#define TGradientCalculator_Class

#include "GrayImage.h"
#include "MonoImage.h"
#include "GaussFilter.h"

/// Gradient Calculator
class TGradientCalculator
{
public:
	TGradientCalculator(TMonoImage* pMask);
	TGradientCalculator(BOOL NormalizeFlag, TMonoImage* pMask);
	TGradientCalculator(BOOL NormalizeFlag, TMonoImage* pMask, int nGradOperator, double* dbGradOptParam);
	TGradientCalculator(BOOL NormalizeFlag, TMonoImage* pMask, int nGradOperator, double* dbGradOptParam, int nMultiGradOpt);

// Attributes
public:
	// Constant
	enum GradOperatorType {GAUSSDEV = 0, SOBEL = 1, PREWITT = 2};
	enum MultiGradType {BDVF = 0, L2NORM = 1};

	GradOperatorType m_nGradOperator;
	MultiGradType m_nMultiGradOpt;

	double* m_dbGradOptParam;
	BOOL m_bNormalized;

	TMonoImage* m_pMask;

// Operations
public:
	TGrayImage<FLOAT>* Gradient(TImage** ppSrcImg, int nFeatureCnt);
	TGrayImage<FLOAT>* Gradient(TGrayImage<FLOAT>** ppHDiff, TGrayImage<FLOAT>** ppVDiff, TImage** ppSrcImg, int nFeatureCnt);
	TGrayImage<FLOAT>* GradientFilter(TGrayImage<FLOAT>* pSrcImg, char cDirection);

// Implementation
public:
	virtual ~TGradientCalculator();
};

#endif

