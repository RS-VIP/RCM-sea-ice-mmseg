/** \file
 *
 * \brief TGaussFilter class.
 *
 */

#ifndef TGaussFilter_Class
#define TGaussFilter_Class


#include "Filter.h"
#include "GrayImage.h"

/// Gauss filter based on TFilter class (due to its separability)
class TGaussFilter : public TFilter
{
// Constructor
public:
	TGaussFilter(double s);
	TGaussFilter(double s, int Width);
	TGaussFilter(int nScale, int nSearchWinSize, double dbSigmaD);

// Data member
protected:
	double Sigma;
	int m_nScale;

// Member function
public:
	TFilter* ObtainDerivative();
	TFilter* Obtain2Derivative();
	TGrayImage<FLOAT>** Filtering(TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt);

private:
	int FindWidth();
	void ConstructGauss(int Width);

// Constructor
public:
	virtual ~TGaussFilter();

};

#endif
