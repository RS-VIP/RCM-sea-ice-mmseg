// GuassFilter.cpp : implementation of the TGaussFilter class
//
// Gauss filter based on Tfilter class (due to its separability)
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "GaussFilter.h"

#define Pi 3.1415926535897
#define GAUSSIAN_DIEOFF 0.0001  // for determining the bandwidth of Gaussian


//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TGaussFilter::TGaussFilter(double s)
	: TFilter()
{
	Sigma = s;
	int Width = FindWidth();
	ConstructGauss(Width);
	m_nScale = 1;
}


//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TGaussFilter::TGaussFilter(double s, int Width)
	: TFilter()
{
	Sigma = s;
	ConstructGauss(Width);
	m_nScale = 1;
}

//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TGaussFilter::TGaussFilter(int nScale, int nSearchWinSize, double dbSigmaD)
	: TFilter()
{
	Sigma = dbSigmaD;
	nSearchWinSize /= 2;
	ConstructGauss(nSearchWinSize);
	m_nScale = nScale;
}

//
// Destructor
/////////////////////////////////////////////////////////////////////////////
TGaussFilter::~TGaussFilter()
{
}


//
// Implementation of the filter cannot be of infinite length. Find where to
// truncate
// 
// Return
//  the bandwidth
/////////////////////////////////////////////////////////////////////////////
int TGaussFilter::FindWidth()
{
	double ss = 2 * Sigma * Sigma;
	int i;
	for (i = 0; i <= 30; i++)  // we assume that the maximum length is 30
	{
		if (exp(- i * i / ss) < GAUSSIAN_DIEOFF)
			break;
	}

	if (i <= 1)
		return 1;
	else
		return i - 1;
}


//
// Construct the discrete Gaussian filter
// 
// Param in
//  Width: width of the filter
/////////////////////////////////////////////////////////////////////////////
void TGaussFilter::ConstructGauss(int Width)
{
	int Length = Width * 2 + 1;
	double* Data = new double[Length];
	int i;
	for (i = 0; i <= Width; i++)
		Data[Width + i] = 0;  // Initialize with 0

	// Generate one-side Gaussian (including the origin)
	double ss = 2 * Sigma * Sigma;
	double den = 3 * sqrt(Pi * ss);  // denominator
	for (i = 0; i <= Width; i ++)
	{
		// We will average values at t-.5, t, t+.5
		Data[Width + i] += exp(- (i - 0.5) * (i - 0.5) / ss);
		Data[Width + i] += exp(- i * i / ss);
		Data[Width + i] += exp(- (i + 0.5) * (i + 0.5) / ss);
		Data[Width + i] /= den;
	}

	// Generate the other side Gaussian
	for (i = 0; i < Width; i ++)
		Data[i] = Data[Length - 1 - i];

	ConstructFilter(Data, Length, Width);
	delete[] Data;
}


//
// Obtain the discrete Gaussian derivative filter
// 
// Return
//  the filter
/////////////////////////////////////////////////////////////////////////////
TFilter* TGaussFilter::ObtainDerivative()
{
	int nWidth = m_nLength / 2;
	double* Data = new double[m_nLength];

	// Generate one-side Gaussian derivative (including the origin)
	double ss = Sigma * Sigma;
	double den = sqrt(2 * Pi) * Sigma;  // denominator
	int i;
	for (i = 0; i <= nWidth; i++)
//		Data[nWidth + i] = -i / ss * exp(- i * i / (2 * ss));
		Data[nWidth + i] = -i / ss * exp(- i * i / (2 * ss)) / den;

	// Generate the other side Gaussian derivative
	for (i = 0; i < nWidth; i ++)
		Data[i] = -Data[m_nLength - 1 - i];

	TFilter* Ret = new TFilter(Data, m_nLength);
	delete[] Data;
	return Ret;
}


//
// Obtain the second derivative of Gauss
// 
// Return
//  the filter
/////////////////////////////////////////////////////////////////////////////
TFilter* TGaussFilter::Obtain2Derivative()
{
	int nWidth = m_nLength / 2;
	double* Data = new double[m_nLength];

	// Generate one-side Gaussian derivative (including the origin)
	double ss = Sigma * Sigma;
	double den = sqrt(2 * Pi) * Sigma;  // denominator
	int i;
	for (i = 0; i <= nWidth; i ++)
		Data[nWidth + i] = (i * i / ss - 1) * exp(- i * i / (2 * ss)) / (den * ss);

	// Generate the other side Gaussian derivative
	for (i = 0; i < nWidth; i ++)
		Data[i] = Data[m_nLength - 1 - i];

	TFilter* Ret = new TFilter(Data, m_nLength);
	delete[] Data;
	return Ret;
}

TGrayImage<FLOAT>** TGaussFilter::Filtering(TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt)
{
	int i, j;
	
	TGrayImage<FLOAT>** ppRetImg = new TGrayImage<FLOAT>*[nFeatureCnt];
	for ( i = 0; i < nFeatureCnt; i++ )
		ppRetImg[i] = new TGrayImage<FLOAT>(ppSrcImg[i]);

	TGrayImage<FLOAT>** ppFiltImg;

	for ( i = 0; i < m_nScale; i++ )
	{
		ppFiltImg = Conv2(ppRetImg, nFeatureCnt);

		for ( j = 0; j < nFeatureCnt; j++ )
		{
			ppRetImg[j]->ConvertFrom(*(ppFiltImg[j]));
			delete ppFiltImg[j];
		}
		delete[] ppFiltImg;
	}

	return ppRetImg;
}
