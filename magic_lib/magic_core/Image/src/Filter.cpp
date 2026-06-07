// Filter.cpp : implementation of the TFilter class
//
// Filter operation
/////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include "Filter.h"

//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TFilter::TFilter()
{
	m_nLength = 0;
	m_nOffset = 0;
	m_dbResponse = 0;
}


// Constructor
/////////////////////////////////////////////////////////////////////////////
TFilter::TFilter(double FilterData[], int FilterLen, int FilterOff)
{
	m_dbResponse = 0;
	ConstructFilter(FilterData, FilterLen, FilterOff);
}


//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TFilter::TFilter(double FilterData[], int FilterLen)
{
	m_dbResponse = 0;
	ConstructFilter(FilterData, FilterLen, FilterLen / 2);
}


//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TFilter::TFilter(TFilter& Source)
{
	m_dbResponse = 0;
	ConstructFilter(Source.m_dbResponse, Source.m_nLength, Source.m_nOffset);
}


//
// Destructor
/////////////////////////////////////////////////////////////////////////////
TFilter::~TFilter()
{
	if (m_dbResponse)
		delete[] m_dbResponse; //MOD_QIN
}


//
// Construct the filter
/////////////////////////////////////////////////////////////////////////////
void TFilter::ConstructFilter(double FilterData[], int FilterLen,
                              int FilterOff)
{
	if (m_dbResponse)
		delete[] m_dbResponse; //MOD_QIN

	m_nLength = FilterLen;
	m_nOffset = FilterOff;
	m_dbResponse = new double[m_nLength];
	for (int i = 0; i < m_nLength; i ++)
		m_dbResponse[i] = FilterData[i];
}


//On univariate images
// 2D convolution
/////////////////////////////////////////////////////////////////////////////
TGrayImage<FLOAT>* TFilter::Conv2(TGrayImage<FLOAT>* pSrcImg)
{
	TGrayImage<FLOAT>* ColumnFiltered = Conv(pSrcImg, true);  // filtering horizontally
	TGrayImage<FLOAT>* Ret = Conv(ColumnFiltered, false);  // filtering vertically
	delete ColumnFiltered;
	return Ret;
}


// 1D convolution
// 
// Param in
//  pSrcImg: univariate source image
//  IsCol: the flag indicating if the convolution is column-wise or row-wise
//
// Return
//  the convolution result
/////////////////////////////////////////////////////////////////////////////
TGrayImage<FLOAT>* TFilter::Conv(TGrayImage<FLOAT>* pSrcImg, BOOL IsCol)
{
	int nImgWidth = pSrcImg->Width();
	int nImgHeight = pSrcImg->Height();
	TGrayImage<FLOAT>* Ret = new TGrayImage<FLOAT>(nImgWidth, nImgHeight);

	double Sum;
	int iA;
	if (IsCol)
	{
		// Column-wise
		for (int Row = 0; Row < nImgHeight; Row++)
		{
			for (int Col = 0; Col < nImgWidth; Col ++)
			{
				Sum = 0.0;
				for (int i = 0; i < m_nLength; i ++)
				{
					// We wrap the boundary data in mirror way
					iA = Map(Col - i + m_nOffset, nImgWidth, PADMIRROR);
					if (iA >= 0)
						Sum += m_dbResponse[i] * (*pSrcImg)(Row, iA);
				}
				(*Ret)(Row, Col) = (float) Sum;
			}
		}
	}
	else
	{
		// Row-wise
		for (int Col = 0; Col < nImgWidth; Col ++)
		{
			for (int Row = 0; Row < nImgHeight; Row ++)
			{
				Sum = 0.0;
				for (int i = 0; i < m_nLength; i ++)
				{
					// We wrap the boundary data in mirror way
					iA = Map(Row - i + m_nOffset, nImgHeight, PADMIRROR);
					if (iA >= 0)
						Sum += m_dbResponse[i] * (*pSrcImg)(iA, Col);
				}
				(*Ret)(Row, Col) = (float) Sum;
			}
		}
	}

	return Ret;
}


//On multivariate images
// 2D convolution when the filter is completely separable
/////////////////////////////////////////////////////////////////////////////
TGrayImage<FLOAT>** TFilter::Conv2(TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt)
{
	TGrayImage<FLOAT>** ColumnFiltered = Conv(ppSrcImg, nFeatureCnt, true);  // filtering horizontally
	TGrayImage<FLOAT>** Ret = Conv(ColumnFiltered, nFeatureCnt, false);  // filtering vertically
	for (int i = 0; i < nFeatureCnt; i++)
		if (ColumnFiltered[i]) delete ColumnFiltered[i];
	if (ColumnFiltered) delete[] ColumnFiltered;

	return Ret;
}


// 1D convolution
// 
// Param in
//  ppSrcImg: multivariate source image
//  IsCol: the flag indicating if the convolution is column-wise or row-wise
//
// Return
//  the convolved result
/////////////////////////////////////////////////////////////////////////////
TGrayImage<FLOAT>** TFilter::Conv(TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt, BOOL IsCol)
{
	int i, j;
	int nImgWidth = ppSrcImg[0]->Width();
	int nImgHeight = ppSrcImg[0]->Height();
	TGrayImage<FLOAT>** Ret = new TGrayImage<FLOAT>*[nFeatureCnt];
	for (i = 0; i < nFeatureCnt; i++)
		Ret[i] = new TGrayImage<FLOAT>(nImgWidth, nImgHeight);
	double* Sum = new double[nFeatureCnt];
	int iA;
	if (IsCol)
	{
		// Column-wise
		for (int Row = 0; Row < nImgHeight; Row ++)
		{
			for (int Col = 0; Col < nImgWidth; Col ++)
			{
				for (j = 0; j < nFeatureCnt; j++)
					Sum[j] = 0;
				for (i = 0; i < m_nLength; i ++)
				{
					// We wrap the boundary data in mirror way
					iA = Map(Col - i + m_nOffset, nImgWidth, PADMIRROR);
					if (iA >= 0)
						for (j = 0; j < nFeatureCnt; j++) 
						    Sum[j] += m_dbResponse[i] * (*ppSrcImg[j])(Row, iA);
				}
				for (j = 0; j < nFeatureCnt; j++)
				    (*Ret[j])(Row, Col) = (float) Sum[j];
			}
		}
	}
	else
	{
		// Row-wise
		for (int Col = 0; Col < nImgWidth; Col ++)
		{
			for (int Row = 0; Row < nImgHeight; Row ++)
			{
				for (j = 0; j < nFeatureCnt; j++)
					Sum[j] = 0;
				for (i = 0; i < m_nLength; i ++)
				{
					// We wrap the boundary data in mirror way
					iA = Map(Row - i + m_nOffset, nImgHeight, PADMIRROR);
					if (iA >= 0)
						for (j = 0; j < nFeatureCnt; j++) 
						    Sum[j] += m_dbResponse[i] * (*ppSrcImg[j])(iA, Col);
				}
				for (j = 0; j < nFeatureCnt; j++)
				    (*Ret[j])(Row, Col) = (float) Sum[j];
			}
		}
	}
    delete[] Sum;
	return Ret;
}


//
// Map the extended signal back into the original range for bordering problems
// in filtering
// 
// Param in
//  x: position
//  Range: original range of the signal
//
// Return
//  mapped position
/////////////////////////////////////////////////////////////////////////////
int TFilter::Map(int x, int Range, int Mode)
{
	switch (Mode)
	{
	case PADMIRROR: //MOD_QIN
		while (x < 0 || x >= Range)
			if (x < 0) 
				x = -x - 1;
			else if (x >= Range)
				x = Range - 1 - (x - Range);
		return x;
	case PADZERO:
		if (x < 0 || x >= Range)
			return 0;
		else
			return x;
	default:
		return -1;
	}
}


//
// Scale the filter
/////////////////////////////////////////////////////////////////////////////
void TFilter::operator*=(FLOAT Value)
{
	for (int i = 0; i < m_nLength; i ++)
		m_dbResponse[i] *= Value;
}


//
// Shift the filter left
/////////////////////////////////////////////////////////////////////////////
void TFilter::operator<<(int Delta)
{
	if (Delta > 0)
		m_nOffset += Delta;
	if (m_nOffset >= m_nLength) //MOD_QIN
		m_nOffset = m_nLength - 1;
}


//
// Shift the filter right
/////////////////////////////////////////////////////////////////////////////
void TFilter::operator>>(int Delta)
{
	if (Delta > 0)
		m_nOffset -= Delta;
	if (m_nOffset < 0) //MOD_QIN
		m_nOffset = 0;
}


//
// Get the norm of the filter (element summation)
/////////////////////////////////////////////////////////////////////////////
double TFilter::GetNorm()
{
	double Norm = 0;
	for (int i = 0; i < m_nLength; i ++)
		Norm += m_dbResponse[i];
	return Norm;
}
