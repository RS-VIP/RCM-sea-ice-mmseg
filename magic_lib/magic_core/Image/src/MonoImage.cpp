// MonoImage.cpp : implementation of the TMonoImage class
//
// image of 2 intensity levels
/////////////////////////////////////////////////////////////////////////////

#include "MonoImage.h"


//
// Constructor
//
TMonoImage::TMonoImage()
	: TGrayImage<BYTE>()
{
}

//
// Constructor
//
TMonoImage::TMonoImage(int nWidth, int nHeight)
	: TGrayImage<BYTE>(nWidth, nHeight)
{
}

//
// Constructor
//
TMonoImage::TMonoImage(TImage* pImg)
	: TGrayImage<BYTE>(pImg)
{

	// threshold and assign the two levels to 0 and 255
	double dbMax, dbMin;
	TImage::GetIntensityRange(dbMax, dbMin, 0);
	double dbMid = (dbMax + dbMin) / 2;
	for (int i = 0; i < m_nHeight; i ++)
	{
		for (int j = 0; j < m_nWidth; j ++)
		{
			if ((*this)(i, j) < dbMid)
				m_ppValues[i][j] = LOW;
			else
				m_ppValues[i][j] = HIGH;
		}
	}
}
//
// Construct a mono image based on a user-defined threshold.
//
TMonoImage::TMonoImage(TImage* pImg, int threshold)
	: TGrayImage<BYTE>(pImg)
{
	for (int i = 0; i < m_nHeight; i++)
	{
		for (int j = 0; j < m_nWidth; j++)
		{
			if((*this)(i, j) < threshold){
				m_ppValues[i][j] = LOW;
			} else {
				m_ppValues[i][j] = HIGH;
			}
		}
	}
}


//
// Destructor
//
TMonoImage::~TMonoImage()
{
}
