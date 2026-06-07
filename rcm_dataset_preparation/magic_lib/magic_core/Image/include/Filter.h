/** \file
 *
 * \brief TFilter class
 *
 * \author Qiyao Yu
 *
 * \date 2000
 *
 */

#ifndef TFilter_Class
#define TFilter_Class


#include "GrayImage.h"


static double SQRT3 = 1.73205080756887729352745;
static double MM_SQRT2 = 1.41421356237309504880168872420969808;


static double cHLeGall[] = {
	-1 / 8.0,
	 2 / 8.0,
	 6 / 8.0,
	 2 / 8.0,
	-1 / 8.0
};
static double gHLeGall[] = {
	 1 / 2.0,
	 2 / 2.0,
	 1 / 2.0
};


static double cHHaar[] = {
	MM_SQRT2 * 1.0 / 2.0,
	MM_SQRT2 * 1.0 / 2.0
};


static double cHDaubechies_4[] = {
	MM_SQRT2 * (1 + SQRT3) / 8.0,
	MM_SQRT2 * (3 + SQRT3) / 8.0,
	MM_SQRT2 * (3 - SQRT3) / 8.0,
	MM_SQRT2 * (1 - SQRT3) / 8.0
};


static double cHDaubechies_6[] = {
	 0.3327,
	 0.8069,
	 0.4599,
	-0.1350,
	-0.0854,
	 0.0352
};


static double cHDaubechies_10[] = {
	 0.1601,
	 0.6083,
	 0.7243,
	 0.1384,
	-0.2423,
	-0.0322,
	 0.0776,
	-0.0062,
	-0.0126,
	 0.0033
};


static double cHDaubechies_12[] = {
	 0.1115407433501095,
	 0.4946238903984533,
	 0.7511339080210959,
	 0.3152503517091982,
	-0.2262646939654400,
	-0.1297668675672625,
	 0.0975016055873225,
	 0.0275228655303053,
	-0.0315820393184862,
	 0.0005538422011614,
	 0.0047772575119455,
	-0.0010773010853085
};


static double cHDaubechies_20[] = {
	 0.026670057901,
	 0.188176800078,
	 0.527201188932,
	 0.688459039454,
	 0.281172343661,
	-0.249846424327,
	-0.195946274377,
	 0.127369340336,
	 0.093057364604,
	-0.071394147166,
	-0.029457536822,
	 0.033212674059,
	 0.003606553567,
	-0.010733175483,
	 0.001395351747,
	 0.001992405295,
	-0.000685856695,
	-0.000116466855,
	 0.000093588670,
	-0.000013264203
};


static double cHBurtAdelson[] = {
	MM_SQRT2 * -1.0 / 20.0,
	MM_SQRT2 *  5.0 / 20.0,
	MM_SQRT2 * 12.0 / 20.0,
	MM_SQRT2 *  5.0 / 20.0,
	MM_SQRT2 * -1.0 / 20.0
};
static double gHBurtAdelson[] = {
	MM_SQRT2 *  -3.0 / 280.0,
	MM_SQRT2 * -15.0 / 280.0,
	MM_SQRT2 *  73.0 / 280.0,
	MM_SQRT2 * 170.0 / 280.0,
	MM_SQRT2 *  73.0 / 280.0,
	MM_SQRT2 * -15.0 / 280.0,
	MM_SQRT2 *  -3.0 / 280.0
};


static double cHAntonini[] = {
	 0.026749,
	-0.016864,
	-0.078223,
	 0.26686,
	 0.60295,
	 0.26686,
	-0.078223,
	-0.016864,
	 0.026749
};
static double gHAntonini[] = {
	-0.091272,
	-0.057544,
	 0.59128,
	 1.11508,
	 0.59128,
	-0.057544,
	-0.091272
};

/// Filter operation
class TFilter
{
public:
	TFilter();
	TFilter(double FilterData[], int FilterLen, int FilterOff);
	TFilter(double FilterData[], int FilterLen);
	TFilter(TFilter& Source);

// Attributes
public:
	double* m_dbResponse;
	int m_nLength;
	int m_nOffset;

// Operations
public:
	/// On univariate image
	TGrayImage<FLOAT>* Conv2(TGrayImage<FLOAT>* pSrcImg);
	/// On univariate image
	TGrayImage<FLOAT>* Conv(TGrayImage<FLOAT>* pSrcImg, BOOL IsCol);

	/// On multivariate image
	TGrayImage<FLOAT>** Conv2(TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt);
	/// On multivariate image
	TGrayImage<FLOAT>** Conv(TGrayImage<FLOAT>** ppSrcImg,int nFeatureCnt, BOOL IsCol);

	//TFilter* Conv(TFilter* fltr);  //This program is wrong!! QIN_20090604

	/// Map the extended signal back into the original range for bordering problems in filtering
	int Map(int x, int Range, int Mode);

	void operator*=(FLOAT Value);
	void operator<<(int Delta);
	void operator>>(int Delta);

	double GetNorm();

protected:
	void ConstructFilter(double FilterData[], int FilterLen, int FilterOff);

// Implementation
public:

	virtual ~TFilter();

};

#endif
