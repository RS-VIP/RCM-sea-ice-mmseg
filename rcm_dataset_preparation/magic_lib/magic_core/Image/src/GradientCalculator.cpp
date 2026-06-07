// GradientCalculator.h : implementation of the TGradientCalculator class
//
// Gradient Calculator
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "GradientCalculator.h"

double xSOBEL2D[3][3] = {
		{1/4, 0, -1/4},
		{2/4, 0, -2/4},
		{1/4, 0, -1/4}
	};
double ySOBEL2D[3][3] = {
		{1/4, 2/4, 1/4},
		{0, 0, 0},
		{-1/4, -2/4, -1/4}
	};

double xPREWITT2D[3][3] = {
		{1/3, 0, -1/3},
		{1/3, 0, -1/3},
		{1/3, 0, -1/3}
	};
double yPREWITT2D[3][3] = {
		{1/3, 1/3, 1/3},
		{0, 0, 0},
		{-1/3, -1/3, -1/3}
	};

double xSOBEL2Dx[3] = {1/2, 0, -1/2};
double xSOBEL2Dy[3] = {1/2, 2/2, 1/2};
double ySOBEL2Dx[3] = {1/2, 2/2, 1/2};
double ySOBEL2Dy[3] = {1/2, 0, -1/2};
	
double xPREWITT2Dx[3] = {1/SQRT3, 0, -1/SQRT3};
double xPREWITT2Dy[3] = {1/SQRT3, 1/SQRT3, 1/SQRT3};
double yPREWITT2Dx[3] = {1/SQRT3, 1/SQRT3, 1/SQRT3};
double yPREWITT2Dy[3] = {1/SQRT3, 0, -1/SQRT3};

//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TGradientCalculator::TGradientCalculator(TMonoImage* pMask)
{
	m_bNormalized = 1;
	m_nGradOperator = GAUSSDEV;
	*m_dbGradOptParam = 1; 
	m_nMultiGradOpt = BDVF;
	m_pMask = pMask;
}

//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TGradientCalculator::TGradientCalculator(BOOL NormalizeFlag, TMonoImage* pMask)
{
	m_bNormalized = NormalizeFlag;
	m_nGradOperator = GAUSSDEV;
	*m_dbGradOptParam = 1; 
	m_nMultiGradOpt = BDVF;
	m_pMask = pMask;
}

//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TGradientCalculator::TGradientCalculator(BOOL NormalizeFlag, TMonoImage* pMask, int nGradOperator, double* dbGradOptParam)
{
	m_bNormalized = NormalizeFlag;
	m_nGradOperator = (GradOperatorType) nGradOperator;
	m_dbGradOptParam = dbGradOptParam;
	m_nMultiGradOpt = BDVF;
	m_pMask = pMask;
}

//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TGradientCalculator::TGradientCalculator(BOOL NormalizeFlag, TMonoImage* pMask, int nGradOperator, double* dbGradOptParam, int nMultiGradOpt)
{
	m_bNormalized = NormalizeFlag;
	m_nGradOperator = (GradOperatorType) nGradOperator;
	m_dbGradOptParam = dbGradOptParam;
	m_nMultiGradOpt = (MultiGradType) nMultiGradOpt;
	m_pMask = pMask;
}

//
// Destructor
/////////////////////////////////////////////////////////////////////////////
TGradientCalculator::~TGradientCalculator()
{
}

/*
// Pre-process the multivariate image on per image base
TGrayImage<FLOAT>** TGradientCalculator::MultiImgPreProc(TImage** ppImgSrc, int nFeatCnt, int* nFilterType,double** dbParamSet)
{
	TGrayImage<FLOAT>** pRet = new TGrayImage<FLOAT>*[nFeatCnt];
	for (int k = 0; k < nFeatCnt; k++)
	{   
		TGrayImage<FLOAT>* pImg = new TGrayImage<FLOAT>(ppImgSrc[k]);
		switch (nFilterType[k])
		{
			case GAUSSIAN: TGaussFilter* pGauss = new TGaussFilter(dbParamSet[k][0], 0, 0);
                           pRet[k] = pGauss->Conv2(pImg);
						   delete pImg;
	                       delete pGauss;
						   break;
			case BILATERAL: TBilatFilter* pBilat = new TBilatFilter(dbParamSet[k][0],dbParamSet[k][1],(int) dbParamSet[k][2]);
			                pRet[k] = pBilat->Filtering(pImg);
							delete pImg;
							delete pBilat;
							break;
			default: pRet[k] = pImg;
		}
	}
	return pRet;
}
*/

// Gradient filter
// Param in: 
// pSrcImg: univariate source image
// cGradComp: gradient component
TGrayImage<FLOAT>* TGradientCalculator::GradientFilter(TGrayImage<FLOAT>* pSrcImg, char cGradComp)
{   
	switch (m_nGradOperator)
	{
		case GAUSSDEV: switch (cGradComp)
					   {
							case 'X': {
			                          TGaussFilter* pGaussFilter = new TGaussFilter(*m_dbGradOptParam);
			                          TFilter* pGaussDevFilter = pGaussFilter->ObtainDerivative();
								      TGrayImage<FLOAT>* pAveRet = pGaussDevFilter->Conv(pSrcImg, true);
					                  TGrayImage<FLOAT>* pRet = pGaussFilter->Conv(pAveRet, false);
									  delete pGaussFilter;
					                  delete pGaussDevFilter;
									  delete pAveRet;
									  return pRet;
									  }
							case 'Y': {
								      TGaussFilter* pGaussFilter = new TGaussFilter(*m_dbGradOptParam);
			                          TFilter* pGaussDevFilter = pGaussFilter->ObtainDerivative();
					                  TGrayImage<FLOAT>* pAveRet = pGaussDevFilter->Conv(pSrcImg, false);
								      TGrayImage<FLOAT>* pRet = pGaussFilter->Conv(pAveRet, true);
									  delete pGaussFilter;
					                  delete pGaussDevFilter;
									  delete pAveRet;
									  return pRet;
									  }
							default: {
								      cout << "Error: No such gradient component!" << endl;
								      exit(1);
									 }
					   }
					   break;
		case SOBEL: switch (cGradComp)
					{
						case 'X': {
							      TFilter* pxFilter = new TFilter(xSOBEL2Dx, 3);
							      TFilter* pyFilter = new TFilter(xSOBEL2Dy, 3);
							      TGrayImage<FLOAT>* pColRet = pxFilter->Conv(pSrcImg, true);
					              TGrayImage<FLOAT>* pRet = pyFilter->Conv(pColRet, false);
								  delete pxFilter;
					              delete pyFilter;
					              delete pColRet;
								  return pRet;
								  }
						case 'Y': {
							      TFilter* pxFilter = new TFilter(ySOBEL2Dx, 3);
							      TFilter* pyFilter = new TFilter(ySOBEL2Dy, 3);
							      TGrayImage<FLOAT>* pColRet = pxFilter->Conv(pSrcImg, true);
					              TGrayImage<FLOAT>* pRet = pyFilter->Conv(pColRet, false);
								  delete pxFilter;
					              delete pyFilter;
					              delete pColRet;
								  return pRet;
								  }
						default: {
							     cout << "Error: No such gradient component!" << endl;
								 exit(1);
								 }
					 }
					 break;
		case PREWITT: switch (cGradComp)
					{
		                case 'X': {
							      TFilter* pxFilter = new TFilter(xPREWITT2Dx, 3);
							      TFilter* pyFilter = new TFilter(xPREWITT2Dy, 3);
							      TGrayImage<FLOAT>* pColRet = pxFilter->Conv(pSrcImg, true);
					              TGrayImage<FLOAT>* pRet = pyFilter->Conv(pColRet, false);
								  delete pxFilter;
					              delete pyFilter;
								  delete pColRet;
								  return pRet;
								  }
						case 'Y': {
							      TFilter* pxFilter = new TFilter(yPREWITT2Dx, 3);
							      TFilter* pyFilter = new TFilter(yPREWITT2Dy, 3);
							      TGrayImage<FLOAT>* pColRet = pxFilter->Conv(pSrcImg, true);
					              TGrayImage<FLOAT>* pRet = pyFilter->Conv(pColRet, false);
								  delete pxFilter;
					              delete pyFilter;
								  delete pColRet;
								  return pRet;
								  }
						default: {
							     cout << "Error: No such gradient component!" << endl;
								 exit(1);
								 }
					 }
					 break;
		default: cout << "Error: No such gradient operator!" << endl;
			     exit(1);
	}
}

//
// Calculate the gradient magnitude on multivariate images
/////////////////////////////////////////////////////////////////////////////
TGrayImage<FLOAT>* TGradientCalculator::Gradient(TImage** ppSrcImg, int nFeatureCnt)
{
	TGrayImage<FLOAT>** ppDiffH = new TGrayImage<FLOAT>*[nFeatureCnt];
	TGrayImage<FLOAT>** ppDiffV = new TGrayImage<FLOAT>*[nFeatureCnt];
	TGrayImage<FLOAT>* pRet = Gradient(ppDiffH, ppDiffV, ppSrcImg, nFeatureCnt);
    if (pRet == NULL)
		return NULL;
	for (int k = 0; k < nFeatureCnt; k++) 
	{
		if (ppDiffH[k]) delete ppDiffH[k];
	    if (ppDiffV[k]) delete ppDiffV[k];
	}
	if (ppDiffH) delete[] ppDiffH;
	if (ppDiffV) delete[] ppDiffV;

	return pRet;
}


//
// Calculate the gradient on multivariate images
//  ppDiffH: the horizontal difference image
//  ppDiffV: the vertical difference image
/////////////////////////////////////////////////////////////////////////////
TGrayImage<FLOAT>* TGradientCalculator::Gradient(TGrayImage<FLOAT>** ppDiffH, TGrayImage<FLOAT>** ppDiffV,
                                                TImage** ppSrcImg, int nFeatureCnt)
{
	if (nFeatureCnt <= 0)
	{
		cout << "Error: The image dimension should be above 0!" << endl;
	    exit(1);
	}
	else if (nFeatureCnt == 1)
	{
		TGrayImage<FLOAT>* pImg = new TGrayImage<FLOAT>(ppSrcImg[0]);
	    int nImgWidth = pImg->Width();
	    int nImgHeight = pImg->Height();
	    ppDiffH[0] = GradientFilter(pImg,'X');  // horizontal difference
		ppDiffV[0] = GradientFilter(pImg,'Y');  // vertical difference

		// Compute the magnitude
	    float fMaxVal = 0; // maximum gradient value within the mask or whole image (if no mask).
		float fTotalMaxVal = 0; // maximum gradient value in the whole image.
	    int i, j;
		for (i = 0; i < nImgHeight; i ++)
			for (j = 0; j < nImgWidth; j ++)
			{
				(*pImg)(i, j) = (FLOAT) sqrt((*ppDiffH[0])(i, j) * (*ppDiffH[0])(i, j) +
		                             (*ppDiffV[0])(i, j) * (*ppDiffV[0])(i, j));
					
				if ( fTotalMaxVal < (*pImg)(i, j) ) 
					fTotalMaxVal = (*pImg)(i, j);

				// calculate maximum gradient value within the mask or whole image (if no mask).
				if ( !m_pMask || ( m_pMask && (*m_pMask)(i, j) ) )
				{
					if ( fMaxVal < (*pImg)(i, j) ) 
						fMaxVal = (*pImg)(i, j);
				}			
			}
		if (fMaxVal == 0 || fTotalMaxVal == 0)
			return NULL;
		if (m_bNormalized)
		{
			for (i = 0; i < nImgHeight; i ++)
				for (int j = 0; j < nImgWidth; j ++)
				{
					if ( m_pMask && !(*m_pMask)(i, j) )
					{
						// this ensures that all gradient values (even those outside mask) are < 1.
						(*pImg)(i, j) = (*pImg)(i, j) / fTotalMaxVal;
					}
					else
						(*pImg)(i, j) = (*pImg)(i, j) / fMaxVal;
				}
		}
		return pImg;
	}
	else
	{   
		int k;
		for (k = 0; k < nFeatureCnt; k++)   // MK  consider for bug 
		{
			TGrayImage<FLOAT>* pImg = new TGrayImage<FLOAT>(ppSrcImg[k]);
	        ppDiffH[k] = GradientFilter(pImg,'X');  // horizontal difference
		    ppDiffV[k] = GradientFilter(pImg,'Y');  // vertical difference
			delete pImg;
		}
         
        
	    // Compute the magnitude
		TGrayImage<FLOAT>* pImg = new TGrayImage<FLOAT>(ppSrcImg[0]);
		int nImgWidth = pImg->Width(); // An IMPORTANT assumption: all component images have the same size COM_QIN
	    int nImgHeight = pImg->Height();
	    int i, j;
	    float fMaxVal = 0;
		float fTotalMaxVal = 0;
	    for (i = 0; i < nImgHeight; i ++)
			for (j = 0; j < nImgWidth; j ++)
			{
				switch (m_nMultiGradOpt)
				{
					case BDVF: {
						       float dx, dy, p = 0, t = 0, q = 0;
					           for (k = 0; k < nFeatureCnt; k++)
							   {
								   dx = (*ppDiffH[k])(i, j);
						           dy = (*ppDiffV[k])(i, j);
						           p += dx * dx;
						           t += dx * dy;
						           q += dy * dy;
							   }
					           (*pImg)(i, j) = (float) sqrt (0.5 * (p + q + sqrt((p + q) * (p + q) - 4 * (p * q - t * t))));
							   }
							   break;
					case L2NORM: {
						         float sum = 0;
							     for (int k = 0; k < nFeatureCnt; k++)
							     sum += (*ppDiffH[k])(i, j) * (*ppDiffH[k])(i, j) + (*ppDiffV[k])(i, j) * (*ppDiffV[k])(i, j);
                                 (*pImg)(i, j) = (float) sqrt(sum);
								 }
								 break;
					default: {
						     cout << "ERROR: No such multi-dimensional gradient computation method!" << endl;
							 exit(1);
							 }
				}

				if ( fTotalMaxVal < (*pImg)(i, j) ) 
					fTotalMaxVal = (*pImg)(i, j);

				if ( !m_pMask || ( m_pMask && (*m_pMask)(i, j) ) )
				{
					if ( fMaxVal < (*pImg)(i, j) ) 
						fMaxVal = (*pImg)(i, j);
				}
			}
			if (m_bNormalized)
			{
				for (i = 0; i < nImgHeight; i ++)
					for (int j = 0; j < nImgWidth; j ++)
					{
						if ( m_pMask && !(*m_pMask)(i, j) )
							(*pImg)(i, j) = (*pImg)(i, j) / fTotalMaxVal;
						else
							(*pImg)(i, j) = (*pImg)(i, j) / fMaxVal;
					}
			}
			return pImg;
	}
}
