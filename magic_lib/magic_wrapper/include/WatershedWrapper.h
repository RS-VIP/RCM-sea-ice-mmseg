//
// Created by max on 2021-12-09.
//

#include "ImageConversion.h"
#include "Watershed.h"

#ifndef ICETOOLS_WATERSHEDWRAPPER_H
#define ICETOOLS_WATERSHEDWRAPPER_H

py::array_t<int> GradientWatershed(py::array_t<float> aSrcImg, int nFeatureCnt)
{
    TGrayImage<float>** ppSrcImg = PyArray2GrayImageArray(aSrcImg);

    TImage** ppGradFeat = new TImage * [nFeatureCnt];
    for (int i = 0; i < nFeatureCnt; i++) {
        ppGradFeat[i] = ppSrcImg[i];
    }

    double dbGradOptParam = 1;
    TGradientCalculator* pGradientCal = new TGradientCalculator(1,
                                                                0,
                                                                TGradientCalculator::GAUSSDEV,
                                                                &dbGradOptParam,
                                                                TGradientCalculator::BDVF);

    TGrayImage<float>* pGradRet = pGradientCal->Gradient(ppGradFeat, nFeatureCnt);

    delete[] ppGradFeat;
    delete pGradientCal;

    TGrayImage<FLOAT>* pFloatGradient;
    pFloatGradient = new TGrayImage<FLOAT>(pGradRet);
    (*pFloatGradient) *= 255;


    TGrayImage<int>* pIntGradient = new TGrayImage<int>(pFloatGradient);	// for the subsequent watershed seg.
    TWaterShed* pWaterShed = new TWaterShed(TWaterShed::NEIGHBOR8);

    TGrayImage<int>* pMap = pWaterShed->Segment(pIntGradient, 0, (TFilter*)NULL);

    delete pFloatGradient;
    delete pIntGradient;
    delete pWaterShed;

    py::array_t<int> aMap = GrayImage2PyArray(pMap);
    return aMap;
}

py::array_t<int> Watershed(py::array_t<float> aGradMag, py::array_t<int> Mask)
{
    TGrayImage<float>* pFloatGradient = PyArray2GrayImage(aGradMag);
    (*pFloatGradient) *= 255;

    TGrayImage<int>* pLandmask = PyArray2GrayImage(Mask);
    TMonoImage* pMask = new TMonoImage(pLandmask);
    delete pLandmask;

    TGrayImage<int>* pIntGradient = new TGrayImage<int>(pFloatGradient);	// for the subsequent watershed seg.
    TWaterShed* pWaterShed = new TWaterShed(TWaterShed::NEIGHBOR8);

    TGrayImage<int>* pMap = pWaterShed->Segment(pIntGradient, pMask, (TFilter*)NULL);

    delete pFloatGradient;
    delete pIntGradient;
    delete pWaterShed;
    delete pMask;

    py::array_t<int> aMap = GrayImage2PyArray(pMap);
    return aMap;
}

py::array_t<float> BlockAverage(py::array_t<float> aSrcImg, int nFeatureCnt, int nBlkSize)
{
    TGrayImage<float>** ppSrcImg = PyArray2GrayImageArray(aSrcImg);

    TGrayImage<float>** ppRet = new TGrayImage<float> * [nFeatureCnt];
    for (int i = 0; i < nFeatureCnt; i++) {
        ppRet[i] = ppSrcImg[i];
        ppRet[i]->BlockAverage(nBlkSize);
    }
    py::array_t<float> aGrad = GrayImageArray2PyArray(ppRet, nFeatureCnt);
    for (int i = 0; i < nFeatureCnt; i++) {
        delete ppRet[i];
    }
    delete[] ppRet;
    return aGrad;
}


#endif //ICETOOLS_WATERSHEDWRAPPER_H
