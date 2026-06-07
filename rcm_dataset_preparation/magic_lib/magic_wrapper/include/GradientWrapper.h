//
// Created by max on 2021-12-09.
//

#include "ImageConversion.h"
#include "GradientCalculator.h"

#ifndef ICETOOLS_GRADIENTWRAPPER_H
#define ICETOOLS_GRADIENTWRAPPER_H

py::array_t<float> Gradient(py::array_t<float> aSrcImg)
{
    if (aSrcImg.ndim() < 2 || aSrcImg.ndim() > 3){
        throw std::runtime_error("Invalid input image!!!");
    }

    int nFeatureCnt;
    TImage** ppGradFeat;
    if (aSrcImg.ndim() == 2)
    {
        nFeatureCnt = 1;
        TGrayImage<float>* pSrcImg = PyArray2GrayImage(aSrcImg);

        ppGradFeat = new TImage * [nFeatureCnt];
        for (int i = 0; i < nFeatureCnt; i++) {
            ppGradFeat[i] = pSrcImg;
        }
    }
    else
    {
        nFeatureCnt = aSrcImg.shape(2);
        TGrayImage<float>** ppSrcImg = PyArray2GrayImageArray(aSrcImg);

        ppGradFeat = new TImage * [nFeatureCnt];
        for (int i = 0; i < nFeatureCnt; i++) {
            ppGradFeat[i] = ppSrcImg[i];
        }
    }

    double dbGradOptParam = 1;
    TGradientCalculator* pGradientCal = new TGradientCalculator(1,
                                                                0,
                                                                TGradientCalculator::GAUSSDEV,
                                                                &dbGradOptParam,
                                                                TGradientCalculator::BDVF);

    TGrayImage<float>* pGradRet = pGradientCal->Gradient(ppGradFeat, nFeatureCnt);

    for (int i = 0; i < nFeatureCnt; i++) {
        delete ppGradFeat[i];
    }
    delete[] ppGradFeat;
    delete pGradientCal;

    py::array_t<float> aGrad = GrayImage2PyArray(pGradRet);
    delete pGradRet;
    return aGrad;
}

#endif //ICETOOLS_GRADIENTWRAPPER_H
