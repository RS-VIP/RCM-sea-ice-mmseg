//
// Created by max on 2021-07-05.
// ImageConversion.h: functions for converting TGrayImages to numpy arrays for the MAGIC python interface.
//

#ifndef MAGIC_PORT_IMAGECONVERSION_H
#define MAGIC_PORT_IMAGECONVERSION_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py=pybind11;

#include "GrayImage.h"
#include "Vector.h"

/// convert a py::array to a univariate TGrayImage. Caller is responsible for cleaning up memory!!!
template <typename T> TGrayImage<T>* PyArray2GrayImage(py::array_t<T> pyArray){
    int H = pyArray.shape(0), W=pyArray.shape(1);
    TGrayImage<T>* pGrayImg = new TGrayImage<T>(W,H);

    auto ra = pyArray.template unchecked<2>();

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            (*pGrayImg)(i,j) = ra(i,j);
        }
    }
    return pGrayImg;
}

/// convert a py::array to an array of TGrayImages (i.e., a multivariate image). Caller is responsible for cleaning up memory!!!
template <typename T> TGrayImage<T>** PyArray2GrayImageArray(py::array_t<T> pyArray){
    int H = pyArray.shape(0), W=pyArray.shape(1);
    int C = pyArray.ndim() == 3 ? pyArray.shape(2): 1;
    TGrayImage<T>** ppGrayImg = new TGrayImage<T>*[C];
    auto ra = pyArray.template unchecked<3>();

    for (int k=0; k<C; k++){
        ppGrayImg[k] = new TGrayImage<T>(W,H);

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                (*ppGrayImg[k])(i,j) = ra(i,j,k);
            }
        }
    }

    return ppGrayImg;
}

/// convert a univariate TGrayImage to a py::array.
template <typename T> py::array_t<T> GrayImage2PyArray(TGrayImage<T>* pGrayImg){
    int H = pGrayImg->Height(), W = pGrayImg->Width();

    size_t arr_size = H*W;
    T* buf = new T[arr_size];

    // Create a Python object that will free the allocated
    // memory when destroyed:
    py::capsule free_when_done(buf, [](void *f) {
        T *buf = reinterpret_cast<T *>(f);
        std::cerr << "freeing memory @ " << f << "\n";
        delete[] buf;
    });

    py::array_t<T, py::array::c_style> pyArray({H,W}, buf, free_when_done);

    auto ra = pyArray.template mutable_unchecked<2>();

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            ra(i,j) = pGrayImg->GetPixelData(j,i);
        }
    }
    return pyArray;
}

/// convert an array of TGrayImages (i.e., a multivariate image) to a py::array.
template <typename T> py::array_t<T> GrayImageArray2PyArray(TGrayImage<T>** ppGrayImg, int nChan){
    int H = ppGrayImg[0]->Height(), W = ppGrayImg[0]->Width();

    size_t arr_size = H*W*nChan;
    T* buf = new T[arr_size];
    //    size_t tsz = sizeof(T);

    // Create a Python object that will free the allocated
    // memory when destroyed:
    py::capsule free_when_done(buf, [](void *f) {
        T *buf = reinterpret_cast<T *>(f);
        std::cerr << "freeing memory @ " << f << "\n";
        delete[] buf;
    });

    py::array_t<T, py::array::c_style> pyArray({H,W, nChan},
                                               buf,
                                               free_when_done);

    auto ra = pyArray.template mutable_unchecked<3>();

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            for (int k=0; k<nChan; k++)
                ra(i,j,k) = (ppGrayImg[k])->GetPixelData(j,i);
        }
    }
    return pyArray;
}

template <typename T> py::array_t<T> VectorVector2PyArray(std::vector<TVector*> pVecVec){
    int H = pVecVec.size();
    int W = pVecVec[0]->GetElementCnt();
    py::array_t<T, py::array::c_style> pyArray({H,W});

    auto ra = pyArray.template mutable_unchecked<2>();

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            ra(i,j) = (T)(*pVecVec[i])(j);
        }
    }
    return pyArray;
};

template <typename T>  std::vector<TVector*> PyArray2VectorVector(py::array_t<T> pyArray){
    int H = pyArray.shape(0), W=pyArray.shape(1);

    std::vector<TVector*> vRet;

    auto ra = pyArray.template unchecked<2>();

    for (int i = 0; i < H; i++) {
        TVector* pVector = new TVector(W, true);
        for (int j = 0; j < W; j++) {
            (*pVector)(j) = (double) ra(i,j);
        }
        vRet.push_back(pVector);
    }
    return vRet;
};

#endif //MAGIC_PORT_IMAGECONVERSION_H
