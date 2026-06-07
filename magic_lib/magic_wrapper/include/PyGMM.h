//
// Created by max on 2021-08-16.
//

#ifndef MAGIC_PORT_PYGMM_H
#define MAGIC_PORT_PYGMM_H

#include <tuple>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "Vector.h"
#include "GaussianRegionClustering.h"


namespace py=pybind11;

class GMM{
public:

    // default constructor
    GMM(int nClassCnt, int nFeatureCnt);
    ~GMM();

    int m_nClassCnt;
    int m_nFeatureCnt;

    /// means of the clusters
    TVector** m_ppMu;
    /// covariance inversion of the clusters
    TMatrix** m_ppCovInv;
    /// prior probs of the clusters
    double* m_dbPrior;
    /// determinant of covariance
    double* m_dbCovDet;
    /// Logs of standard deviation and priors, for speed consideration
    double* m_dbLnSigma;
    /// Logs of standard deviation and priors, for speed consideration
    double* m_dbLnPrior;

    // slope of linear relationship between mean and incidence angle.
    TVector** m_ppIncSlope;

    void CopyFromGaussClassifier(TGaussianRegionClustering* pClassifier);
};

class PyGMM{
    /** a class for sending GMM parameters to python */
public:
    PyGMM(int nClass,
          int nFeat,
          py::array_t<float> Mu,
          py::array_t<float> CovInv,
          py::array_t<float> Prior,
          py::array_t<float> CovDet,
          py::array_t<float> LnSigma,
          py::array_t<float> LnPrior,
          py::array_t<float> IncSlope); // constructor from numpy arrays

    PyGMM(GMM* gmm); // constructor from gmm class

    PyGMM(py::dict GMMParams);

    GMM* togmm(); // can't be accessed through pybind interface bc gmm not a python type
    py::dict todict();

    int m_nClass;
    int m_nFeat;

    py::array_t<float> m_Mu;
    py::array_t<float> m_CovInv;
    py::array_t<float> m_Prior;
    py::array_t<float> m_CovDet;
    py::array_t<float> m_LnSigma;
    py::array_t<float> m_LnPrior;
    py::array_t<float> m_IncSlope;

    py::array_t<float> GetMu(){ return m_Mu; }
    py::array_t<float> GetCovInv(){ return m_CovInv; }
    py::array_t<float> GetPrior(){ return m_Prior; }
    py::array_t<float> GetCovDet(){ return m_CovDet; }
    py::array_t<float> GetLnSigma(){ return m_LnSigma; }
    py::array_t<float> GetLnPrior(){ return m_LnPrior; }
    py::array_t<float> GetIncSlope(){ return m_IncSlope; }
};

#endif //MAGIC_PORT_PYGMM_H
