//
// Created by max on 2021-08-16.
//

#include "PyGMM.h"

/* GMM IMPLEMENTATION */

GMM::GMM(int nClassCnt, int nFeatureCnt){
    m_nClassCnt = nClassCnt;
    m_nFeatureCnt = nFeatureCnt;
    m_ppMu = new TVector*[nClassCnt];
    m_ppCovInv = new TMatrix*[nClassCnt];
    m_ppIncSlope = new TVector*[nClassCnt];
    for ( int i = 0; i < nClassCnt; i++ )
    {
        m_ppMu[i] = new TVector(nFeatureCnt, false);
        m_ppMu[i]->SetZero();
        m_ppIncSlope[i]= new TVector(nFeatureCnt, false);
        m_ppIncSlope[i]->SetZero();
        m_ppCovInv[i] = new TMatrix(nFeatureCnt, nFeatureCnt);
        m_ppCovInv[i]->SetZero();
    }
    m_dbPrior = new double[nClassCnt];
    m_dbCovDet = new double[nClassCnt];
    m_dbLnSigma = new double[nClassCnt];
    m_dbLnPrior = new double[nClassCnt];


    // hard coded initialization, I might get rid of this eventually but for now it doesn't matter because this should
    // always be overwritten.
    double glevl = 255. / (m_nClassCnt + 2);
    for (int i=0; i<m_nClassCnt;i++){
        for (int j=0; j<m_nFeatureCnt; j++){
            (*m_ppMu[i])(j) = (i+1)*glevl;
            for(int k=0; k<m_nFeatureCnt; k++){
                if (j == k)
                    (*m_ppCovInv[i])(j,k) = 0.005;
                else
                    (*m_ppCovInv[i])(j,k) = 0.0;

            }
        }
        m_dbCovDet[i] = sqrt(m_ppCovInv[i]->Inv()->Determinant());
        m_dbLnSigma[i] = log(m_dbCovDet[i]);
        m_dbPrior[i] = 0.5;
        m_dbLnPrior[i] = log(0.5);
    }
}

GMM::~GMM(){
    if (m_ppMu)
    {
        for ( int i = 0; i < m_nClassCnt; i++ ) {
            delete m_ppMu[i];
            delete m_ppIncSlope[i];
        }

        delete[] m_ppMu;
        delete[] m_ppIncSlope;
    }
    if (m_ppCovInv)
    {
        for ( int i = 0; i < m_nClassCnt; i++ )
            delete m_ppCovInv[i];
        delete[] m_ppCovInv;
    }
    delete[] m_dbPrior;
    delete[] m_dbCovDet;
    delete[] m_dbLnSigma;
    delete[] m_dbLnPrior;
}

void GMM::CopyFromGaussClassifier(TGaussianRegionClustering *pClassifier)
{

    TVector** ppMu = pClassifier->GetMean();
    TMatrix** ppCovInv = pClassifier->GetCovInv();
    double* dbPrior = pClassifier->GetPrior();
    for (int i=0; i<m_nClassCnt;i++){
        for (int j=0; j<m_nFeatureCnt; j++){
            (*m_ppMu[i])(j) = (*ppMu[i])(j);
            for(int k=0; k<m_nFeatureCnt; k++){
                (*m_ppCovInv[i])(j,k) = (*ppCovInv[i])(j,k);
            }
        }
        m_dbCovDet[i] = sqrt(m_ppCovInv[i]->Inv()->Determinant());
        m_dbLnSigma[i] = log(m_dbCovDet[i]);
        m_dbPrior[i] = dbPrior[i];
        m_dbLnPrior[i] = log(dbPrior[i]);
        m_ppIncSlope[i]->SetZero();
    }
}


/* PyGMM IMPLEMENTATION */

PyGMM::PyGMM(GMM* gmm)
{
    m_nClass = gmm->m_nClassCnt;
    m_nFeat = gmm->m_nFeatureCnt;

    py::array_t<float, py::array::c_style> Mu({m_nClass, m_nFeat});
    py::array_t<float, py::array::c_style> CovInv({m_nClass, m_nFeat, m_nFeat});
    py::array_t<float, py::array::c_style> Prior(m_nClass);
    py::array_t<float, py::array::c_style> CovDet(m_nClass);
    py::array_t<float, py::array::c_style> LnSigma(m_nClass);
    py::array_t<float, py::array::c_style> LnPrior(m_nClass);
    py::array_t<float, py::array::c_style> IncSlope({m_nClass, m_nFeat});

    auto rMu = Mu.mutable_unchecked<2>();
    for (int i=0;i<m_nClass;i++){
        for(int j=0; j< m_nFeat; j++){
            rMu(i,j) = (float)(*(gmm->m_ppMu[i]))(j);
        }
    }

    auto rCovInv = CovInv.mutable_unchecked<3>();
    for (int i=0;i<m_nClass;i++){
        for(int j=0; j< m_nFeat; j++){
            for(int k=0; k<m_nFeat; k++){
                rCovInv(i,j,k) = (float)(*(gmm->m_ppCovInv[i]))(j,k);
            }
        }
    }

    auto rPrior = Prior.mutable_unchecked<1>();
    for (int i=0;i<m_nClass;i++){
        rPrior(i) = (float) gmm->m_dbPrior[i];
    }

    auto rCovDet = CovDet.mutable_unchecked<1>();
    for (int i=0;i<m_nClass;i++){
        rCovDet(i) = (float) gmm->m_dbCovDet[i];
    }

    auto rLnSigma = LnSigma.mutable_unchecked<1>();
    for (int i=0;i<m_nClass;i++){
        rLnSigma(i) = (float) gmm->m_dbLnSigma[i];
    }

    auto rLnPrior = LnPrior.mutable_unchecked<1>();
    for (int i=0;i<m_nClass;i++){
        rLnPrior(i) = (float) gmm->m_dbLnPrior[i];
    }

    auto rInc = IncSlope.mutable_unchecked<2>();
    for (int i=0;i<m_nClass;i++){
        for(int j=0; j< m_nFeat; j++){
            rInc(i,j) = (float)(*(gmm->m_ppIncSlope[i]))(j);
        }
    }

    m_Mu = Mu;
    m_CovInv = CovInv;
    m_Prior = Prior;
    m_CovDet = CovDet;
    m_LnSigma = LnSigma;
    m_LnPrior = LnPrior;
    m_IncSlope = IncSlope;
}

GMM* PyGMM::togmm() {
    GMM* gmm = new GMM(m_nClass, m_nFeat);

    auto rMu = m_Mu.unchecked<2>();
    auto rCovInv = m_CovInv.unchecked<3>();
    auto rPrior = m_Prior.unchecked<1>();
    auto rCovDet = m_CovDet.unchecked<1>();
    auto rLnSigma = m_LnSigma.unchecked<1>();
    auto rLnPrior = m_LnPrior.unchecked<1>();
    auto rInc = m_IncSlope.unchecked<2>();

    for (int c=0; c<m_nClass;c++){

        TVector* vMuc = new TVector(m_nFeat, true);
        vMuc->SetZero();
        TVector* vInc = new TVector(m_nFeat, true);
        vInc->SetZero();

        TMatrix* vCic = new TMatrix(m_nFeat, m_nFeat);
        vCic->SetZero();

        for (int i=0; i<m_nFeat; i++){
            (*vMuc)(i) = (double) rMu(c,i);
            (*vInc)(i) = (double) rInc(c,i);
            for(int j=0; j<m_nFeat; j++){
                (*vCic)(i,j) = (double) rCovInv(c,i,j);
            }
        }

        gmm->m_ppMu[c] = vMuc;
        gmm->m_ppCovInv[c] = vCic;
        gmm->m_dbPrior[c] = (double) rPrior(c);
        gmm->m_dbCovDet[c] = (double) rCovDet(c);
        gmm->m_dbLnSigma[c] = (double) rLnSigma(c);
        gmm->m_dbLnPrior[c] = (double) rLnPrior(c);
        gmm->m_ppIncSlope[c] = vInc;
    }
    return gmm;
}



PyGMM::PyGMM(int nClass, int nFeat, py::array_t<float> Mu, py::array_t<float> CovInv, py::array_t<float> Prior,
             py::array_t<float> CovDet, py::array_t<float> LnSigma, py::array_t<float> LnPrior, py::array_t<float> IncSlope) {
    m_nClass = nClass;
    m_nFeat = nFeat;
    m_Mu = Mu;
    m_CovInv = CovInv;
    m_Prior = Prior;
    m_CovDet = CovDet;
    m_LnSigma = LnSigma;
    m_LnPrior = LnPrior;
    m_IncSlope = IncSlope;
}

PyGMM::PyGMM(py::dict GMMParams){

    // todo: add error checking to this function
    m_nClass = py::int_(GMMParams["nClass"]);
    m_nFeat = py::int_(GMMParams["nFeat"]);
    m_Mu = py::array_t<float>(GMMParams["Mu"]);
    m_CovInv = py::array_t<float>(GMMParams["CovInv"]);
    m_Prior = py::array_t<float>(GMMParams["Prior"]);
    m_CovDet = py::array_t<float>(GMMParams["CovDet"]);
    m_LnSigma = py::array_t<float>(GMMParams["LnSigma"]);
    m_LnPrior = py::array_t<float>(GMMParams["LnPrior"]);
    m_IncSlope = py::array_t<float>(GMMParams["IncSlope"]);
}

py::dict PyGMM::todict(){
    py::dict GMMParams;

    GMMParams["nClass"] = m_nClass;
    GMMParams["nFeat"] = m_nFeat;
    GMMParams["Mu"] = m_Mu;
    GMMParams["CovInv"] = m_CovInv;
    GMMParams["Prior"] = m_Prior;
    GMMParams["CovDet"] = m_CovDet;
    GMMParams["LnSigma"] = m_LnSigma;
    GMMParams["LnPrior"] = m_LnPrior;
    GMMParams["IncSlope"] = m_IncSlope;

    return GMMParams;
}