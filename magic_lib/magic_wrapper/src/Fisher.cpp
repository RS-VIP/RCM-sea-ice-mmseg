//
// Created by max on 2021-08-16.
//

#include "Fisher.h"

double ComputeMinFisher(GMM* gmm)
{
    int nFeatureCnt = gmm->m_nFeatureCnt;
    int nClassCnt = gmm->m_nClassCnt;
    double dbMin = 100000000;  // a big enough number
    double dbVal;
    int i, j;
    if (nFeatureCnt == 1){
        for (i = 0; i < nClassCnt - 1; i ++){
            for (j = i + 1; j < nClassCnt; j ++){
                dbVal = (*gmm->m_ppMu[i])(0) - (*gmm->m_ppMu[j])(0);
                dbVal *= dbVal;
                dbVal /= (1 / (*gmm->m_ppCovInv[i])(0, 0) * (gmm->m_dbPrior[i]/(gmm->m_dbPrior[i]+gmm->m_dbPrior[j])) +
                        1 / (*gmm->m_ppCovInv[j])(0, 0) * (gmm->m_dbPrior[j]/(gmm->m_dbPrior[i]+gmm->m_dbPrior[j])));

                if (dbVal < dbMin)
                    dbMin = dbVal;
            }
        }
    }else if (nFeatureCnt > 1) //QIN_ADD
        {
        TMatrix* WSMat1 = new TMatrix(nFeatureCnt, nFeatureCnt);
        TMatrix* WSMat2 = new TMatrix(nFeatureCnt, nFeatureCnt);
        TVector* BSVec = new TVector(nFeatureCnt, false);
        TVector* BSVec_T = new TVector(nFeatureCnt, false);
        TMatrix *INVWSMat;
        for ( i = 0; i < nClassCnt - 1; i++ )
        {
            for ( j = i + 1; j < nClassCnt; j++ )
            {
                for ( int k = 0; k < nFeatureCnt; k++ )
                    (*BSVec)(k) = ((*gmm->m_ppMu[i])(k) - (*gmm->m_ppMu[j])(k));

                gmm->m_ppCovInv[i]->Inv(WSMat1);
                *WSMat1 *= (gmm->m_dbPrior[i]/(gmm->m_dbPrior[i]+gmm->m_dbPrior[j]));
                gmm->m_ppCovInv[j]->Inv(WSMat2);
                *WSMat2 *= (gmm->m_dbPrior[j]/(gmm->m_dbPrior[i]+gmm->m_dbPrior[j]));
                (*WSMat1) += (*WSMat2);

                INVWSMat = WSMat1->Inv();

                BSVec->Mul(BSVec_T, INVWSMat);
                dbVal = BSVec_T->VecMul(BSVec);

                if (dbVal < 0)
                {
                    // printf("Fisher calculation wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    exit(1);
                }

                if (dbVal < dbMin)
                    dbMin = dbVal;
            }
        }
        delete WSMat1;
        delete WSMat2;
        delete INVWSMat;
        delete BSVec;
        delete BSVec_T;
        }
    else
    {
        //cout << "Error: The number of features should be above 0!" << endl;
        exit(1);
    }

    return dbMin;
}

double ComputeFisher(PyGMM* pygmm){
    GMM* gmm = pygmm->togmm();
    double fish = ComputeMinFisher(gmm);
    return fish;
}