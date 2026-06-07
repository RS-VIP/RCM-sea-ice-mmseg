//
// Created by max on 2021-08-13.
//

#ifndef MAGIC_PORT_TEXTUREEXTRACT_H
#define MAGIC_PORT_TEXTUREEXTRACT_H

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "TextureFeatExtract.h"
#include "ImageConversion.h"

//py::array_t<float> extract_glcitr(py::array_t<float> aSrcImg);
py::array_t<float> extract_glcitr(py::array_t<float> aSrcImg, py::dict glc_param);
py::array_t<float> extract_gabor(py::array_t<float> aSrcImg);

py::array_t<float> CapsuleTest();
py::array_t<float> NoCapsuleTest();
py::array_t<float> test_array_capsule();

class gabor_param {
public:
    gabor_param();

    // todo: python conversion functions implementation
    gabor_param(py::dict param_in);
    //    py::dict gabor_param_dict();

    // variables for Gabor texture tab
    double m_dbPostLamda;
    double m_dblow_rad;
    double m_dbhigh_rad;
    double m_dbinc_rad;
    double m_dblow_ang;
    double m_dbhigh_ang;
    double m_dbinc_ang;

    bool m_bNormFlagGabor;
    bool m_bUnitLamdaFlag;
    bool m_bLinSpaceFlag;
    bool m_bPostFltFlag;
    bool m_bRealOnlyFlag;
    bool m_bDCRemFlag;

    bool m_bSaveImages;
    int m_nGaborFeatCnt;
};

class glcitr_param {
public:
    glcitr_param();
    glcitr_param(py::dict param_in);
    py::dict glcitr_param_dict();

    // variable for GLCITR Tab
    int m_nDist;
    int *m_nSpat;
    int m_qGray;
    int m_pixDist;
    int m_Xwinsize;
    int m_Ywinsize;

    bool m_basm;
    bool m_bcon;
    bool m_bcor;
    bool m_bdis;
    bool m_bent;
    bool m_bhom;
    bool m_binv;
    bool m_bmu;
    bool m_bstd;

    bool m_bnormflagGLC;
    bool m_bpadflag;
    bool m_bdirInv;  //direction invariant

    int m_nGlcitrFeatCnt;
};



#endif //MAGIC_PORT_TEXTUREEXTRACT_H
