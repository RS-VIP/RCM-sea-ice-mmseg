//
// Created by max on 2021-08-13.
//

#include "TextureWrapper.h"


glcitr_param::glcitr_param()
{
    // default glcitr parameters
    m_nDist = 0;
    *m_nSpat;
    m_qGray = 64;
    m_pixDist = 1;
    m_Xwinsize = 25;
    m_Ywinsize = 25;

    m_basm = true;
    m_bcon = true;
    m_bcor = true;
    m_bdis = true;
    m_bent = true;
    m_bhom = true;
    m_binv = true;
    m_bmu = true;
    m_bstd = true;

    m_bnormflagGLC = true;
    m_bpadflag = true;
    m_bdirInv = true;  //direction invariant

    // m_bSaveImages = true;
    // default config has 4 radii, 6 angular orientations
    m_nGlcitrFeatCnt = 0; // initialize to 0 for some reason??
}

glcitr_param::glcitr_param(py::dict param_in) {
    // default glcitr parameters
    m_nDist = 0;
    *m_nSpat;
    m_qGray = py::int_(param_in["QGRAY"]); //64
    m_pixDist = py::int_(param_in["PIXDIST"]);//1
    m_Xwinsize = py::int_(param_in["XWINSZ"]);//25
    m_Ywinsize = py::int_(param_in["YWINSZ"]);//25

    m_basm = py::bool_(param_in["ASM"]);
    m_bcon = py::bool_(param_in["CON"]);
    m_bcor = py::bool_(param_in["COR"]);
    m_bdis = py::bool_(param_in["DIS"]);
    m_bhom = py::bool_(param_in["HOM"]);
    m_binv = py::bool_(param_in["INV"]);
    m_bent = py::bool_(param_in["ENT"]);
    m_bmu  = py::bool_(param_in["MU"]);
    m_bstd = py::bool_(param_in["STD"]);

    m_bnormflagGLC = py::bool_(param_in["NORM"]);
    m_bpadflag = true;
    m_bdirInv = true;  //direction invariant

    // m_bSaveImages = true;
    // default config has 4 radii, 6 angular orientations
    m_nGlcitrFeatCnt = 0; // initialize to 0 for some reason??
}


py::array_t<float> extract_glcitr(py::array_t<float> aSrcImg, py::dict glc_param_dict={})
{
    // univariate input image
    TGrayImage<float>* pSrcImg = PyArray2GrayImage(aSrcImg);

    glcitr_param* gp;
    if (glc_param_dict.empty()) gp = new glcitr_param(); // default parameters
    else gp = new glcitr_param(glc_param_dict);

    TGrayImage<FLOAT>** ppFeatImgGlcitr = NULL;
    TTextureFeatExtract* pFeatExtractor = new TTextureFeatExtract();

    int usemsr[TTextureFeatExtract::NUM_TEX_MEASURES] = {0};
    int nMsr = 0;
    if (gp->m_basm)
    {
        usemsr[nMsr] = TTextureFeatExtract::ASM;
        nMsr++;
    }
    if (gp->m_bcon)
    {
        usemsr[nMsr] = TTextureFeatExtract::CON;
        nMsr++;
    }
    if (gp->m_bcor)
    {
        usemsr[nMsr] = TTextureFeatExtract::COR;
        nMsr++;
    }
    if (gp->m_bdis)
    {
        usemsr[nMsr] = TTextureFeatExtract::DIS;
        nMsr++;
    }
    if (gp->m_bent)
    {
        usemsr[nMsr] = TTextureFeatExtract::ENT;
        nMsr++;
    }
    if (gp->m_bhom)
    {
        usemsr[nMsr] = TTextureFeatExtract::HOM;
        nMsr++;
    }
    if (gp->m_binv)
    {
        usemsr[nMsr] = TTextureFeatExtract::INV;
        nMsr++;
    }
    if (gp->m_bmu)
    {
        usemsr[nMsr] = TTextureFeatExtract::MU;
        nMsr++;
    }
    if (gp->m_bstd)
    {
        usemsr[nMsr] = TTextureFeatExtract::STD;
        nMsr++;
    }
    int temp_nSpat[8] = {gp->m_pixDist, 0, gp->m_pixDist, gp->m_pixDist, 0, gp->m_pixDist, -gp->m_pixDist, -gp->m_pixDist};
    gp->m_nSpat = temp_nSpat;

    if (gp->m_bdirInv){
        // isotropic features - all directions averaged
        gp->m_nDist = 1;
    }else{
        // directional features
        gp->m_nDist = 4;
    }

    ppFeatImgGlcitr = pFeatExtractor->Glcitr(gp->m_bpadflag, pSrcImg, gp->m_qGray,  usemsr, nMsr, gp->m_Xwinsize, gp->m_Ywinsize, gp->m_nSpat, gp->m_nDist, gp->m_bnormflagGLC, gp->m_bdirInv, "");

    gp->m_nGlcitrFeatCnt = nMsr*(gp->m_nDist);  //nDist;

    py::array_t<float> aFeatResult = GrayImageArray2PyArray(ppFeatImgGlcitr, gp->m_nGlcitrFeatCnt);


    delete pSrcImg;
    for(int i=0; i<gp->m_nGlcitrFeatCnt; i++) delete ppFeatImgGlcitr[i];
    delete ppFeatImgGlcitr;
    delete gp;
    delete pFeatExtractor;
    return aFeatResult;
}



/// convert a univariate TGrayImage to a py::array.
// template <typename T> py::array_t<T> GrayImage2PyArray(TGrayImage<T>* pGrayImg){
//     int H = pGrayImg->Height(), W = pGrayImg->Width();
//     py::array_t<T, py::array::c_style> pyArray({H,W});

//     auto ra = pyArray.template mutable_unchecked<2>();

//     for (int i = 0; i < H; i++) {
//         for (int j = 0; j < W; j++) {
//             ra(i,j) = pGrayImg->GetPixelData(j,i);
//         }
//     }
//     return pyArray;
// }

py::array_t<float> CapsuleTest(){

    float* buf = new float[1000*1000*10];

    for(int i=0; i<1000*1000*10; i++)
        buf[i] = 0.0;
    
    // Create a Python object that will free the allocated
    // memory when destroyed:
    py::capsule free_when_done(buf, [](void *f) {
        float *buf = reinterpret_cast<float*>(f);
        std::cerr << "freeing memory @ " << f << "\n";
        delete[] buf;
    });

    py::array_t<float, py::array::c_style> arr(     {1000,1000,10},
                                                    buf,
                                                    free_when_done );
    return arr;                                  
}

py::array_t<float> NoCapsuleTest(){

    float* buf = new float[1000*1000*10];

    for(int i=0; i<1000*1000*10; i++)
        buf[i] = 0.0;

    py::array_t<float, py::array::c_style> arr(     {1000,1000,10},
                                                    buf);
    return arr;                                  
}


py::array_t<float> test_array_capsule(){
    // py::array_t<float> arr = CapsuleTest();
    py::array_t<float> arr = NoCapsuleTest();
    return arr;
}


gabor_param::gabor_param()
{
    // default parameters
    m_dbPostLamda = 2.0/3.0;
    m_dblow_rad = 6.0;
    m_dbhigh_rad = 9.0;
    m_dbinc_rad = 1.0;
    m_dblow_ang = 0.0;
    m_dbhigh_ang = 179.0;
    m_dbinc_ang = 30.0;

    m_bNormFlagGabor = false;
    m_bUnitLamdaFlag = false;
    m_bLinSpaceFlag = false;
    m_bPostFltFlag = false;
    m_bRealOnlyFlag = false;
    m_bDCRemFlag = true;

    m_bSaveImages = false;
    // default config has 4 radii, 6 angular orientations
    m_nGaborFeatCnt = 0; // initialize to 0, this value is passed by reference to the feature extraction function and
    // incremented as each new feature is computed
}

py::array_t<float> extract_gabor(py::array_t<float> aSrcImg)
{
    // NOTICE: there is something weird going on with the gabor features extracted with this function, the results seem to
    // depend on the size of the input image... not sure if I am passing in wrong parameters or what, I don't have a good reference
    // for how the GaborBovik function is supposed to be used. Anyway I dont really need gabor features so I'll leave it for now

    // univariate input image
    TGrayImage<float>* pSrcImg = PyArray2GrayImage(aSrcImg);

    TTextureFeatExtract* pFeatExtractor = new TTextureFeatExtract();

    gabor_param* gp = new gabor_param();
    const char* bandname = "/home/max/code/cpp/magic_port/magic_py/";

    TGrayImage<float>** ppFeatImgGabor = nullptr;

    ppFeatImgGabor = pFeatExtractor->GaborBovik(pSrcImg, gp->m_bSaveImages, gp->m_bRealOnlyFlag, gp->m_bPostFltFlag, gp->m_dbPostLamda, gp->m_dblow_rad, gp->m_dbhigh_rad,
                                                gp->m_dbinc_rad, gp->m_bLinSpaceFlag, gp->m_dblow_ang, gp->m_dbhigh_ang, gp->m_dbinc_ang, gp->m_bUnitLamdaFlag, gp->m_bDCRemFlag,
                                                gp->m_nGaborFeatCnt, gp->m_bNormFlagGabor, bandname);

    py::array_t<float> aFeatResult = GrayImageArray2PyArray(ppFeatImgGabor, gp->m_nGaborFeatCnt);

    delete pSrcImg;
    delete gp;
    delete ppFeatImgGabor;
    delete pFeatExtractor;

    return aFeatResult;
}
