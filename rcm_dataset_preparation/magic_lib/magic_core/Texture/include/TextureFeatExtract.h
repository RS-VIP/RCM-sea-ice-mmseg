/** \file
 *
 * \brief TTextureFeatExtract class.
 *
 * \todo TODO: Document this file.
 *
 */

#ifndef TTextureFeatExtract_Class
#define TTextureFeatExtract_Class

#include "GrayImage.h"
#include "matrix_gabor.h"
#include "SignalProcessing.h"

#define MAXFILENAMELEN 255

/// Class of texture feature extraction
/**
 * \todo TODO: Document this class.
 */
class TTextureFeatExtract
{
public:
	TTextureFeatExtract();

//Attribute
public:
	/// Measure identifyers
	enum glctxtid
	{
	  ASM,
	  ENT,
	  CON,
	  DIS,
	  HOM,
	  INV,
	  MU,
	  STD,
	  COR,
      NUM_TEX_MEASURES
	};

	/// Intermediary measures that must be carried to calculate real measures
	enum tmpglctxtid
	{
	  STDTEMP,
	  CORTEMP,
	  NUM_TEMP_MEASURES
	};

//Operations
public:
	//On univariate images
	//Ross Jan 2016
	//First function is overloaded to allow progressbar to be updated 
	// Max 2021: Gross, get this GUI nonsense out of here!!
	//TGrayImage<FLOAT>** Glcitr(System::Windows::Forms::ProgressBar^ progressbar, bool b_saveImages, BOOL bPadflag, TGrayImage<FLOAT>* pSrcImg, int qGray, int* usemsr, int nMsr, int nWinXSize, int nWinYSize, int* nSpat, int nDist, int nNormFlag, bool bdirInv=true, const char* source_band_name = "");
	TGrayImage<FLOAT>** Glcitr(BOOL bPadflag, 
		TGrayImage<FLOAT>* pSrcImg,
		int qGray, int* usemsr,
		int nMsr,
		int nWinXSize,
		int nWinYSize,
		int* nSpat,
		int nDist,
		int nNormFlag,
		bool bdirInv=true,
		const char* source_band_name = "");//The number of output feature images equals to nMsr*nDist 
	
	TGrayImage<FLOAT>** GaborFltrBank(BOOL bRealOnlyFlag,
		TGrayImage<FLOAT>* pSrcImg, 
		int side,
		double Ul, 
		double Uh,
		int nScaleNum, 
		int nOrientNum,
		BOOL bDCRemFlag,
		int nNormFlag); //The number of output feature images equals to nScaleNum*nOrientNum
	
	TGrayImage<FLOAT>** GaborBovik(TGrayImage<FLOAT>* pSrcImg,
		bool saveImages,
		BOOL bRealOnlyFlag,
		BOOL bPostFltFlag,
		double dbPostLamda,
		double low_rad, 
		double high_rad, 
		double inc_rad, 
		BOOL bLinSpaceFlag,
		double low_ang, 
		double high_ang,
		double inc_ang,
		BOOL bUnitLamdaFlag, 
		BOOL bDCRemFlag,
		int &nGaborCnt,
		int nNormFlag,
		const char* source_band_name = "");

	//On multivariate images
private:
	void TextureMeasures(float** ftrsout, int** img, int nXsize, int nYSize, int gnum, 
                         int nFlszX, int nFlszY, int nSpatX, int nSpatY, 
                         int* usemsr, int nMeas, int bDirInv);

	void moveout(float** pfOutBuf, int* iOutIndex, float* measures, int nMeas, int* msrbool);

	void build_dirinv(float* measures, int** glcm, float* tstats,
					  int** ppbyData, int gnum, int nFlszX, int nFlszY, int nSpat,
					  float*** lk, int* msrbool, int bDirInv);

	void build_single(float* measures, int** glcm, float* tstats,
                      int** ppbyData, int gnum, int nFlszX, int nFlszY, int nSpatX, int nSpatY,
					  float*** lk, int* msrbool, int bDirInv);

	void update_dirinv(float* measures, int** glcm, float* tstats, 
					   int** ppbyData, int nFlszX, int nFlszY, int nSpat,
					   int nPixOffset,   float*** lk, int* msrbool);

	void update_single(float* measures, int** glcm, float* tstats, 
                       int** ppbyData, int nFlszX, int nFlszY, int nSpatX, int nSpatY,
                       int nPixOffset,   float*** lk, int* msrbool);


	void findrange(int* startX, int* stopX, int* startY, int* stopY, int nFlszX, int nFlszY, int nSpatX, int nSpatY);

	void GaborFilteredImg(Matrix *FilteredImg_real, Matrix *FilteredImg_imag, Matrix *img, int side, double Ul, double Uh, int scale, int orientation, int flag);
	
	void Gabor(Matrix *Gr, Matrix *Gi, int s, int n, double Ul, double Uh, int scale, int orientation, int flag);

	double GaborFreq(double u, double v, double sigma2, double freq, double lambda2, double coeff);

	double GaborPost(double u, double v, double sigma2, double lambda2, double coeff, double postlamda);
public:
	void SetSaveFlag(bool bSaveGabor, bool bSaveGlcitr, char* sLocation);// added by sharif to set global flag for saving
	

private:
	bool m_bSaveGabor;// added by sharif
	bool m_bSaveGlictr;// added by sharif
	const char* m_SaveLocation; // Folder location to save features


public:
	virtual ~TTextureFeatExtract();
};

#endif
