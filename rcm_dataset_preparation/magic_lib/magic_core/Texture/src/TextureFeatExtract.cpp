#pragma warning( disable : 4244 )

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<vector>
#include "TextureFeatExtract.h"
#include <fstream> 
#include <ctime>

//Threshold for detecting a zero value 10^(-12)
#define ZERO_THRESH 1e-12

//using namespace System::Windows::Forms;

// Constructor
/////////////////////////////////////////////////////////////////////////////
TTextureFeatExtract::TTextureFeatExtract()
{
	m_bSaveGabor = false;// added by sharif
	m_bSaveGlictr = false;// added
	m_SaveLocation = "";

}

// Destructor
/////////////////////////////////////////////////////////////////////////////
TTextureFeatExtract :: ~TTextureFeatExtract()
{
}


void TTextureFeatExtract::SetSaveFlag(bool bSaveGabor, bool bSaveGlcitr, char* sLocation)// added by sharif to set global flag for saving feature
{
	m_bSaveGabor = bSaveGabor;
	m_bSaveGlictr = bSaveGlcitr;
	m_SaveLocation = sLocation;
}


/***************** Operations on univariate images ******************/
/********************************************************************/

// Calculates Haralick features using the faster iterative method
// Param in:

// Return:
// Note: the number of output feature images equals to nMsr*nDist 
/////////////////////////////////////////////////////////////////////////////
//Ross Jan 2016 Overloaded this function
// Max 2021 Got rid of the overloaded version, sorry Ross
/*
TGrayImage<FLOAT>** TTextureFeatExtract :: Glcitr(System::Windows::Forms::ProgressBar^ progressbar, bool b_saveImages, BOOL bPadflag, TGrayImage<FLOAT>* pSrcImg, int qGray, int* usemsr, int nMsr, int nWinXSize, int nWinYSize, int* nSpat, int nDist, int nNormFlag, bool bdirInv, const char* source_band_name)
{
	BOOL bFeatImgOut = 0;
	char statistics[9][4];
	float statLimits[9][2];

	int i, j, k, m;

	TGrayImage<int>* pImg = new TGrayImage<int>(pSrcImg);

	std::vector<TGrayImage<FLOAT>*> images;
	std::vector<string> bandnames;
	char featimgfile_bil[MAXFILENAMELEN];
	if (m_SaveLocation != NULL && b_saveImages)
	{
		strcpy_s(featimgfile_bil, MAXFILENAMELEN, m_SaveLocation);
		strcat_s(featimgfile_bil, MAXFILENAMELEN, "Glc_");
	}

	if (bPadflag){ 
		int nHalfYSize = (int) (nWinYSize - 1)/2;
		int nHalfXSize = (int) (nWinXSize - 1)/2;

		pImg->Pad( nHalfYSize, nHalfYSize, nHalfXSize, nHalfXSize, PADMIRROR);

		nWinYSize = 2*nHalfYSize + 1;
		nWinXSize = 2*nHalfXSize + 1;
	}

	int nImgHeight = pImg->Height();
	int nImgWidth = pImg->Width();
	int feat_width = nImgWidth - nWinXSize + 1;
	int feat_height = nImgHeight - nWinYSize + 1;
	int numSamps = feat_width*feat_height;

	TGrayImage<int>* bmpImg = new TGrayImage<int>(feat_width, feat_height);


	//initialize resulting feature images
	TGrayImage<FLOAT>** pRet = new TGrayImage<FLOAT>*[nMsr*nDist];
	for (i = 0; i < nMsr*nDist; i++)
		pRet[i] = new TGrayImage<FLOAT>(nImgWidth - nWinXSize + 1, nImgHeight - nWinYSize + 1);

	float** ftrs = new float*[numSamps];

	//Progressbar feature implemented by Ross in Jan 2016
	progressbar->Enabled = true;
	progressbar->PerformStep();

	for ( i = 0; i < numSamps; i++ )
	{
		ftrs[i] = new float[nMsr];
	}

	progressbar->Value = 5;
	progressbar->PerformStep();

	int** img = pImg->GetData();

	int minGrayVal = 100000000, maxGrayVal = -100000000;

	for ( i = 0; i < nImgHeight; i++ ){
		for (j = 0; j < nImgWidth; j++ ){
			if ( img[i][j] < minGrayVal )
				minGrayVal = img[i][j];
			else if ( img[i][j] > maxGrayVal )
				maxGrayVal = img[i][j];
		}
	}
	progressbar->Value = 10;
	progressbar->PerformStep();

	//The following two lines about maxGrayVal and minGrayVal are added by Zhijie for the SAR SeaIce Imagery due to the images are not consistently histogram stretched.
	//These two lines means that we are not gonna stretch the input image dynamically with the max and min value in the image, because this may change the information of the texture signal. 
	//However this is wierd and hardcoded. Attention has to be paid to these two lines, because first we may wanna stretch the image dynamically and second the image values may be within the range from 0 to 1 rather than from 0 to 255. 
	maxGrayVal = 255;
	minGrayVal = 0;

	//setup max and min theoretical bounds on GLCM feature measures
	for ( i = 0; i < nMsr; i++ ){
		switch ( usemsr[i] ){
		case ASM: strcpy_s(statistics[i], 4, "ASM");
			statLimits[ASM][0] = (float) 1.0 / ( (float) (qGray*qGray) );
			statLimits[ASM][1] = 1.0;
			break;
		case ENT: strcpy_s(statistics[i], 4, "ENT");
			statLimits[ENT][0] = 0.0;
			statLimits[ENT][1] = 2*log( (float) qGray);
			break;
		case CON: strcpy_s(statistics[i], 4, "CON");
			statLimits[CON][0] = 0.0;
			statLimits[CON][1] =  (float) ((qGray - 1)*(qGray - 1));
			break;
		case DIS: strcpy_s(statistics[i], 4, "DIS");
			statLimits[DIS][0] = 0.0;
			statLimits[DIS][1] = fabs ( (float) (qGray - 1) );
			break;
		case HOM: strcpy_s(statistics[i], 4, "HOM");
			statLimits[HOM][0] = (float) ( 1.0 / ( 1.0 + (float) ( (qGray - 1)*(qGray - 1) ) ) );
			statLimits[HOM][1] = 1.0;
			break;
		case INV: strcpy_s(statistics[i], 4, "INV");
			statLimits[INV][0] = (float) ( 1.0 / ( 1.0 + fabs ( (float) (qGray - 1) ) ) );
			statLimits[INV][1] = 1.0;
			break;
		case MU:  strcpy_s(statistics[i], 4, "MU");
			statLimits[MU][0] = 0.0;
			statLimits[MU][1] = (float) qGray;
			break;
		case STD: strcpy_s(statistics[i], 4, "STD");
			statLimits[STD][0] = 0.0;
			statLimits[STD][1] = (float) ( ( (float) qGray + 1.0 )*( (float) qGray - 1.0 ) / 12.0 );
			break;
		case COR: strcpy_s(statistics[i], 4, "COR");
			statLimits[COR][0] = -1.0;
			statLimits[COR][1] = 1.0;
			break;
		default: {
			fprintf(stderr,"*** ERROR: The specified statistic is not in the list!\n");
			exit(1);
				 }
		}
	}
	progressbar->Value = 15;
	progressbar->PerformStep();

	if ( ( qGray  < 1 ) || ( qGray > ( maxGrayVal - minGrayVal + 1 ) ) )
	{
		MessageBox::Show("Invalid grey level quantization.  It should be an integer between 1 and the source bitdepth (usually 256).  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf(stderr,"*** ERROR: (%d) is wrong quantized grey level! It should be an interger between 1 and %d!\n",qGray, maxGrayVal - minGrayVal + 1);
		exit (1); 
	}

	for ( i = 0; i < nDist; i++ ){
		if ( nSpat[i] > nWinXSize - 1 ){
			MessageBox::Show("Specified inter-pixel distance exceeds X direction windows size.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
			fprintf(stderr,"Error: Specified inter-pixel distance: %d exceeds X direction windows size: %d !\n", nSpat[i], nWinXSize); 
			exit(1);
		}
		if ( nSpat[i] > nWinYSize - 1 ){ 
			MessageBox::Show("Specified inter-pixel distance exceeds Y direction windows size.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
			fprintf(stderr,"Error: Specified inter-pixel distance: %d exceeds Y direction windows size: %d !\n", nSpat[i], nWinYSize); 
			exit(1);
		}
	}


	//bin the source image to the GLCM quantization levels
	float dbQBin = (float) ( qGray - 1 ) / ( maxGrayVal - minGrayVal );
	for ( i = 0; i < nImgHeight; i++ )
		for (j = 0; j < nImgWidth; j++ )
			img[i][j]= (int) ( dbQBin * ( img[i][j] - minGrayVal ) );

	progressbar->Value = 20;
	progressbar->Step = 30 / nDist / nMsr + 1;

	for ( m = 0; m < nDist; m++ ){
		//extract GLCM features
		TextureMeasures(ftrs, img, nImgWidth, nImgHeight, qGray, nWinXSize, nWinYSize, nSpat[m*2], nSpat[m*2+1], usemsr, nMsr, bdirInv);

		for ( k = 0; k < nMsr; k++ ){
			progressbar->PerformStep();

			//copy ftrs to pRet.  Normalize if signified
			for ( i = 0; i < feat_height; i++ ){
				for ( j = 0; j < feat_width; j++ ){
					if (nNormFlag){
						//normalize to [0,1] based on theoretical statistical bounds
						if ( fabs( statLimits[usemsr[k]][1] - statLimits[usemsr[k]][0] ) > ZERO_THRESH ){
							(*pRet[m*nMsr + k])(i, j) = (ftrs[i*feat_width + j][k] - statLimits[usemsr[k]][0]) / (statLimits[usemsr[k]][1] - statLimits[usemsr[k]][0]);
						}else{
							(*pRet[m*nMsr + k])(i, j) = 0;
						}
					}else{
						(*pRet[m*nMsr + k])(i, j) = ftrs[i*feat_width + j][k];
					}
				}
			}

			if (!this->m_bSaveGlictr)
				continue;


			char band_descriptor[MAXFILENAMELEN];
			char full_file_path[MAXFILENAMELEN];
			if (b_saveImages)
			{

				if (bdirInv)
					sprintf(band_descriptor, "GLCM_meas=%s_step=%i_winX=%d_winY=%d_qnt=%d_srcimg=%s", statistics[k], nSpat[m * 2], nWinXSize, nWinYSize, qGray, source_band_name);
				else
					sprintf(band_descriptor, "GLCM_meas=%s_stepX=%i_stepY=%i_winX=%d_winY=%d_qnt=%d_srcimg=%s", statistics[k], nSpat[m * 2], nSpat[m * 2 + 1], nWinXSize, nWinYSize, qGray, source_band_name);

				if (CanShowStatus())//send filename to GUI
					DisplayStatus(band_descriptor,-1);

				//push result to .bil image stack
				images.push_back(pRet[m*nMsr + k]);					
				string bandname(band_descriptor);
				bandnames.push_back(bandname);
			
				//save to csv file

				sprintf(full_file_path, "%s%s.csv", m_SaveLocation, band_descriptor);
				pRet[m*nMsr + k]->WritePixelValuesToFile(full_file_path);
			}

			//normalize to [0,255] based on image max,min.  Then save as bmp.
			float normmin= 99999999999;
			float normmax=-99999999999;
			for (i = 0; i < feat_height*feat_width; i ++ ){
				if (normmin > ftrs[i][k])
					normmin = ftrs[i][k];
				if (normmax < ftrs[i][k])
					normmax = ftrs[i][k];
			}

			//If flagged to save images to disk space, do so
			if (b_saveImages)
			{
				for (i = 0; i < feat_height; i++) {
					for (j = 0; j < feat_width; j++) {
						if (normmax - normmin > ZERO_THRESH) {
							float zero_to_one = (ftrs[i*feat_width + j][k] - normmin) / (normmax - normmin);
							(*bmpImg)(i, j) = (int)(255 * zero_to_one + 0.5);
						}
						else {
							(*bmpImg)(i, j) = 0;
						}
					}
				}
				sprintf(full_file_path, "%s%s.bmp", m_SaveLocation, band_descriptor);
				bmpImg->Save(full_file_path);
			}

		}
	}

	if (this->m_bSaveGlictr && b_saveImages){
		//save all feature bands to multiband bil file as 32 bit floats
		MultivariateReadWrite<FLOAT> mrw;
		char featimgfile_bil[MAXFILENAMELEN];
		time_t now = time(0);  // get time now
		struct tm  tstruct;
		tstruct = *localtime(&now);
		char time_buf[200];
		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d_%H.%M.%S", &tstruct);
		sprintf(featimgfile_bil, "%sGLCMBank_%s.bil", m_SaveLocation, time_buf);
		mrw.Write(string(featimgfile_bil), images, ImageHeaderInfo::IMG_FLOAT32, bandnames);

	}
	delete bmpImg;

	progressbar->Step = 1;
	int valueLeft = 100 - progressbar->Value;
	int whenStep = numSamps / valueLeft;
	int steps = 0;

	// Free memory;
	delete pImg;
	for ( i = 0; i < numSamps; i++ )  //this takes super long
	{
		if (ftrs[i]) 
			delete[] ftrs[i];
		if (i > whenStep * steps)
		{
			progressbar->PerformStep();
			steps++;
		}
	}
	if (ftrs) delete[] ftrs;

	progressbar->Value = 0;

	return pRet;
}
*/

TGrayImage<FLOAT>** TTextureFeatExtract :: Glcitr(BOOL bPadflag, TGrayImage<FLOAT>* pSrcImg, int qGray, int* usemsr, int nMsr, int nWinXSize, int nWinYSize, int* nSpat, int nDist, int nNormFlag, bool bdirInv, const char* source_band_name)
{
	BOOL bFeatImgOut = 0;

	char statistics[9][4];
	float statLimits[9][2];

	int i, j, k, m;

	TGrayImage<int>* pImg = new TGrayImage<int>(pSrcImg);


	std::vector<TGrayImage<FLOAT>*> images;
	std::vector<string> bandnames;
	char featimgfile_bil[MAXFILENAMELEN];

	strcpy( featimgfile_bil, m_SaveLocation);
	strcat( featimgfile_bil, "Glc_");

	if (bPadflag){ 
		int nHalfYSize = (int) (nWinYSize - 1)/2;
		int nHalfXSize = (int) (nWinXSize - 1)/2;

		pImg->Pad( nHalfYSize, nHalfYSize, nHalfXSize, nHalfXSize, PADMIRROR);

		nWinYSize = 2*nHalfYSize + 1;
		nWinXSize = 2*nHalfXSize + 1;
	}

	int nImgHeight = pImg->Height();
	int nImgWidth = pImg->Width();
	int feat_width = nImgWidth - nWinXSize + 1;
	int feat_height = nImgHeight - nWinYSize + 1;
	int numSamps = feat_width*feat_height;

	TGrayImage<int>* bmpImg = new TGrayImage<int>(feat_width, feat_height);


	//initialize resulting feature images
	TGrayImage<FLOAT>** pRet = new TGrayImage<FLOAT>*[nMsr*nDist];  
	for ( i = 0; i < nMsr*nDist; i++ )
		pRet[i] = new TGrayImage<FLOAT>( nImgWidth - nWinXSize + 1, nImgHeight - nWinYSize + 1 );

	float** ftrs = new float*[numSamps];

	for ( i = 0; i < numSamps; i++ )
	{
		ftrs[i] = new float[nMsr];
	}

	int** img = pImg->GetData();

	int minGrayVal = 100000000, maxGrayVal = -100000000;

	for ( i = 0; i < nImgHeight; i++ ){
		for (j = 0; j < nImgWidth; j++ ){
			if ( img[i][j] < minGrayVal )
				minGrayVal = img[i][j];
			else if ( img[i][j] > maxGrayVal )
				maxGrayVal = img[i][j];
		}
	}

	//The following two lines about maxGrayVal and minGrayVal are added by Zhijie for the SAR SeaIce Imagery due to the images are not consistently histogram stretched.
	//These two lines means that we are not gonna stretch the input image dynamically with the max and min value in the image, because this may change the information of the texture signal. 
	//However this is wierd and hardcoded. Attention has to be paid to these two lines, because first we may wanna stretch the image dynamically and second the image values may be within the range from 0 to 1 rather than from 0 to 255. 
	maxGrayVal = 255;
	minGrayVal = 0;

	//setup max and min theoretical bounds on GLCM feature measures
	for ( i = 0; i < nMsr; i++ ){
		switch ( usemsr[i] ){
		case ASM: strcpy(statistics[i], "ASM");
			statLimits[ASM][0] = (float) 1.0 / ( (float) (qGray*qGray) );
			statLimits[ASM][1] = 1.0;
			break;
		case ENT: strcpy(statistics[i], "ENT");
			statLimits[ENT][0] = 0.0;
			statLimits[ENT][1] = 2*log( (float) qGray);
			break;
		case CON: strcpy(statistics[i], "CON");
			statLimits[CON][0] = 0.0;
			statLimits[CON][1] =  (float) ((qGray - 1)*(qGray - 1));
			break;
		case DIS: strcpy(statistics[i], "DIS");
			statLimits[DIS][0] = 0.0;
			statLimits[DIS][1] = fabs ( (float) (qGray - 1) );
			break;
		case HOM: strcpy(statistics[i], "HOM");
			statLimits[HOM][0] = (float) ( 1.0 / ( 1.0 + (float) ( (qGray - 1)*(qGray - 1) ) ) );
			statLimits[HOM][1] = 1.0;
			break;
		case INV: strcpy(statistics[i], "INV");
			statLimits[INV][0] = (float) ( 1.0 / ( 1.0 + fabs ( (float) (qGray - 1) ) ) );
			statLimits[INV][1] = 1.0;
			break;
		case MU:  strcpy(statistics[i], "MU");
			statLimits[MU][0] = 0.0;
			statLimits[MU][1] = (float) qGray;
			break;
		case STD: strcpy(statistics[i], "STD");
			statLimits[STD][0] = 0.0;
			statLimits[STD][1] = (float) ( ( (float) qGray + 1.0 )*( (float) qGray - 1.0 ) / 12.0 );
			break;
		case COR: strcpy(statistics[i], "COR");
			statLimits[COR][0] = -1.0;
			statLimits[COR][1] = 1.0;
			break;
		default: {
			fprintf(stderr,"*** ERROR: The specified statistic is not in the list!\n");
			exit(1);
				 }
		}
	}

	if ( ( qGray  < 1 ) || ( qGray > ( maxGrayVal - minGrayVal + 1 ) ) )
	{
		//MessageBox::Show("Invalid grey level quantization.  It should be an integer between 1 and the source bitdepth (usually 256).  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf(stderr,"*** ERROR: (%d) is wrong quantized grey level! It should be an interger between 1 and %d!\n",qGray, maxGrayVal - minGrayVal + 1);
		exit (1); 
	}

	for ( i = 0; i < nDist; i++ ){
		if ( nSpat[i] > nWinXSize - 1 ){
			//MessageBox::Show("Specified inter-pixel distance exceeds X direction windows size.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
			fprintf(stderr,"Error: Specified inter-pixel distance: %d exceeds X direction windows size: %d !\n", nSpat[i], nWinXSize); 
			exit(1);
		}
		if ( nSpat[i] > nWinYSize - 1 ){ 
			//MessageBox::Show("Specified inter-pixel distance exceeds Y direction windows size.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
			fprintf(stderr,"Error: Specified inter-pixel distance: %d exceeds Y direction windows size: %d !\n", nSpat[i], nWinYSize); 
			exit(1);
		}
	}


	//bin the source image to the GLCM quantization levels
	float dbQBin = (float) ( qGray - 1 ) / ( maxGrayVal - minGrayVal );
	for ( i = 0; i < nImgHeight; i++ )
		for (j = 0; j < nImgWidth; j++ )
			img[i][j]= (int) ( dbQBin * ( img[i][j] - minGrayVal ) );


	for ( m = 0; m < nDist; m++ ){
		//extract GLCM features
		TextureMeasures(ftrs, img, nImgWidth, nImgHeight, qGray, nWinXSize, nWinYSize, nSpat[m*2], nSpat[m*2+1], usemsr, nMsr, bdirInv);

		for ( k = 0; k < nMsr; k++ ){
			//copy ftrs to pRet.  Normalize if signified
			for ( i = 0; i < feat_height; i++ ){
				for ( j = 0; j < feat_width; j++ ){
					if (nNormFlag){
						//normalize to [0,1] based on theoretical statistical bounds
						if ( fabs( statLimits[usemsr[k]][1] - statLimits[usemsr[k]][0] ) > ZERO_THRESH ){
							(*pRet[m*nMsr + k])(i, j) = (ftrs[i*feat_width + j][k] - statLimits[usemsr[k]][0]) / (statLimits[usemsr[k]][1] - statLimits[usemsr[k]][0]);
						}else{
							(*pRet[m*nMsr + k])(i, j) = 0;
						}
					}else{
						(*pRet[m*nMsr + k])(i, j) = ftrs[i*feat_width + j][k];
					}
				}
			}

			if (!this->m_bSaveGlictr)
				continue;

			char band_descriptor[MAXFILENAMELEN];
			char full_file_path[MAXFILENAMELEN];
			if (bdirInv)
				sprintf(band_descriptor, "GLCM_meas=%s_step=%i_winX=%d_winY=%d_qnt=%d_srcimg=%s", statistics[k], nSpat[m*2], nWinXSize, nWinYSize, qGray, source_band_name);
			else
				sprintf(band_descriptor, "GLCM_meas=%s_stepX=%i_stepY=%i_winX=%d_winY=%d_qnt=%d_srcimg=%s", statistics[k], nSpat[m*2], nSpat[m*2+1], nWinXSize, nWinYSize, qGray, source_band_name);

			//if (CanShowStatus())//send filename to GUI
			//	DisplayStatus(band_descriptor,-1);

			//push result to .bil image stack
			images.push_back(pRet[m*nMsr + k]);					
			string bandname(band_descriptor);
			bandnames.push_back(bandname);

			//TODO: Remove these as comment
			//These lines were commented out by Michael due to the how large the csv files were and since they aren't needed for
			//ground truthing and there isn't much space on the D:/ drive, I had to temporarily get rid of them. If they haven't been already,
			//they can be commented out.

			//save to csv file
			//sprintf(full_file_path,"%s%s.csv", m_SaveLocation, band_descriptor);
			//pRet[m*nMsr + k]->WritePixelValuesToFile(full_file_path);


			//normalize to [0,255] based on image max,min.  Then save as bmp.
			float normmin= 99999999999;
			float normmax=-99999999999;
			for (i = 0; i < feat_height*feat_width; i ++ ){
				if (normmin > ftrs[i][k])
					normmin = ftrs[i][k];
				if (normmax < ftrs[i][k])
					normmax = ftrs[i][k];
			}

			for ( i = 0; i < feat_height; i++ ){
				for ( j = 0; j < feat_width; j++ ){
					if (normmax-normmin > ZERO_THRESH){
						float zero_to_one = (ftrs[i*feat_width + j][k] - normmin) / (normmax-normmin);
						(*bmpImg)(i,j) = (int) (255*zero_to_one + 0.5);
					}else{
						(*bmpImg)(i,j) = 0;
					}
				}
			}

			sprintf(full_file_path,"%.230s%.20s.bmp", m_SaveLocation, band_descriptor);
			// todo: fix this
//			bmpImg->Save(full_file_path);
		}
	}

	if (this->m_bSaveGlictr){
//		//save all feature bands to multiband bil file as 32 bit floats
//		MultivariateReadWrite<FLOAT> mrw;
//		char featimgfile_bil[MAXFILENAMELEN];
//		time_t now = time(0);  // get time now
//		struct tm  tstruct;
//		tstruct = *localtime(&now);
//		char time_buf[200];
//		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d_%H.%M.%S", &tstruct);
//		sprintf(featimgfile_bil, "%sGLCMBank_%s.bil", m_SaveLocation, time_buf);
//		mrw.Write(string(featimgfile_bil), images, ImageHeaderInfo::IMG_FLOAT32, bandnames);

	}
	delete bmpImg;

	// Free memory;
	delete pImg;
	for ( i = 0; i < numSamps; i++ )  //this takes super long
		if (ftrs[i]) delete[] ftrs[i];
	if (ftrs) delete[] ftrs;

	return pRet;
}


// Calculates the outputs of the Gabor filter bank, which is designed by giving the range of the spatial-frequency
// (Ul and Uh) and the total number of scales (nScaleNum) and orientations (nOrientNum) used to partition the spectrum.
// Param in:

// Return:
TGrayImage<FLOAT>** TTextureFeatExtract :: GaborFltrBank(BOOL bRealOnlyFlag, TGrayImage<FLOAT>* pSrcImg, int side, double Ul, double Uh, int nScaleNum, int nOrientNum, BOOL bDCRemFlag, int nNormFlag) //The number of output feature images equals to nScaleNum*nOrientNum
{
	BOOL bFeatImgOut = 0;
	int i, j, k, m;
	Matrix *img, *F_r, *F_i;
	float dbMax, dbMin;
	int nImgHeight, nImgWidth;

	std::vector<TGrayImage<FLOAT>*> images;
	std::vector<string> bandnames;

	nImgHeight = pSrcImg->Height();
	nImgWidth = pSrcImg->Width();

	// Dynamic memory allocation
	CreateMatrix(&img, nImgHeight, nImgWidth);
	for ( i = 0; i < nImgHeight; i++ )
		for ( j = 0; j < nImgWidth; j++ ) 
			img->data[i][j] = (double) ( (*pSrcImg)(i, j) );

	CreateMatrix(&F_r, nImgHeight*nScaleNum, nImgWidth*nOrientNum);
	CreateMatrix(&F_i, nImgHeight*nScaleNum, nImgWidth*nOrientNum);

	TGrayImage<FLOAT>** pRet = new TGrayImage<FLOAT>*[nScaleNum*nOrientNum];
	for ( i = 0; i < nScaleNum*nOrientNum; i++ )
		pRet[i] = new TGrayImage<FLOAT>(nImgWidth, nImgHeight);

	GaborFilteredImg(F_r, F_i, img, side, Ul, Uh, nScaleNum, nOrientNum, bDCRemFlag);

	dbMax = -100000000;
	dbMin = 100000000;
	for ( i = 0; i < nScaleNum*nImgHeight; i++ ){
		for ( j = 0; j < nOrientNum*nImgWidth; j++ ){
			if ( bRealOnlyFlag ){
				if ( dbMax < fabs( (float) F_r->data[i][j]) )
					dbMax = fabs( (float) F_r->data[i][j] );
				if ( dbMin > fabs( (float) F_r->data[i][j]) )
					dbMin = fabs( (float) F_r->data[i][j] );	
			}
			else{
				float tmp = (float) sqrt( F_r->data[i][j]*F_r->data[i][j] + F_i->data[i][j]*F_i->data[i][j] );
				if ( dbMax < tmp )
					dbMax = tmp;
				if ( dbMin > tmp )
					dbMin = tmp;
			}
		}
	}


	for ( i = 0; i < nScaleNum; i++ ){
		for ( j = 0; j < nOrientNum; j++ ){
			if ( nNormFlag == 1 ){
				for ( k = 0; k < nImgHeight; k++ ){
					for ( m = 0; m < nImgWidth; m++ ){
						if ( bRealOnlyFlag ){
							if ( fabs( dbMax ) > ZERO_THRESH )
								(*pRet[i*nOrientNum + j])(k, m) =  fabs( (float) F_r->data[i*nImgHeight + k][j*nImgWidth + m]) / dbMax; // Real part only
							else
								(*pRet[i*nOrientNum + j])(k, m) = 0;
						}
						else{
							if ( fabs( dbMax ) > ZERO_THRESH )
								(*pRet[i*nOrientNum + j])(k, m) = ( (float) sqrt( ( F_r->data[i*nImgHeight + k][j*nImgWidth + m] )*( F_r->data[i*nImgHeight + k][j*nImgWidth + m] ) + ( F_i->data[i*nImgHeight + k][j*nImgWidth + m] )*( F_i->data[i*nImgHeight + k][j*nImgWidth + m] ) ) ) / dbMax;
							else
								(*pRet[i*nOrientNum + j])(k, m) = 0;
						}
					}
				}
			}
			else if ( nNormFlag == 2 ){
				for ( k = 0; k < nImgHeight; k++ ){
					for ( m = 0; m < nImgWidth; m++ ){
						if ( bRealOnlyFlag ){
							if ( fabs( dbMax - dbMin ) > ZERO_THRESH )
								(*pRet[i*nOrientNum + j])(k, m) = ( fabs( (float) F_r->data[i*nImgHeight + k][j*nImgWidth + m] ) - dbMin ) / ( dbMax - dbMin ); // Real part only
							else
								(*pRet[i*nOrientNum + j])(k, m) = 0;
						}
						else{
							if ( fabs( dbMax - dbMin ) > ZERO_THRESH )
								(*pRet[i*nOrientNum + j])(k, m) = ( ( (float) sqrt( ( F_r->data[i*nImgHeight + k][j*nImgWidth + m] )*( F_r->data[i*nImgHeight + k][j*nImgWidth + m] ) + ( F_i->data[i*nImgHeight + k][j*nImgWidth + m] )*( F_i->data[i*nImgHeight + k][j*nImgWidth + m] ) ) ) - dbMin ) / ( dbMax - dbMin );
							else
								(*pRet[i*nOrientNum + j])(k, m) = 0;
						}
					}
				}
			}
			else{
				for ( k = 0; k < nImgHeight; k++ ){
					for ( m = 0; m < nImgWidth; m++ ){
						if ( bRealOnlyFlag )
							(*pRet[i*nOrientNum + j])(k, m) = fabs( (float) F_r->data[i*nImgHeight + k][j*nImgWidth + m] ); // Real part only
						else
							(*pRet[i*nOrientNum + j])(k, m) = (float) sqrt( ( F_r->data[i*nImgHeight + k][j*nImgWidth + m] )*( F_r->data[i*nImgHeight + k][j*nImgWidth + m] ) + ( F_i->data[i*nImgHeight + k][j*nImgWidth + m] )*( F_i->data[i*nImgHeight + k][j*nImgWidth + m] ) );
					}
				}
			}

			if (!bFeatImgOut)
				continue;


			//char tbuf[20];
			TGrayImage<int>* pFeatImg = new TGrayImage<int>( nImgWidth, nImgHeight );

			if (nNormFlag){		
				for ( k = 0; k < nImgHeight; k++ )
					for ( m = 0; m < nImgWidth; m++ )
						(*pFeatImg)(k, m) = (int) ( 255*(*pRet[i*nOrientNum + j])(k, m) + 0.5 );
			}
			else{	
				double dbMaxVal = -100000000, dbMinVal = 100000000;
				for ( k = 0; k < nImgHeight; k++ ){
					for ( m = 0; m < nImgWidth; m++ ){
						if ( dbMaxVal < (*pRet[i*nOrientNum + j])(k, m) )
							dbMaxVal = (*pRet[i*nOrientNum + j])(k, m);
						if ( dbMinVal > (*pRet[i*nOrientNum + j])(k, m) )
							dbMinVal = (*pRet[i*nOrientNum + j])(k, m);
					}
				}

				for ( k = 0; k < nImgHeight; k++ ){
					for ( m = 0; m < nImgWidth; m++ ){
						if ( fabs( dbMaxVal - dbMinVal ) > ZERO_THRESH )
							(*pFeatImg)(k, m) = (int) ( 255 * ( (*pRet[i*nOrientNum + j])(k, m) - dbMinVal ) / ( dbMaxVal - dbMinVal ) + 0.5 );
						else
							(*pFeatImg)(k, m) = (int) 0;
					}
				}
			}

			//save each Gabor filter result as separate bmp file
			char featimgfile_f[MAXFILENAMELEN] = "featimg";
			sprintf(featimgfile_f,"Scale%dOrient%d.bmp",i,j);
			// pFeatImg->Save( featimgfile_f );
			delete pFeatImg;

			//format feature outputs for output to a single .bil file
			images.push_back(pRet[i*nOrientNum + j]);
			char featimgbandname[MAXFILENAMELEN];
			sprintf(featimgbandname,"Scale%dOrient%d.bmp",i,j);
			bandnames.push_back(string(featimgbandname));		

		}
	}
	if (bFeatImgOut){
//		//save all feature bands to multiband bil file as 32 bit floats
//		MultivariateReadWrite<FLOAT> mrw;
//		char featimgfile_bil[MAXFILENAMELEN];
//		sprintf(featimgfile_bil, "Gabor.bil");
//		mrw.Write(string(featimgfile_bil), images, ImageHeaderInfo::IMG_FLOAT32, bandnames);
	}

	// Free Memory
	FreeMatrix(img);
	FreeMatrix(F_r);
	FreeMatrix(F_i);

	return pRet;
}


//Much of what goes on in this function is based off this paper: 
//Clausi, David A., and M. E. M. E. Jernigan. "Designing Gabor filters for optimal texture separability." Pattern Recognition 33.11 (2000): 1835-1849.

/*
TGrayImage<FLOAT>** TTextureFeatExtract :: GaborBovik(TGrayImage<FLOAT> *pSrcImg, bool saveImages, BOOL bRealOnlyFlag, BOOL bPostFltFlag, double dbPostLamda, 
	double low_rad, double high_rad, double inc_rad, BOOL bLinSpaceFlag, double low_ang, double high_ang, double inc_ang, 
	BOOL bUnitLamdaFlag, BOOL bDCRemFlag, int &nGaborCnt, int nNormFlag, System::Windows::Forms::ProgressBar^ progressbar, const char* source_band_name)
*/

TGrayImage<FLOAT>** TTextureFeatExtract :: GaborBovik(TGrayImage<FLOAT> *pSrcImg, bool saveImages, BOOL bRealOnlyFlag, BOOL bPostFltFlag, double dbPostLamda, 
	double low_rad, double high_rad, double inc_rad, BOOL bLinSpaceFlag, double low_ang, double high_ang, double inc_ang, 
	BOOL bUnitLamdaFlag, BOOL bDCRemFlag, int &nGaborCnt, int nNormFlag, const char* source_band_name) // int &nGaborCnt,
{

	BOOL bFeatImgOut = 0;

	int	row, col;	    /* image row and column indices */
	int nRows, nCols;   /* number of input image rows and cols*/
	int nNewRows, nNewCols;   /* number of image rows and cols after padding before FFT*/
	int pRows, pCols;   /* padded rows and cols*/
	int tpRows, tpCols; /* temperary variables */

	float	**imgR;	    /* buffer for input images*/
	float	**imgI;
	float	**filterR;	/* buffer for filter responses */
	float	**filterI;

	double decibels =  6.0;
	double	lambda, lambda2;	/* aspect ratio; used to modify sigma in x direction */
	double	sigma, sigma2;	/* standard deviation of Gaussian in spatial domain */
	double	filter_val;	/* filter value at a particular point */
	double	alpha, gain;
	double band_rad, band_ang;
	double	u,v;
	double	up, vp;//, ux, uy, vx, vy, uxy,vxy, uxp, vxp, uyp, vyp, uxyp, vxyp;
	double	theta, freq;
	int	num_rad;	/* number of filter radial center frequencies */
	int	num_ang;	/* number of filter angles */
	double 	ag, rd, ratio;
	int	count = 0;
	int	i, crows, ccols;
	double	costheta, sintheta, Pi, coeff;

	float dbMax, dbMin; /* maximum and minimum filter response for normalization purpose */

	//for saving multiband .bil file
	std::vector<TGrayImage<FLOAT>*> images;
	std::vector<string> bandnames;

	nRows = pSrcImg->Height();
	nCols = pSrcImg->Width();
	nNewRows = nRows;
	nNewCols = nCols;

	/* create image buffers */
	imgR = (float **) malloc( nRows*sizeof( float* ) );
	if ( imgR == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}
	imgR[0] = (float *) calloc( nRows*nCols, sizeof( float ) );
	if ( imgR[0] == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}
	imgI = (float **) malloc( nRows*sizeof( float* ) );
	if ( imgI == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}
	imgI[0] = (float *) calloc( nRows*nCols, sizeof( float ) );
	if ( imgI[0] == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}

	for ( i = 1; i < nRows; i++ )
	{
		imgR[i] = imgR[i-1] + nCols;
		imgI[i] = imgI[i-1] + nCols;
	}

	for ( row = 0; row < nRows; row++ )
		for (col = 0; col < nCols; col++ )
			imgR[row][col] = (*pSrcImg)(row, col);  // Here, the maximum and minimum pixel values are 255 and 0

	fft(imgR, imgI, nNewRows, nNewCols, pRows, pCols);

	// Initialize Gabor filter bank parameters
	gain = pow(10.0,-decibels/20.0);
	alpha = sqrt(log(1.0/gain)/2.0);

	// Default value settings
	if ( low_rad < 0.0 ) low_rad = 1.0; 
	if ( high_rad < 0.0 ) high_rad = log(nNewRows/2.0) / log(2.0); 
	if ( inc_rad < 0.0 ) inc_rad = 1.0;

	/* if high_rad is smaller than low_rad, end program */
	/* if high_rad is greater than the image size, end program */
	/* if low_rad is less than 1.0, end program */
	if ( high_rad < low_rad ) 
	{
		//MessageBox::Show("High radius needs to be less than low radius.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf(stderr,"Incorrect parameter selection (high_rad, low_rad).\n");
		fprintf(stderr,"Program terminated.\n");
		exit(1);
	}

	if ( bLinSpaceFlag && ( high_rad > ((double) (nNewRows/2.0)) ) || (!bLinSpaceFlag) && ( high_rad > (double) (log(nNewRows/2.0)/log(2.0))) ) { //QIN_WHY: (double) log(nNewRows/2.0)/log(2.0) 
		//MessageBox::Show("Highest radial centre frequency greater than the max frequency.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf(stderr,"Highest radial centre frequency (%lf) greater than the max frequency.\n", high_rad);
		fprintf(stderr,"Program terminated.\n");
		exit(1);
	}

	if ( low_rad < 1.0 ) {
		//MessageBox::Show("Lowest radial centre frequency cannot be less than 1.0.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf(stderr,"Lowest radial centre frequency cannot be less than 1.0.\n");
		fprintf(stderr,"Program terminated.\n");
		exit(0);
	}

	num_rad = (int)((high_rad - low_rad)/inc_rad + 1.0);

	// Default value settings
	if ( low_ang < 0.0 || high_ang < 0.0 || inc_ang < 0.0){
		low_ang = 0; 
		high_ang = 179; 
		inc_ang = 30;
	}


	/* orientation increment must be a positive number less than 180 degrees */
	if (inc_ang > 180.0 || inc_ang <= 0.0)  {
		//MessageBox::Show("Angle increment must be a positive number less than 180 degrees.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf(stderr,"inc_ang must be a positive number less than 180 degrees.\n");
		fprintf(stderr,"Program terminated.\n");
		exit(1);
	}

	if ( ( high_ang < low_ang ) ) {
		//MessageBox::Show("High angle must be less than low angle.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf(stderr,"Incorrect parameter selection (high_ang, low_ang).\n");
		fprintf(stderr,"Program terminated.\n");
		exit(1);
	}

	num_ang = (int)((high_ang - low_ang)/inc_ang + 1.0);

	fprintf(stderr,"Number of radii and angular bandwidths selected: %5d %5d\n", num_rad, num_ang);

	if (bLinSpaceFlag)
		fprintf(stderr,"Selected filter radii (cpi): ");
	else
		fprintf(stderr,"Selected filter radii (oct): ");
	for (i = 0; i < num_rad; ++i) fprintf(stderr,"%5.1f ",low_rad + i*inc_rad);
	fprintf(stderr,"\n");

	fprintf(stderr,"Selected filter orientations (deg): ");
	for (i = 0;i < num_ang; ++i) fprintf(stderr,"%5.1f ", low_ang + i*inc_ang);
	fprintf(stderr,"\n");

	band_rad = inc_rad;
	ratio = ( pow( 2.0, band_rad ) - 1.0 ) / ( pow( 2.0, band_rad ) + 1.0 );
	Pi = 4.0*atan(1.0); // the returned value of atan() is radians
	coeff = -2.0*Pi*Pi;

	if (bUnitLamdaFlag) 
	{
		lambda = 1.0; lambda2 = 1.0;

		/* keep band_ang in degrees */
		band_ang = 180.0/Pi*2.0*atan(ratio);
	}
	else 
		band_ang = inc_ang;

	if (bLinSpaceFlag)
		fprintf(stderr,"Filter  Angle   Freq   Sigma  Lambda  RadBand(cpi)  AngBand\n");
	else
		fprintf(stderr,"Filter  Angle   Freq   Sigma  Lambda  RadBand(oct)  AngBand\n");

	crows = nNewRows/2; ccols = nNewCols/2;

	if (bDCRemFlag) 
		imgR[crows][ccols] = 0.0;

	nGaborCnt = num_rad*num_ang;

	/* create result buffers */
	filterR = (float **) malloc( nNewRows*sizeof( float* ) );
	if ( filterR == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}
	filterR[0] = (float *) calloc( nNewRows*nNewCols, sizeof( float ) );
	if ( filterR[0] == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}
	filterI = (float **) malloc( nNewRows*sizeof( float* ) );
	if ( filterI == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}
	filterI[0] = (float *) calloc( nNewRows*nNewCols, sizeof( float ) );
	if ( filterI[0] == NULL ) 
	{
		//MessageBox::Show("Out of memory.  Terminating.", "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
		fprintf( stderr, "Fatal Error - out of memory.\n");
		exit(-1);
	}
	for ( i = 1; i < nNewRows; i++ )
	{
		filterR[i] = filterR[i-1] + nNewCols;
		filterI[i] = filterI[i-1] + nNewCols;
	}

	/* create return images*/
	TGrayImage<FLOAT>** pRet = new TGrayImage<FLOAT>*[nGaborCnt];
	for ( i = 0; i < nGaborCnt; i++ )
		pRet[i] = new TGrayImage<FLOAT>(nCols, nRows);

	for ( rd = low_rad; rd <= high_rad; rd += inc_rad ){ 
		for ( ag = low_ang; ag <= high_ang; ag += inc_ang ){

			theta = ag*Pi/180.0;
			costheta = cos(theta); 
			sintheta = sin(theta);

			if (bLinSpaceFlag)
				freq = rd;
			else
				freq=pow( 2.0, (rd - 1.0) )*sqrt(2.0);

			sigma = alpha/(Pi*freq*tan(band_ang/2.0*Pi/180.0));
			sigma2 = sigma*sigma;


			if (!bUnitLamdaFlag){
				lambda = alpha/(Pi*freq*sigma*ratio); 
				lambda2 = lambda*lambda;
			}


			//Apply Gabor filter in spatial-frequency domain
			for( row = 0; row < nNewRows; row++ ){ 
				for ( col = 0; col < nNewCols; col++ ){
					u = (double) (col - ccols);
					v = (double) (crows - row);
					up = u*costheta + v*sintheta;
					vp = -u*sintheta + v*costheta;

					filter_val = GaborFreq(up,vp, sigma2, freq, lambda2, coeff);

					filterR[row][col] = imgR[row][col]*filter_val;
					filterI[row][col] = imgI[row][col]*filter_val;
				}
			}
			ifft(filterR, filterI, nNewRows, nNewCols);



			if (!bRealOnlyFlag){
				//store magnitude of complex response in filterR
				rec2circ(filterR, filterI, nNewRows, nNewCols);
				for( row = 0; row < nNewRows; row++ ) 
					for ( col = 0; col < nNewCols; col++ ) 
						filterI[row][col] = 0.0;
			}
			else{
				for( row = 0; row < nNewRows; row++ ) 
					for ( col = 0; col < nNewCols; col++ ){
						filterR[row][col] = fabs( filterR[row][col] );
						filterI[row][col] = 0.0;
					}
			}	


			//Apply post filter smoothing in spatial-frequency domain
			if (bPostFltFlag){
				fft(filterR, filterI, nNewRows, nNewCols, tpRows, tpCols);

				for( row = 0; row < nNewRows; row++ ){ 
					for ( col = 0; col < nNewCols; col++ ){
						u = (double) (col - ccols);
						v = (double) (crows - row);
						up = u*costheta + v*sintheta;
						vp = -u*sintheta + v*costheta;

						filter_val = GaborPost(up, vp, sigma2, lambda2, coeff, dbPostLamda);

						filterR[row][col] *= filter_val;
						filterI[row][col] *= filter_val;
					}
				}

				ifft(filterR, filterI, nNewRows, nNewCols);
				rec2circ(filterR, filterI, nNewRows, nNewCols);
			}


			//copy filter response to return image
			for( row = pRows; row < nRows + pRows; row++ ) 
				for ( col = pCols; col < nCols + pCols; col++ ) 
					(*pRet[count])(row - pRows, col - pCols) = filterR[row][col];

			count++;
			fprintf(stderr," %3d   %5.1lf   %5.2lf  %7.6lf  %6.3lf  %6.3lf   %6.3lf\n", count, ag, freq, sigma ,lambda, band_rad, band_ang);
		}
	}

	//progressbar->Value = 30;

	free(filterR[0]);
	free(filterR);
	free(filterI[0]);
	free(filterI);
	free(imgR[0]);
	free(imgR);
	free(imgI[0]);
	free(imgI);

	dbMax = -100000000;
	dbMin = 100000000;
	for ( i = 0; i < count; i++ ){
		for ( row = 0; row < nRows; row++ ){
			for ( col = 0; col < nCols; col++ ){
				if ( dbMax < (*pRet[i])(row, col) )
					dbMax = (*pRet[i])(row, col);
				if ( dbMin > (*pRet[i])(row, col) )
					dbMin = (*pRet[i])(row, col);	
			}
		}
	}

	//progressbar->Step = 70 / num_rad / num_ang + 1;

	for ( i = 0; i < num_rad; i++ ){
		for ( int j = 0; j < num_ang; j++ ){
			//progressbar->PerformStep();
			if ( nNormFlag == 1 ){  //divide all responses by max response
				for ( row = 0; row < nRows; row++ ){
					for ( col = 0; col < nCols; col++ ){
						if ( fabs( dbMax ) > ZERO_THRESH )
							(*pRet[i*num_ang + j])(row, col) =  (*pRet[i*num_ang + j])(row, col) / dbMax;
						else
							(*pRet[i*num_ang + j])(row, col) = 0;
					}
				}
			}
			else if ( nNormFlag == 2 ){  //scale responses to between 0 and 1
				for ( row = 0; row < nRows; row++ ){
					for ( col = 0; col < nCols; col++ ){
						if ( fabs( dbMax - dbMin ) > ZERO_THRESH )
							(*pRet[i*num_ang + j])(row, col) = ( (*pRet[i*num_ang + j])(row, col) - dbMin ) / ( dbMax - dbMin );
						else
							(*pRet[i*num_ang + j])(row, col) = 0;
					}
				}
			}


			//if ( bFeatImgOut ) this is flag is replced by global flag m_bSaveGabor by sharif
			if (true)
				continue;

			char featimgfile_f[MAXFILENAMELEN];// = "Gabr_featimg";
			char band_descriptor[MAXFILENAMELEN];
			sprintf(band_descriptor, "Gabr_freqoctIDX=%d_angIDX=%d_postlambda=%2.2f_srcimg=%s", i, j, dbPostLamda, source_band_name);

			//Michael - added boolean saveImages, since we do not want to save the images to disk when we do feature sampling for feature selection,
			//this gives us the option to only generate the images within MAGIC
			if (saveImages)
			{
				TGrayImage<int>* pFeatImg = new TGrayImage<int>(nCols, nRows);

				if (nNormFlag) {
					//scale to (0,255)
					for (row = 0; row < nRows; row++)
						for (col = 0; col < nCols; col++)
							(*pFeatImg)(row, col) = (int)(255 * (*pRet[i*num_ang + j])(row, col) + 0.5);
				}
				else {
					double dbMaxVal = -100000000, dbMinVal = 100000000;
					for (row = 0; row < nRows; row++) {
						for (col = 0; col < nCols; col++) {
							if (dbMaxVal < (*pRet[i*num_ang + j])(row, col))
								dbMaxVal = (*pRet[i*num_ang + j])(row, col);
							if (dbMinVal > (*pRet[i*num_ang + j])(row, col))
								dbMinVal = (*pRet[i*num_ang + j])(row, col);
						}
					}

					//scale to (0,255)
					for (row = 0; row < nRows; row++) {
						for (col = 0; col < nCols; col++) {
							if (fabs(dbMaxVal - dbMinVal) > ZERO_THRESH)
								(*pFeatImg)(row, col) = (int)(255 * ((*pRet[i*num_ang + j])(row, col) - dbMinVal) / (dbMaxVal - dbMinVal) + 0.5);
							else
								(*pFeatImg)(row, col) = (int)0;
						}
					}
				}

				sprintf(featimgfile_f, "%.230s%.20s.bmp", m_SaveLocation, band_descriptor);
				// pFeatImg->Save(featimgfile_f);
				//if (CanShowStatus())// SHARIF, send filename to GUI
					//DisplayStatus(band_descriptor, -1);
				delete pFeatImg;
			}

			//save csv file
			if (saveImages)
			{
				sprintf(featimgfile_f, "%.230s%.20s.csv", m_SaveLocation, band_descriptor);
				pRet[i*num_ang + j]->WritePixelValuesToFile(featimgfile_f);  //save to txt file
			}
			//format feature outputs for output to a single .bil file
			images.push_back(pRet[i*num_ang + j]);
			bandnames.push_back(string(band_descriptor));			
		}
	}

	if (this->m_bSaveGabor && saveImages){
//		//save all feature bands to multiband bil file as 32 bit floats
//		MultivariateReadWrite<FLOAT> mrw;
//		char featimgfile_bil[MAXFILENAMELEN];
//		time_t now = time(0);  // get time now
//		struct tm  tstruct;
//		tstruct = *localtime(&now);
//		char time_buf[200];
//		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d_%H.%M.%S", &tstruct);
//		sprintf(featimgfile_bil, "%sGaborBank_%s.bil", m_SaveLocation, time_buf);
//		mrw.Write(string(featimgfile_bil), images, ImageHeaderInfo::IMG_FLOAT32, bandnames);
	}

	//progressbar->Value = 0;
	//progressbar->Step = 1;

	return pRet;
}

/***************** Operations on univariate images ******************/
/********************************************************************/


/************************* Intermediary Operations **************************/

/*** TextureMeasures : Calculates Haralick features using the faster iterative method ***/
//
//Param in: ftrsout - extracted feature measures
//			img     - input quantized image array
//			nXSize  - image width
//			nYSize  - image height
//			gnum    - # of quantized grey levels in image
//          nFlszX  - desired windows width.  [1 to image width]
//          nFlszY  - desired windows height.  [1 to image height]
//          nSpatX  - interpixel space of delta x columns. [1 to nFlszX-1]
//                  - interpixel distance if direction invariant
//          nSpatY  - interpixel space of delta y rows. [1 to nFlszY-1]
//          usemsr  - list of desired features to calculate.  Pass as an array
//                   of numbers according to the enum
//          nMeas   - number of feature measures to output
//          bDirInv - flag for directional invariance
//                  - 1- true, 0-false
//
//Param out: feature vectors for all pixels of the image
/****************************************************************************************/
void TTextureFeatExtract :: TextureMeasures(float** ftrsout, int** img, int nXSize, int nYSize, int gnum, 
	int nFlszX, int nFlszY, int nSpatX, int nSpatY, 
	int* usemsr, int nMeas, int bDirInv)
{
	int numpairs; //number of pixel pairs in window, single direction
	int doublepairs; //number of pixel pairs in window, symmetrical directions
	int nFXHalf, nFYHalf; // half of the window size
	// bools indicating whether to calculate the corresponding feature
	int msrbool[NUM_TEX_MEASURES]={0};    
	int i, j, k, lkindex, lkindex2, row, col, line, pixel; //counters
	int iOutIndex;
	int ftroutindex;
	int totalrows;

	//dynamic memory
	int** ppbyData; //current strip of image being processed by the windows
	float** pfOutBuf; //array of feature vectors for current row
	int** glcm; //the grey level co-ocurrence matrix
	float* tempmeas; //temporary statistics for COR
	float* measures; //array of current measures data for current window
	float** lktable[NUM_TEX_MEASURES] = {NULL}; //array of pointers to lookup tables

	/////////////////////// init parameters
	// map list of desired features to array of bools
	for (i = 0; i < nMeas; i++)
		msrbool[usemsr[i]] = 1;

	//the following measures are dependent on others being calculated as well
	if (msrbool[COR])
		msrbool[STD] = 1;

	if (msrbool[STD])
		msrbool[MU] = 1;

	/////////////////////// calculate the number of pixel pairs in the window
	if (bDirInv) //sum of #of pairs for each orientation
		numpairs = (nFlszX - abs(nSpatX)) * (nFlszY - abs(0)) +
		(nFlszX - abs(nSpatX)) * (nFlszY - abs(nSpatX)) + 
		(nFlszX - abs(0)) * (nFlszY - abs(nSpatX)) + 
		(nFlszX - abs(-nSpatX)) * (nFlszY - abs(nSpatX));
	else
		numpairs = (nFlszX - abs(nSpatX)) * (nFlszY - abs(nSpatY));

	doublepairs = numpairs * 2;

	/////////////////////// Allocate Memory
	if (msrbool[ASM] || msrbool[ENT]) //init the matrix if one is needed
	{
		glcm = (int**) calloc( gnum, sizeof(int*) );
		if ( !( glcm ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate glcm\n" ); 
			exit(1);
		}
		for ( i = 0; i < gnum; i++ )
			glcm[i] = (int*) calloc( gnum, sizeof(int) );
	}
	else
		glcm = NULL;

	ppbyData = (int**) calloc( nFlszY, sizeof(int*) );
	if ( !( ppbyData ) ) 
	{
		fprintf( stderr, "Not enough memory to allocate ppbyData\n" ); 
		exit(1);
	}
	for ( i = 0; i < nFlszY; i++ )
		ppbyData[i] = (int*) calloc( nXSize, sizeof(int) );

	pfOutBuf = (float**) calloc( nXSize - nFlszX + 1, sizeof(float*) );
	if ( !( pfOutBuf ) ) 
	{
		fprintf( stderr, "Not enough memory to allocate pfOutBuf\n" ); 
		exit(1);
	}
	for ( i = 0; i < nXSize - nFlszX + 1; i++ )
		pfOutBuf[i] = (float*) calloc( nMeas, sizeof(float) );

	measures = (float*) calloc( NUM_TEX_MEASURES, sizeof(float) );
	tempmeas = (float*) calloc( NUM_TEMP_MEASURES, sizeof(float) );

	// allocate lookup tables
	if (msrbool[ASM])
	{
		lktable[ASM] = (float**) calloc( doublepairs + 1, sizeof(float*) );
		if ( !( lktable[ASM] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[ASM]\n" ); 
			exit(1);
		}
		for ( i = 0; i <= doublepairs; i++ )
			lktable[ASM][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[ENT])
	{
		lktable[ENT] = (float**) calloc( doublepairs + 1, sizeof(float*) );
		if ( !( lktable[ENT] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[ENT]\n" ); 
			exit(1);
		}
		for ( i = 0; i <= doublepairs; i++ )
			lktable[ENT][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[DIS])
	{
		lktable[DIS] = (float**) calloc( gnum, sizeof(float*) );
		if ( !( lktable[DIS] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[DIS]\n" ); 
			exit(1);
		}
		for ( i = 0; i < gnum; i++ )
			lktable[DIS][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[CON])
	{
		lktable[CON] = (float**) calloc( gnum, sizeof(float*) );
		if ( !( lktable[CON] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[CON]\n" ); 
			exit(1);
		}
		for ( i = 0; i < gnum; i++ )
			lktable[CON][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[HOM])
	{
		lktable[HOM] = (float**) calloc( gnum, sizeof(float*) );
		if ( !( lktable[HOM] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[HOM]\n" ); 
			exit(1);
		}
		for ( i = 0; i < gnum; i++ )
			lktable[HOM][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[INV])
	{
		lktable[INV] = (float**) calloc( gnum, sizeof(float*) );
		if ( !( lktable[INV] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[INV]\n" ); 
			exit(1);
		}
		for ( i = 0; i < gnum; i++ )
			lktable[INV][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[STD])
	{
		lktable[STD] = (float**) calloc( gnum, sizeof(float*) );
		if ( !( lktable[STD] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[STD]\n" ); 
			exit(1);
		}
		for ( i = 0; i < gnum; i++ )
			lktable[STD][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[MU])
	{
		lktable[MU] = (float**) calloc( 2*(gnum - 1) + 1, sizeof(float*) );
		if ( !( lktable[MU] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[MU]\n" ); 
			exit(1);
		}
		for ( i = 0; i <= 2*(gnum - 1); i++ )
			lktable[MU][i] = (float*) calloc( 1, sizeof(float) );
	}
	if (msrbool[COR])
	{
		lktable[COR] = (float**) calloc( gnum, sizeof(float*) );
		if ( !( lktable[COR] ) ) 
		{
			fprintf( stderr, "Not enough memory to allocate lktable[COR]\n" ); 
			exit(1);
		}
		for ( i = 0; i < gnum; i++ )
			lktable[COR][i] = (float*) calloc( gnum, sizeof(float) );
	}

	////////////////////// generate lookups
	//ASM ENT
	if (msrbool[ENT])
		lktable[ENT][0][0] = 0; // 0 log 0 = undefined = 0

	for (lkindex = 0; lkindex <= doublepairs; lkindex++)
	{
		if (msrbool[ASM]) //up to doublepairs inclusive
			lktable[ASM][lkindex][0] = ((float)lkindex / (float)doublepairs)
			* ((float)lkindex / (float)doublepairs);

		//up to doublepairs inclusive. 0 log 0 = undefined = 0
		if ( msrbool[ENT] & ( lkindex > 0 ) ) 
			lktable[ENT][lkindex][0] = (float)lkindex / (float)doublepairs 
			* log((float)lkindex / (float)doublepairs);
	}

	//DIS CON HOM INV STD
	for (lkindex = 0; lkindex <= gnum-1; lkindex++)
	{
		if (msrbool[CON]) // up to gnum
			lktable[CON][lkindex][0] = ((float)lkindex*lkindex) / numpairs;

		if (msrbool[DIS]) // up to gnum
			lktable[DIS][lkindex][0] = (float)lkindex/numpairs;

		if (msrbool[HOM]) // up to gnum
			lktable[HOM][lkindex][0] = (float) 1 / (numpairs * (1+lkindex*lkindex));

		if (msrbool[INV]) // up to gnum
			lktable[INV][lkindex][0] = (float) 1 / (numpairs * (1+lkindex));

		if (msrbool[STD]) //up to gnum
			lktable[STD][lkindex][0] = ((float)lkindex*lkindex) / doublepairs;
	}

	// MU COR
	if (msrbool[MU]) //up to 2*gnum - 1
		for (lkindex = 0; lkindex <= 2*(gnum-1); lkindex++)
			lktable[MU][lkindex][0] = (float)lkindex / (float)doublepairs;

	if (msrbool[COR]) //up to (gnum,gnum)
	{
		for (lkindex = 0; lkindex <= gnum-1; lkindex++)
			for (lkindex2 = 0; lkindex2 <= gnum-1; lkindex2++)
				lktable[COR][lkindex][lkindex2] = ((float)lkindex*lkindex2) / numpairs;
	}

	////////////////////// Misc Init
	//size of half a window (floored, since nFlsz will be odd)
	nFYHalf = (nFlszY-1)/2; 
	nFXHalf = (nFlszX-1)/2; 

	////////////////////// load first strip of image
	for (row = 0; row < nFlszY; row++)
		for (col = 0; col < nXSize; col++)
			ppbyData[row][col] = img[row][col];

	ftroutindex = 0;
	totalrows = nYSize - 2*nFYHalf;

	fprintf(stderr,"*** Starting ***\n");
	//////////////////////iterate through image, window by window
	for (line = nFYHalf; line < (nYSize - nFYHalf); line++)
	{
		iOutIndex = 0; //current index into the output buffer strip

		//////////////////// do init for first window of the strip
		if (bDirInv)
			build_dirinv(measures, glcm, tempmeas, 
			ppbyData, gnum, nFlszX, nFlszY, nSpatX,
			lktable, msrbool, bDirInv);
		else
			build_single(measures, glcm, tempmeas, 
			ppbyData, gnum, nFlszX, nFlszY, nSpatX, nSpatY, 
			lktable, msrbool, bDirInv);

		//////////////////// move data to the output buffer strip
		moveout(pfOutBuf, &iOutIndex, measures, nMeas, usemsr);
		//for (int i = 0; i < 9; i++)
		//{
		//	measures[i] = 0;
		//}

		//////////////////// update the other windows of the strip
		for (pixel = 1; pixel < (nXSize - nFlszX + 1); pixel++)
		{
			if (bDirInv) //measures[], glcm[][] and tempmeas[] are updated
				update_dirinv(measures, glcm, tempmeas, 
				ppbyData, nFlszX, nFlszY, nSpatX, 
				pixel, lktable, msrbool);
			else
				update_single(measures, glcm, tempmeas, 
				ppbyData, nFlszX, nFlszY, nSpatX, nSpatY, 
				pixel, lktable, msrbool);

			//////////////////// move data to the output buffer strip
			moveout(pfOutBuf, &iOutIndex, measures, nMeas, usemsr);
			//for (int i = 0; i < 9; i++)
			//{
			//	measures[i] = 0;
			//}
		}

		//////////////////// Output Current Line of Vectors To file Here
		// print_mat_flt(pfOutBuf,0,nXSize-nFlszX,0,nMeas-1,"");

		for (k=0; k <= nXSize-nFlszX; k++)
		{
			for (j = 0; j< nMeas; j++)
				ftrsout[ftroutindex][j] = pfOutBuf[k][j];
			ftroutindex++;
		}

		/////////////////// Update the image strip if there are still rows to be loaded
		if (line + nFYHalf + 1 < nYSize)
		{
			for (row = 0; row < (nFlszY - 1); row++)
			{
				for (col = 0; col < nXSize; col++)
					ppbyData[row][col] = ppbyData[row+1][col];
			}
			for (col = 0; col < nXSize; col++) //get next line of pixels
				ppbyData[nFlszY-1][col] = img[line + nFYHalf + 1][col];
		}

		i = totalrows / 5;
		if ( line != nFYHalf && ( line - nFYHalf ) % i == 0 )
			fprintf(stderr, "* %d of %d rows completed *\n", line - nFYHalf + 1,totalrows);  //QIN_MOD
		else if ( line == nYSize - nFYHalf -1)  //QIN_MOD
			fprintf(stderr, "* %d of %d rows completed *\n", line - nFYHalf + 1, totalrows); //QIN_MOD
	}

	//////////////////////////////// Free Memory
	if (msrbool[ASM] || msrbool[ENT]) //free the matrix if one is needed
	{
		for ( i = 0; i <= gnum - 1; i++ )
			free(glcm[i]);
		free(glcm);
	}

	for ( i = 0; i <= nFlszY - 1; i++ )
		free(ppbyData[i]);
	free(ppbyData);

	for ( i = 0; i <= nXSize - nFlszX; i++ )
		free(pfOutBuf[i]);
	free(pfOutBuf);

	free(measures);
	free(tempmeas);

	// free lookup tables
	if (msrbool[ASM])
	{
		for ( i = 0; i <= doublepairs; i++ )
			free(lktable[ASM][i]);
		free(lktable[ASM]);
	}
	if (msrbool[ENT])
	{
		for ( i = 0; i <= doublepairs; i++ )
			free(lktable[ENT][i]);
		free(lktable[ENT]);
	}
	if (msrbool[DIS])
	{
		for ( i = 0; i <= gnum - 1; i++ )
			free(lktable[DIS][i]);
		free(lktable[DIS]);
	}
	if (msrbool[CON])
	{
		for ( i = 0; i <= gnum - 1; i++ )
			free(lktable[CON][i]);
		free(lktable[CON]);
	}
	if (msrbool[HOM])
	{
		for ( i = 0; i <= gnum - 1; i++ )
			free(lktable[HOM][i]);
		free(lktable[HOM]);
	}
	if (msrbool[INV])
	{
		for ( i = 0; i <= gnum - 1; i++ )
			free(lktable[INV][i]);
		free(lktable[INV]);
	}
	if (msrbool[STD])
	{
		for ( i = 0; i <= gnum - 1; i++ )
			free(lktable[STD][i]);
		free(lktable[STD]);
	}
	if (msrbool[MU])
	{
		for ( i = 0; i <= 2*(gnum - 1); i++ )
			free(lktable[MU][i]);
		free(lktable[MU]);
	}
	if (msrbool[COR])
	{
		for ( i = 0; i <= gnum - 1; i++ )
			free(lktable[COR][i]);
		free(lktable[COR]);
	}
}

/******************** moveout : Move measures data to output *********************/
// Moves the measures data to its location in the output buffer
//*****************************************************************************
void TTextureFeatExtract :: moveout( float** pfOutBuf, int* iOutIndex, float* measures, int nMeas, int* usemsr )
{
	for ( int t = 0; t < nMeas; t++ )
		pfOutBuf[*iOutIndex][t] = measures[usemsr[t]];

	*iOutIndex = *iOutIndex + 1;
}

//******** build_dirinv : Initialize row for all orientations ********************
// Cumulatively calls build_single for all orientations: 0, 45, 90 and 135 degrees
// 
// Param in:  same as build_single
// Param out: same as build_single but results are cumulative for all 4 directions
//*****************************************************************************
void TTextureFeatExtract :: build_dirinv( float* measures, int** glcm, float* tstats,
	int** ppbyData, int gnum, int nFlszX, int nFlszY, int nSpat,
	float*** lk, int* msrbool, int bDirInv )
{
	build_single(measures, glcm, tstats, 
		ppbyData, gnum, nFlszX, nFlszY, nSpat, 0, lk, msrbool, bDirInv);
	build_single(measures, glcm, tstats, 
		ppbyData, gnum, nFlszX, nFlszY, nSpat, nSpat, lk, msrbool, bDirInv);
	build_single(measures, glcm, tstats, 
		ppbyData, gnum, nFlszX, nFlszY, 0, nSpat, lk, msrbool,  bDirInv);
	build_single(measures, glcm, tstats, 
		ppbyData, gnum, nFlszX, nFlszY, -nSpat, nSpat, lk, msrbool, bDirInv);
}

//*************** build_single : Initialize row for single orientation ***********
// Calculates measures and glcm for first window of each row from scratch
//
// Param in:  measures - feature vector of the current window to fill
//            glcm     - glcm of the current window to fill
//            tstats   - intermediary stats vector that needs to be carried
//            lk       - lookup tables
//            msrbool  - array of bools indicating which features to calculate
//            bdirInv  - direction invariant flag
// Param out: measures - updated feature vector of the current window
//            glcm     - built glcm of the current window
//            tstats   - updated intermediary stats vector of the current window
//*****************************************************************************
void TTextureFeatExtract :: build_single( float* measures, int** glcm, float* tstats,
	int** ppbyData, int gnum, int nFlszX, int nFlszY, int nSpatX, int nSpatY,
	float*** lk, int* msrbool, int bDirInv )
{
	int i, j;
	int startc, endc, startr, endr;
	int matold1, matnew1; //ASM, ENT variables
	int matold2, matnew2; 
	int dindex;
	int msrindex, tindex, curr, curc; //counters
	float var;

	//reset values only if single or first orientation of init
	if (!(bDirInv) || (nSpatY == 0) )
	{
		for (msrindex = 0; msrindex < NUM_TEX_MEASURES; msrindex++)
			measures[msrindex] = 0;

		for (tindex = 0; tindex < NUM_TEMP_MEASURES; tindex++)
			tstats[tindex] = 0;

		if (glcm != NULL) // only reset if calculating ASM and/or ENT
		{
			for ( i = 0; i <= gnum - 1; i++ )
				for (j = 0; j <= gnum - 1; j++ )
					glcm[i][j] = 0;
		}
	}

	findrange(&startc, &endc, &startr, &endr, nFlszX, nFlszY, nSpatX, nSpatY);

	for (curr = startr; curr < endr; curr++)
	{
		for (curc = startc; curc < endc; curc++)
		{
			//find the current pair of pixels in image
			i = ppbyData[curr][curc];
			j = ppbyData[curr + nSpatY][curc + nSpatX];

			dindex = abs(i - j);

			//initialize statistics for this pair of pixels
			if (msrbool[CON])
				measures[CON] = measures[CON] + lk[CON][dindex][0];         

			if (msrbool[DIS])
				measures[DIS] = measures[DIS] + lk[DIS][dindex][0];

			if (msrbool[HOM])
				measures[HOM] = measures[HOM] + lk[HOM][dindex][0];         

			if (msrbool[INV])
				measures[INV] = measures[INV] + lk[INV][dindex][0];

			if (msrbool[MU])
				measures[MU] = measures[MU] + lk[MU][i + j][0];

			if (msrbool[STD])
				tstats[STDTEMP] = tstats[STDTEMP] + lk[STD][i][0] + lk[STD][j][0];

			if (msrbool[COR])
				tstats[CORTEMP] = tstats[CORTEMP] + lk[COR][i][j];

			//ASM and ENT: considering both (i,j) and (j,i) pairs at once
			if (msrbool[ASM] || msrbool[ENT])
			{
				matold1 = glcm[i][j];
				matnew1 = matold1 + 1;
				glcm[i][j] = matnew1; //update matrix
				matold2 = glcm[j][i];
				matnew2 = matold2 + 1;
				glcm[j][i] = matnew2; //update matrix
			}

			if (msrbool[ASM])
			{
				measures[ASM] -= (lk[ASM][matold1][0] + lk[ASM][matold2][0]);
				if ( abs(measures[ASM]) <= ZERO_THRESH )
					measures[ASM] = lk[ASM][matnew1][0] + lk[ASM][matnew2][0];
				else
					measures[ASM] += (lk[ASM][matnew1][0] + lk[ASM][matnew2][0]);
			}

			//signs switched because ENT is -ve of summation
			if (msrbool[ENT]) //NOTE: This result is x2 the original PCI code and is "correct"
			{
				measures[ENT] += (lk[ENT][matold1][0] + lk[ENT][matold2][0]);
				if ( abs(measures[ENT]) <= ZERO_THRESH )
					measures[ENT] = - lk[ENT][matnew1][0] - lk[ENT][matnew2][0];
				else
					measures[ENT] -= (lk[ENT][matnew1][0] + lk[ENT][matnew2][0]);
			}
		}
	}

	//STD calculations
	if (msrbool[STD])
	{
		var = tstats[STDTEMP] - measures[MU]*measures[MU];
		if ( var <= ZERO_THRESH )  //QIN_MOD
		{
			var = 0;
			measures[STD] = 0;
		}
		else 
			measures[STD] = sqrt(var);
	}

	//COR calculations
	if (msrbool[COR])
	{
		if (measures[STD] <= ZERO_THRESH)
			measures[COR] = 1;
		else
			measures[COR] = (tstats[CORTEMP] - measures[MU]*measures[MU]) / var;
	}
}

/********************* update_dirinv : update for all orientations **********************/
// Cumulatively calls update_single for all orientations: 0, 45, 90 and 135 degrees
//
// Param in:  same as update_single
// Param out: same as update_single but results are cumulative for all 4 directions
/****************************************************************************************/
void TTextureFeatExtract :: update_dirinv(float* measures, int** glcm, float* tstats, 
	int** ppbyData, int nFlszX, int nFlszY, int nSpat,
	int nPixOffset,   float*** lk, int* msrbool)
{
	update_single(measures, glcm, tstats, 
		ppbyData, nFlszX, nFlszY, nSpat, 0, nPixOffset, lk, msrbool);
	update_single(measures, glcm, tstats, 
		ppbyData, nFlszX, nFlszY, nSpat, nSpat, nPixOffset, lk, msrbool);
	update_single(measures, glcm, tstats, 
		ppbyData, nFlszX, nFlszY, 0, nSpat, nPixOffset, lk, msrbool);
	update_single(measures, glcm, tstats, 
		ppbyData, nFlszX, nFlszY, -nSpat, nSpat, nPixOffset, lk, msrbool);
}

/********************** update_single : update for single orientation **************************/
// Param in:  measures   - feature vector for the previous window
//            nPixOffset - current column position of window
//            glcm       - glcm for the previous window
//            tstats     - intermediary stats vector to be carried
//			  ppbyData   - current row of source image
//			  nFlszX     - window width
//            nFlszY     - window height
//            nSpatX     - interpixel space of delta x columns
//            nSpatY     - interpixel space of delta y rows
//            lk         - lookup tables
//            msrbool    - array of bools indicating which features to calculate
//            bdirInv    - direction invariant flag
// Param out: measures   - updated feature vector for the current window
//            glcm       - updated glcm for the current window
//            tstats     - updated intermediary stats vector for the current window
/*****************************************************************************************/
void TTextureFeatExtract :: update_single( float* measures, int** glcm, float* tstats, 
	int** ppbyData, int nFlszX, int nFlszY, int nSpatX, int nSpatY,
	int nPixOffset,   float*** lk, int* msrbool )
{
	int ir, jr, il, jl;
	int lmatold1, lmatnew1; //ASM, ENT variables for left column
	int lmatold2, lmatnew2;
	int rmatold1, rmatnew1; //ASM, ENT variables for right column
	int rmatold2, rmatnew2;
	int startc, endc, startr, endr;
	int dindexl, dindexr;
	int curr; //counters
	float var;

	findrange(&startc, &endc, &startr, &endr, nFlszX, nFlszY, nSpatX, nSpatY);

	for (curr = startr; curr < endr; curr++)
	{
		//old col
		il = ppbyData[curr][nPixOffset + startc - 1];
		jl = ppbyData[curr+nSpatY][nPixOffset + startc - 1 + nSpatX];  //cooccurring pair of point il

		//new col
		ir = ppbyData[curr][nPixOffset + endc - 1];
		jr = ppbyData[curr + nSpatY][nPixOffset + endc - 1 + nSpatX];

		dindexl = abs(il - jl);
		dindexr = abs(ir - jr);

		//calculate statistics: curstat = curstat - oldcol + newcol
		if (msrbool[CON])
		{
			measures[CON] -= lk[CON][dindexl][0];
			if (abs(measures[CON]) <= ZERO_THRESH)
				measures[CON] = lk[CON][dindexr][0];
			else
				measures[CON] += lk[CON][dindexr][0];
		}

		if (msrbool[DIS])
		{
			measures[DIS] -= lk[DIS][dindexl][0];
			if (abs(measures[DIS]) <= ZERO_THRESH)
				measures[DIS] = lk[DIS][dindexr][0];
			else
				measures[DIS] += lk[DIS][dindexr][0];
		}

		if (msrbool[HOM])
		{
			measures[HOM] -= lk[HOM][dindexl][0];
			if (abs(measures[HOM]) <= ZERO_THRESH)
				measures[HOM] = lk[HOM][dindexr][0];
			else
				measures[HOM] += lk[HOM][dindexr][0];
		}

		if (msrbool[INV])
		{
			measures[INV] -= lk[INV][dindexl][0];
			if (abs(measures[INV]) <= ZERO_THRESH)
				measures[INV] = lk[INV][dindexr][0];
			else
				measures[INV] += lk[INV][dindexr][0];
		}

		if (msrbool[MU])
		{
			measures[MU] -= lk[MU][il+jl][0];
			if (abs(measures[MU]) <= ZERO_THRESH)
				measures[MU] = lk[MU][ir+jr][0];
			else
				measures[MU] += lk[MU][ir+jr][0];
		}

		if (msrbool[STD])
		{
			tstats[STDTEMP] -= (lk[STD][il][0] + lk[STD][jl][0]);
			if (abs(tstats[STDTEMP]) <= ZERO_THRESH)
				tstats[STDTEMP] = lk[STD][ir][0] + lk[STD][jr][0];
			else
				tstats[STDTEMP] += (lk[STD][ir][0] + lk[STD][jr][0]);
		}

		if (msrbool[COR])
		{
			tstats[CORTEMP] -= lk[COR][il][jl];
			if (abs(tstats[CORTEMP]) <= ZERO_THRESH)
				tstats[CORTEMP] = lk[COR][ir][jr];
			else
				tstats[CORTEMP] += lk[COR][ir][jr];
		}

		//ASM, ENT calculations: considering both (i,j) and (j,i) pairs at once
		if (msrbool[ASM] || msrbool[ENT])
		{
			// remove left column
			lmatold1 = glcm[il][jl];
			lmatnew1 = lmatold1 - 1;
			glcm[il][jl] = lmatnew1; //update matrix
			lmatold2 = glcm[jl][il];
			lmatnew2 = lmatold2 - 1;
			glcm[jl][il] = lmatnew2; //update matrix

			// add right column         
			rmatold1 = glcm[ir][jr];
			rmatnew1 = rmatold1 + 1;
			glcm[ir][jr] = rmatnew1; //update matrix
			rmatold2 = glcm[jr][ir];
			rmatnew2 = rmatold2 + 1;
			glcm[jr][ir] = rmatnew2; //update matrix        
		}

		if (msrbool[ASM])
		{
			measures[ASM] -= (lk[ASM][lmatold1][0] + lk[ASM][lmatold2][0]); //left pairs
			if (abs(measures[ASM]) <= ZERO_THRESH)
				measures[ASM] = lk[ASM][lmatnew1][0] + lk[ASM][lmatnew2][0];
			else
				measures[ASM] += (lk[ASM][lmatnew1][0] + lk[ASM][lmatnew2][0]);

			measures[ASM] -= (lk[ASM][rmatold1][0] + lk[ASM][rmatold2][0]); //right pairs
			if (abs(measures[ASM]) <= ZERO_THRESH)
				measures[ASM] = lk[ASM][rmatnew1][0] + lk[ASM][rmatnew2][0];
			else
				measures[ASM] += (lk[ASM][rmatnew1][0] + lk[ASM][rmatnew2][0]);
		}

		if (msrbool[ENT])  //NOTE: This result is x2 the original PCI code and is "correct"
		{
			//signs switched because ENT is -ve of summation

			measures[ENT] += (lk[ENT][lmatold1][0] + lk[ENT][lmatold2][0]); //left pairs
			if (abs(measures[ENT]) <= ZERO_THRESH)
				measures[ENT] = - lk[ENT][lmatnew1][0] - lk[ENT][lmatnew2][0];
			else
				measures[ENT] -= (lk[ENT][lmatnew1][0] + lk[ENT][lmatnew2][0]);

			measures[ENT] += (lk[ENT][rmatold1][0] + lk[ENT][rmatold2][0]); //right pairs
			if (abs(measures[ENT]) <= ZERO_THRESH)
				measures[ENT] = - lk[ENT][rmatnew1][0] - lk[ENT][rmatnew2][0];
			else
				measures[ENT] -= (lk[ENT][rmatnew1][0] + lk[ENT][rmatnew2][0]);
		}
	}

	//STD calculations
	if (msrbool[STD])
	{
		var = tstats[STDTEMP] - measures[MU]*measures[MU];
		if (var <= ZERO_THRESH)
		{
			var = 0;
			measures[STD] = 0;
		}
		else
			measures[STD] = sqrt(var);
	}

	//COR calculations
	if (msrbool[COR])
	{
		if (measures[STD] <= ZERO_THRESH)
			measures[COR] = 1;
		else
			measures[COR] = (tstats[CORTEMP] - measures[MU]*measures[MU]) / var;
	}
}



/***************** findrange : Find internal Start/Stop values *****************/
// Param in: nFlszX - window width
//			 nFlszY - window height
//			 nSpatX - delta x
//			 nSpatY - delta y
// Param out: Starting and ending values for X and Y directions
/*******************************************************************************/
void TTextureFeatExtract :: findrange(int* startX, int* stopX, int* startY, int* stopY,
	int nFlszX, int nFlszY, int nSpatX, int nSpatY)
{
	if (0 > -nSpatX) //start is max of the two
		*startX = 0;
	else
		*startX = -nSpatX;

	if (nFlszX > nFlszX-nSpatX) //end is min of the two
		*stopX = nFlszX-nSpatX;
	else
		*stopX = nFlszX;

	if (0 > -nSpatY) //start is max of the two
		*startY = 0;
	else
		*startY = -nSpatY;

	if (nFlszY > nFlszY-nSpatY) //end is min of the two
		*stopY = nFlszY-nSpatY;
	else
		*stopY = nFlszY;
}

/* --------------------------------------------------------------------------------------
The GaborFilteredImg provides the outputs of the Gabor filter bank
-----------------------------------------------------------------------------------------*/
void TTextureFeatExtract :: GaborFilteredImg(Matrix *FilteredImg_real, Matrix *FilteredImg_imag, Matrix *img, int side, double Ul, double Uh, int scale, int orientation, int flag)
{
	int h, w, xs, ys, border, r1, r2, r3, r4, hei, wid, s, n;
	Matrix *IMG, *IMG_imag, *Gr, *Gi, *Tmp_1, *Tmp_2, *F_1, *F_2, *G_real, *G_imag, *F_real, *F_imag;

	border = side;
	hei = img->height;
	wid = img->width;

	/* FFT2 */
	xs = (int) pow(2.0, ceil(logBase2((double)(img->height+2.0*border))));
	ys = (int) pow(2.0, ceil(logBase2((double)(img->width+2.0*border))));

	CreateMatrix(&IMG, xs, ys);

	r1 = img->width+border;
	r2 = img->width+2*border;
	for (h=0;h<border;h++) {
		for (w=0;w<border;w++)
			IMG->data[h][w] = img->data[border-1-h][border-1-w];
		for (w=border;w<r1;w++)
			IMG->data[h][w] = img->data[border-1-h][w-border];
		for (w=r1;w<r2;w++)
			IMG->data[h][w] = img->data[border-1-h][2*img->width-w+border-1];
	}

	r1 = img->height+border;
	r2 = img->width+border;
	r3 = img->width+2*border;
	for (h=border;h<r1;h++) {
		for (w=0;w<border;w++)
			IMG->data[h][w] = img->data[h-border][border-1-w];
		for (w=border;w<r2;w++)
			IMG->data[h][w] = img->data[h-border][w-border];
		for (w=r2;w<r3;w++)
			IMG->data[h][w] = img->data[h-border][2*img->width-w+border-1];
	}

	r1 = img->height+border;
	r2 = img->height+2*border;
	r3 = img->width+border;
	r4 = img->width+2*border;
	for (h=r1;h<r2;h++) {
		for (w=0;w<border;w++)
			IMG->data[h][w] = img->data[2*img->height-h+border-1][border-1-w];
		for (w=border;w<r3;w++)
			IMG->data[h][w] = img->data[2*img->height-h+border-1][w-border];
		for (w=r3;w<r4;w++)
			IMG->data[h][w] = img->data[2*img->height-h+border-1][2*img->width-w+border-1];
	}

	CreateMatrix(&F_real, xs, ys);
	CreateMatrix(&F_imag, xs, ys);
	CreateMatrix(&IMG_imag, xs, ys);

	Mat_FFT2(F_real, F_imag, IMG, IMG_imag);

	/* ----------- compute the Gabor filtered output ------------- */

	CreateMatrix(&Gr, 2*side+1, 2*side+1);
	CreateMatrix(&Gi, 2*side+1, 2*side+1);
	CreateMatrix(&Tmp_1, xs, ys);
	CreateMatrix(&Tmp_2, xs, ys);
	CreateMatrix(&F_1, xs, ys);
	CreateMatrix(&F_2, xs, ys);
	CreateMatrix(&G_real, xs, ys);
	CreateMatrix(&G_imag, xs, ys);

	for (s=0;s<scale;s++) {
		for (n=0;n<orientation;n++) {
			Gabor(Gr, Gi, s+1, n+1, Ul, Uh, scale, orientation, flag);
			Mat_Copy(F_1, Gr, 0, 0, 0, 0, 2*side, 2*side);
			Mat_Copy(F_2, Gi, 0, 0, 0, 0, 2*side, 2*side);
			Mat_FFT2(G_real, G_imag, F_1, F_2);

			Mat_Product(Tmp_1, G_real, F_real);
			Mat_Product(Tmp_2, G_imag, F_imag);
			Mat_Substract(IMG, Tmp_1, Tmp_2);

			Mat_Product(Tmp_1, G_real, F_imag);
			Mat_Product(Tmp_2, G_imag, F_real);
			Mat_Sum(IMG_imag, Tmp_1, Tmp_2);

			Mat_IFFT2(Tmp_1, Tmp_2, IMG, IMG_imag);

			Mat_Copy(FilteredImg_real, Tmp_1, s*hei, n*wid, 2*side, 2*side, hei+2*side-1, wid+2*side-1);
			Mat_Copy(FilteredImg_imag, Tmp_2, s*hei, n*wid, 2*side, 2*side, hei+2*side-1, wid+2*side-1);
		}
	}

	FreeMatrix(Gr);
	FreeMatrix(Gi);
	FreeMatrix(Tmp_1);
	FreeMatrix(Tmp_2);
	FreeMatrix(F_1);
	FreeMatrix(F_2);
	FreeMatrix(G_real);
	FreeMatrix(G_imag);
	FreeMatrix(F_real);
	FreeMatrix(F_imag);
	FreeMatrix(IMG);
	FreeMatrix(IMG_imag);
}


/* ------------------------------------------------------------------------------------------------------
The Gabor function generates a Gabor filter with the selected index 's' and 'n' (scale and orientation, 
respectively) from a Gabor filter bank. This filter bank is designed by giving the range of spatial 
frequency (Uh and Ul) and the total number of scales and orientations used to partition the spectrum. 

The returned filter is stored in 'Gr' (real part) and 'Gi' (imaginary part).
--------------------------------------------------------------------------------------------------------*/	
void TTextureFeatExtract :: Gabor(Matrix *Gr, Matrix *Gi, int s, int n, double Ul, double Uh, int scale, int orientation, int flag)
{
	double base, a, u0, z, Uvar, Vvar, Xvar, Yvar, X, Y, G, t1, t2, m;
	int x, y, side;

	base = Uh/Ul;
	a = pow(base, 1.0/(double)(scale-1));  // band_rad

	u0 = Uh/pow(a, (double) scale-s); // radius central frequency

	Uvar = (a-1.0)*u0/((a+1.0)*sqrt(2.0*log(2.0))); // 1/(2*Pi*sigmax)

	z = -2.0*log(2.0)*(Uvar*Uvar)/u0;
	Vvar = tan(PI/(2*orientation))*(u0+z)/sqrt(2.0*log(2.0)-z*z/(Uvar*Uvar));

	Xvar = 1.0/(2.0*PI*Uvar);
	Yvar = 1.0/(2.0*PI*Vvar);

	t1 = cos(PI/orientation*(n-1.0));
	t2 = sin(PI/orientation*(n-1.0));

	side = (int) (Gr->height-1)/2;

	for (x=0;x<2*side+1;x++) {
		for (y=0;y<2*side+1;y++) {
			X = (double) (x-side)*t1+ (double) (y-side)*t2;
			Y = (double) -(x-side)*t2+ (double) (y-side)*t1;
			G = 1.0/(2.0*PI*Xvar*Yvar)*pow(a, (double) scale-s)*exp(-0.5*((X*X)/(Xvar*Xvar)+(Y*Y)/(Yvar*Yvar)));
			Gr->data[x][y] = G*cos(2.0*PI*u0*X);
			Gi->data[x][y] = G*sin(2.0*PI*u0*X);
		}
	}

	/* if flag = 1, then remove the DC from the filter */
	if (flag == 1) {
		m = 0;
		for (x=0;x<2*side+1;x++)
			for (y=0;y<2*side+1;y++)
				m += Gr->data[x][y];

		m /= pow((double) 2.0*side+1, 2.0);

		for (x=0;x<2*side+1;x++)
			for (y=0;y<2*side+1;y++)
				Gr->data[x][y] -= m;
	}	
}

/* Function to generate spatial-frequency domain Gabor filter */
double TTextureFeatExtract :: GaborFreq(double u, double v, double sigma2, double freq, double lambda2, double coeff)
{
	double result = 0.0;

	result = exp(coeff*sigma2*( (u-freq)*(u-freq)*lambda2+v*v ) ); // coeff = -2*Pi^2

	return result;
}

/* Function to generate spatial-frequency domain POST Gabor filter */
double TTextureFeatExtract :: GaborPost(double u, double v, double sigma2, double lambda2, double coeff, double postlamda)
{
	double result = 0.0;

	result = exp(coeff*sigma2*( u*u*lambda2 + v*v )/postlamda); // coeff = -2*Pi^2

	return result;
}
