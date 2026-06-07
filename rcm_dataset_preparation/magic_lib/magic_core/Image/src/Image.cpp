// Image.cpp : implementation of the TImage class
//
// Base class for miscellaneous images of different intensity or color levels,
// It's the object for all image operation.
/////////////////////////////////////////////////////////////////////////////

#include "Image.h"
//#include "Bitmap.h"

//
// Constructor
//
TImage::TImage()
{
	m_nWidth = 0;
	m_nHeight = 0;
}

//
// Destructor
//
TImage::~TImage()
{
}

////
//// Load the image from file
////
//BOOL TImage::Load(const char* szFilePath)
//{
//    Mat img = imread(szFilePath, IMREAD_UNCHANGED);
////	// Read the input picture
////	TBitmap *pBmp = new TBitmap(szFilePath, ios::in);
////	TImage* pImg = pBmp->GetImageResource();
////	if (pImg == 0)
////	{
////		cout << "Input image error\n";
////		delete pBmp;
////		return false;
////	}
////
////	ConvertFrom(*pImg);
////	delete pBmp;
////	delete pImg;
////	return true;
//}

////
//// Save the image to file
////
//void TImage::Save(const char* szFilePath)
//{
//	// Read the input picture
//	TBitmap *pBmp = new TBitmap(szFilePath, ios::out);
//	pBmp->SaveImage(this);
//	delete pBmp;
//}

//void TImage::SaveTIFF(const char* szFilePath)
//{
//    TIFF* out = TIFFOpen(szFilePath, "w");
//    if(out){
//
//        int bitspersample = 16;
//        int nsamples=1;
//
//        TIFFSetField(out, TIFFTAG_IMAGELENGTH, m_nHeight);
//        TIFFSetField(out, TIFFTAG_IMAGEWIDTH,  m_nWidth);
//        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, nsamples);
//        TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
//        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bitspersample);
//        TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, m_nWidth));
//
//        uint16* buf = new uint16[m_nHeight*m_nWidth];
//
//        for(int i=0; i<m_nHeight; i++){
//            for(int j=0; j<m_nWidth; j++){
//                for(int n=0; n < bitspersample;) {
//                   ...nevemind I got bored i;ll just use the bitmap saving code for now
//                }
//            }
//
//        }
//    }
//}


//
// Public access of the image width
//
int TImage::Width() const
{
	return m_nWidth;
}

//
// Public access of the image height
//
int TImage::Height() const
{
	return m_nHeight;
}

//
// Retrieve a sub area of the image
//
// Param in
//  nRow, nCol: starting position
//  nWidth, nHeight: width and height of the area to be retrieved
//
TImage* TImage::operator()(int nRow, int nCol, int nWidth, int nHeight) const
{
	if (nRow < 0 || nRow + nHeight > m_nHeight ||
	    nCol < 0 || nCol + nWidth > m_nWidth)
		return 0;

	int nPadOnTop = nRow;
	int nPadOnBottom = m_nHeight - nRow - nHeight;
	int nPadOnLeft = nCol;
	int nPadOnRight = m_nWidth - nCol - nWidth;

	TImage* pResult = Clone();
	pResult->Depad(nPadOnTop, nPadOnBottom, nPadOnLeft, nPadOnRight);
	return pResult;
}

//
// Get value range of the image intensity
// 
// Param out
//  dbMin, dbMax: minimum and max value of image intensity
// 
void TImage::GetIntensityRange(double& dbMax, double& dbMin, TImage* pMask)
{
	GetIntensityRange(dbMax, dbMin, pMask, 0, 0, m_nWidth, m_nHeight);
}

