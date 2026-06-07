/** \file
 *
 * \brief TImage class
 *
 */

#ifndef TImage_Class
#define TImage_Class

#include "DataType.h"
#include <iostream>
//#include "geotiff/xtiffio.h"
//#include "geotiff/geotiff.h"
//#include "tiffio.h"

using namespace std;

#define SHIFTLEFT 0
#define SHIFTRIGHT 1
#define SHIFTUP 2
#define SHIFTDOWN 3

#define PADZERO 0
#define PADLASTLINE 1
#define PADCYCLE 2
#define PADMIRROR 3

/// Base class for images
class TImage
{
public:
	TImage();

// Attributes
public:
	enum enImageType {TRUECOLOR = 0, GRAY = 1, BINARY = 2};

protected:
	/// width and height of the image
	int m_nWidth, m_nHeight; 
	/// image type

// Operations
public:
	int Width() const;
	int Height() const;

	virtual enImageType GetType() const = 0;

	/// Load from specified file.
//	BOOL Load(const char* szFilePath);
	/// Save to specified file.
//	void Save(const char* szFilePath);

// Implementation
public:
	/// Convert from the specifed image
	virtual void ConvertFrom(const TImage& pSrcImg) = 0;
	/// Clone
	virtual TImage* Clone() const = 0;

	/// Retrieve a sub area of the image and return as TImage
	/**
	 * Param in
	 * nRow, nCol: starting position
	 * nWidth, nHeight: width and height of the area to be retrieved
	 */
	virtual TImage* operator()(int nRow, int nCol, int nWidth, int nHeight) const;

	/// Get the pixel intensity at the specified position
	virtual double GetPixelIntensity(int nRow, int nCol) const = 0;

	/// Write pixel value into a buffer
	virtual void WritePixelBytes(BYTE* pBuf, int nRow, int nCol) = 0;

	/// Add value to all pixels.
	virtual TImage& operator+=(float fVal) = 0;
	/// Multiply value to all pixels.
	virtual TImage& operator*=(float fVal) = 0;

	/// Get all absolute pixel values of the current image
	virtual TImage& Absolute() = 0;
	/// Negate the current image
	virtual TImage& Negate() = 0;
	/// Get all inverse pixel values of the current image
	virtual TImage& Inverse() = 0;

	/// Shift images according to specified direction and padding type
	/** 
	 * Param in
	 * nDist: distance of the shift
	 *
	 * nDir: direction of the shift
	 *
	 * nType: padding type
	 *
	 */
	virtual TImage& Shift(int nDistance, int nDirection, int nType) = 0;

	/// Padding images
	/**
	 * Param in
	 * nPadOnTop: distance of padding for top
	 * nPadOnBottom: distance of padding for bottom
	 * nPadOnLeft: distance of padding for left
	 * nPadOnRight: distance of padding for right
	 * nPadType: padding type
	*/	
	virtual TImage& Pad(int nPadOnTop, int nPadOnBottom, int nPadOnLeft, int nPadOnRight, int nPadType) = 0;

	/// Depadding images
	/** 
	 * Param in
	 * nPadOnTop: distance of padding for top
	 * nPadOnBottom: distance of padding for bottom
	 * nPadOnLeft: distance of padding for left
	 * nPadOnRight: distance of padding for right
	 */	
	virtual TImage& Depad(int nPadOnTop, int nPadOnBottom, int nPadOnLeft, int nPadOnRight) = 0;

	/// Resize the image
	virtual TImage& Resize(int nWidth, int nHeight) = 0;

	/// Block Average
	virtual TImage& BlockAverage(int blocksize) = 0;

	/// Get value range of the image intensity for entire image.
	/**
	 * Param in
	 *  pMask: mask image
	 * 
	 * Param out
	 *  dbMax, dbMin: max and minimum value of intensity
	 */	
	void GetIntensityRange(double& dbMax, double& dbMin, TImage* pMask);

	/// Get value range of the image intensity for sub-area.
	/**
	 * Param in
	 *  pMask: mask image
	 *  nStartRow, nStartCol: starting position of the area
	 *  nWidth, nHeight: width and height of the area
	 * 
	 * Param out
	 *  dbMax, dbMin: max and minimum value of intensity
	 */	
	virtual bool GetIntensityRange(double& dbMax, double& dbMin, TImage* pMask,
	                               int nStartRow, int nStartCol, int nWidth, int nHeight) = 0;

	///Get histogram of the image intensity
	/** 
	 * Param in
	 *  dbMin, dbMax: minimum and maximum value of the intensity
	 *  nHistCnt: number of histogram bins
	 *  pMask: mask
	 * 
	 * Return
	 *  the histogram in the form of an integer array
	 */	
	virtual int* GetIntensityHisto(double dbMin, double dbMax, int nHistCnt, TImage* pMask) = 0;

	/// Compute the mean value (DC) of the image	
	virtual float GetIntensityMean(TImage* pMask) = 0;

	/// Dump the gray scale image pixel values
	/**
	 * Param out
	 * os: output stream where to dump to
	 */	
	virtual void dump(ostream& os) {};

	virtual ~TImage();

};

#endif

