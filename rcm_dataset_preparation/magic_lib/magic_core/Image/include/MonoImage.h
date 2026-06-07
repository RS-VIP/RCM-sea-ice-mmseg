/** \file
 *
 * \brief TMonoImage class.
 *
 */

#ifndef TMonoImage_Class
#define TMonoImage_Class


#include "GrayImage.h"

///Image of only 2 intensity levels
class TMonoImage : public TGrayImage<BYTE>
{
public:
	TMonoImage();
	TMonoImage(TImage* pImg);
	TMonoImage(int nWidth, int nHeight);
	TMonoImage(TImage* pImg, int threshold);

public:
	enum enBinaryValue {LOW = 0, HIGH = 255};

// Implementation
public:
	virtual enImageType GetType() const { return TImage::BINARY; }
	virtual ~TMonoImage();

};

#endif
