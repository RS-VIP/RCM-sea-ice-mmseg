/** \file
 *
 * \brief TGrayImage class
 *
 */
#pragma once 

#ifndef TGrayImage_Class
#define TGrayImage_Class

#include <cstring>
#include <assert.h>
#include "Image.h"
#include <iostream>
#include <fstream>
#include<math.h>

/// Base class for gray scale images.
template<typename T> 
class TGrayImage : public TImage
{

//
// Attributes
// 
protected:

	/// pointers to data array
	T** m_ppValues; 

	/// memory for all pixel values
	T* m_pMem;

// 
// Operations
// 
public: 


	/// Constructor
	TGrayImage() : TImage()
	{
		m_pMem = 0;
		m_ppValues = 0;
		m_nWidth = 0;
		m_nHeight = 0;
	};

	/// Constructor with memory allocation
	TGrayImage(int nWidth, int nHeight) : TImage()
	{
		m_pMem = 0;
		m_ppValues = 0;
		m_nWidth = 0;
		m_nHeight = 0;
		Alloc(nWidth, nHeight);
	};

	/// Copy constructor
	TGrayImage(TImage* pImg) : TImage()
	{
		m_pMem = 0;
		m_ppValues = 0;
		m_nWidth = 0;
		m_nHeight = 0;
		ConvertFrom(*pImg);
	};

	float getWidth(){
		return m_nWidth;
	};

	/// Destructor
	virtual ~TGrayImage()
	{
		Free();
	};

	virtual void WritePixelValuesToFile(const char *filename) const
	{
		std::ofstream file = std::ofstream(filename, ios::trunc);
		file.precision(7);
		file << std::scientific;
		
		for (int row = 0; row < m_nHeight; row++){
			for (int col = 0; col < m_nWidth; col++){
				file << GetPixelIntensity(row, col);
				if (col != m_nWidth - 1)
				{
					file << ',';
				}
			}
			file << '\n';
		}

		file.close();
	};

	virtual void ReplacePixelsOfValue(T valueToReplace, T valueToReplaceWith)
	{
		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
			{
				if (m_ppValues[i][j] == valueToReplace)
				{
				m_ppValues[i][j] = valueToReplaceWith;
				}
			}
		}
	}

	//Added by Michael to get value of given pixel
	T GetPixelData(int nCol, int nRow) const
	{
		if (nRow < 0 || nRow >= m_nHeight ||
			nCol < 0 || nCol >= m_nWidth)
			throw std::runtime_error("Array indices out of bounds!!");

		return m_ppValues[nRow][nCol];
	}

	//Michae Stone - May 2016
	virtual void ReplacePixel(int nCol, int nRow, T valueToReplaceWith)
	{
		m_ppValues[nCol][nRow] = valueToReplaceWith;
	}




private:

	/// Allocate memory
	void Alloc(int nWidth, int nHeight)
	{
		if (m_ppValues &&  // already allocated before
			(m_nWidth != nWidth || m_nHeight != nHeight))  // dimension inconsistent
			Free();

		if (!m_ppValues)
		{
			m_nWidth = nWidth;
			m_nHeight = nHeight;
			m_pMem = new T[nWidth * nHeight];
			m_ppValues = new T*[nHeight];
			for (int i = 0; i < nHeight; i ++)
				m_ppValues[i] = m_pMem + i * nWidth;
		}
	};

	/// Free memory
	void Free()
	{
		if (m_pMem)
		{
			delete[] m_pMem;  //QIN_MOD
			m_pMem = 0;
		}
		if (m_ppValues)
		{
			delete[] m_ppValues;  //QIN_MOD
			m_ppValues = 0;
		}
	};


public:

	virtual enImageType GetType() const { return TImage::GRAY; }

	/// Access the pixel value specified by row and column
	T& operator()(int row, int col) const
	{
		if (row < 0 || row >= m_nHeight ||
			col < 0 || col >= m_nWidth)
		{
			//Debug::WriteLine(row + ", " + m_nHeight + "\t" + col + ", " + m_nWidth);
            throw std::runtime_error("Array indices out of bounds!!");
		}
		return m_ppValues[row][col];
	};


	/// Get the pixel value specified by row and column, out of bounds handled by padding type.
	/**
	 * This does not let you change the pixel value. If paddingtype unknown, zero pads.
	 */
	T operator()(int row, int col, int paddingtype)
	{
		if (row < 0 || row >= m_nHeight ||
			col < 0 || col >= m_nWidth)
		{
			switch (paddingtype)
			{
				case PADMIRROR:
					row = MirrorMap(row, m_nHeight);
					col = MirrorMap(col, m_nWidth);
					break;
				case PADZERO:
					return 0;
				default:
					return 0; // zero pad.
			}

		}
		return m_ppValues[row][col];
	};



	/// When x is outside of [0, range), calculate the index that corresponds to mirror padding.
	/**
	 *  Otherwise, this preserves the value of x
	 */
	static int MirrorMap(int x, int Range)
	{
		while (x < 0 || x >= Range)
			if (x < 0) 
				x = -x - 1;
			else if (x >= Range)
				x = Range - 1 - (x - Range);
		return x;
	}

	/// Clip the value
	TGrayImage& Clip(T tMin, T tMax)
	{
		T tTemp;
		for (int i = 0; i < m_nHeight; i ++)
			for (int j = 0; j < m_nWidth; j ++)
		{
			if ((tTemp = m_ppValues[i][j]) < tMin)
				m_ppValues[i][j] = tMin;
			else if (tTemp > tMax)
				m_ppValues[i][j] = tMax;
		}
		return *this;
	};

 
	/// Sets all pixels to the specified value.
	void SetValue(T value)
	{
			for (int i = 0; i < m_nHeight; i ++)
				for (int j = 0; j < m_nWidth; j ++)
			{
				m_ppValues[i][j] = value;
			}
	}

	/// Get the pointer to data
	T** GetData() { return m_ppValues; };

// 
// Implementation
// 
public:

	/// Convert from the specified image
	virtual void ConvertFrom(const TImage &pSrcImg)
	{
		Alloc(pSrcImg.Width(), pSrcImg.Height());
		for (int i = 0; i < m_nHeight; i++)
			for (int j = 0; j < m_nWidth; j++)
			{
				m_ppValues[i][j] = (T)pSrcImg.GetPixelIntensity(i, j);
			}
	};

	/// Load a Tiff image	
	void loadTiff(unsigned char* tiffData,int width,int height, bool bBoundingBox, int startrow, int startcol, int stoprow, int stopcol)
	{     
		if (tiffData)
		{
			if (bBoundingBox)
			{
				if (startrow < 0)
					startrow = 0;
				if (startcol < 0)
					startcol = 0;

				if (stoprow > height-1)
					stoprow = height-1;
				if (stopcol > width-1)
					stopcol = width-1;
			}
			else
			{
				startrow = 0;
				startcol = 0;
				stoprow = height-1;
				stopcol = width-1;
			}

			Alloc(stopcol-startcol+1, stoprow-startrow+1);
			for (int i = startrow; i <= stoprow; i ++)
				for (int j = startcol; j <= stopcol; j ++)
				{
					m_ppValues[i-startrow][j-startcol] = (T) tiffData[i*width+j]; 
				}			
		}
	}

//    void loadGeoTiff(const char* fname, bool bVerbose = false)
//    {
//	    TIFF 	*tif=(TIFF*)0;  /* TIFF-level descriptor*/
//        GTIF	*gtif=(GTIF*)0; /* GeoKey-level descriptor */
//
//        tif=XTIFFOpen(fname,"r");
//        gtif = GTIFNew(tif);
//        if ((!gtif) || (!tif))
//        {
//            fprintf(stderr,"I failed to open your geotiff, sorry buddy.\n");
//            return;
//        }
//        int w=0, h=0;
//        int nSamplesPerPixel=0; // # of channels in the tiff image
//        int nBitsPerSample=0; // bit depth
//        int SampleFormat=0; // Sample format ((1) unsigned int, (2) 2's c. signed, (3) IEEE FP, (4) undefined )
//
//        // Get some info about the tiff
//        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
//        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
//        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &nBitsPerSample);
//        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nSamplesPerPixel);
//        TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &SampleFormat);
////        TIFFGetField(tif, TIFFTA);
//
//        if (SampleFormat==1 && nBitsPerSample==8){
//            // 8 bit integer samples
//            unsigned char* data = NULL;
//            data = (unsigned char*) _TIFFmalloc((tsize_t) (TIFFNumberOfStrips(tif) * TIFFStripSize(tif)));
//            for(tstrip_t strip=0; strip < TIFFNumberOfStrips(tif); strip++)
//                TIFFReadEncodedStrip(tif, strip, data + strip* TIFFStripSize(tif), (tsize_t)-1);
//            Alloc(w, h);
//            for (int i = 0; i <= h-1; i++)
//                for (int j = 0; j <= w-1; j++)
//                {
//                    m_ppValues[i][j] = (T)data[i*w+ j];
//                }
//
//        }
//        else if (SampleFormat==1 && nBitsPerSample==16){
//            // 16 bit integer samples
//            unsigned short* data = NULL;
//            data = (unsigned short*) _TIFFmalloc((tsize_t) (TIFFNumberOfStrips(tif) * TIFFStripSize(tif)));
//            for(tstrip_t strip=0; strip < TIFFNumberOfStrips(tif); strip++)
//                TIFFReadEncodedStrip(tif, strip, data + strip* TIFFStripSize(tif), (tsize_t)-1);
//            Alloc(w, h);
//            for (int i = 0; i <= h-1; i++)
//                for (int j = 0; j <= w-1; j++)
//                {
//                    m_ppValues[i][j] = (T)data[i*w+ j];
//                }
//
//        }
//        else
//        {
//            fprintf(stderr,"Unsupported tiff file format, right now I can only do single-channel tiffs of integer dtype.\n");
//            return;
//        }
//        return;
//    }

	void loadTiff(unsigned short* tiffData16bit, int width, int height, bool bBoundingBox, int startrow, int startcol, int stoprow, int stopcol)
	{
		if (tiffData16bit)
		{
			if (bBoundingBox)
			{
				if (startrow < 0)
					startrow = 0;
				if (startcol < 0)
					startcol = 0;

				if (stoprow > height - 1)
					stoprow = height - 1;
				if (stopcol > width - 1)
					stopcol = width - 1;
			}
			else
			{
				startrow = 0;
				startcol = 0;
				stoprow = height - 1;
				stopcol = width - 1;
			}

			Alloc(stopcol - startcol + 1, stoprow - startrow + 1);
			for (int i = startrow; i <= stoprow; i++)
				for (int j = startcol; j <= stopcol; j++)
				{
					m_ppValues[i - startrow][j - startcol] = (T)tiffData16bit[i*width + j];
				}
		}
	}

//    void loadTiff(unsigned char* tiffData8bit, int width, int height, bool bBoundingBox, int startrow, int startcol, int stoprow, int stopcol)
//    {
//        if (tiffData8bit)
//        {
//            if (bBoundingBox)
//            {
//                if (startrow < 0)
//                    startrow = 0;
//                if (startcol < 0)
//                    startcol = 0;
//
//                if (stoprow > height - 1)
//                    stoprow = height - 1;
//                if (stopcol > width - 1)
//                    stopcol = width - 1;
//            }
//            else
//            {
//                startrow = 0;
//                startcol = 0;
//                stoprow = height - 1;
//                stopcol = width - 1;
//            }
//
//            Alloc(stopcol - startcol + 1, stoprow - startrow + 1);
//            for (int i = startrow; i <= stoprow; i++)
//                for (int j = startcol; j <= stopcol; j++)
//                {
//                    m_ppValues[i - startrow][j - startcol] = (T)tiffData8bit[i*width + j];
//                }
//        }
//    }

    template <typename FromT>
    static TGrayImage<T>* ConvertFrom(TGrayImage<FromT>* cImg, int rowStart, int colStart, int iwidth, int iheight)
    {
        if (((colStart + iwidth - 1) > cImg->Width()) || ((rowStart + iheight -1) > cImg->Height()) || rowStart < 0 || colStart < 0)
	    {
		    printf("Dimension mismatch");
		    TGrayImage<T>* subImg =  new TGrayImage<T>(0,0);
		    return subImg;
	    } 

	    TGrayImage<T>* subImage = new TGrayImage<T>(iwidth,iheight);

	    T** subImgdata = subImage->GetData();
	    FromT** Imgdata = cImg->GetData();

	    for (int i = 0; i < (iheight); i++)
	    {
		    for ( int j = 0; j < (iwidth); j++)
		    {
			    subImgdata[i][j] = (T)(Imgdata[i + rowStart][j + colStart]);
		    }

	    }
	    return subImage;

    }

	virtual TImage* Clone() const
	{
		TGrayImage* pRet = new TGrayImage(m_nWidth, m_nHeight);
		memcpy(pRet->m_pMem, m_pMem, m_nWidth * m_nHeight * sizeof(T));
		return pRet;
	};

	/// Get the pixel intensity at the specified position
	virtual double GetPixelIntensity(int nRow, int nCol) const
	{
		return m_ppValues[nRow][nCol];
	};

	/// Write pixel value into a buffer
	virtual void WritePixelBytes(BYTE* pBuf, int nRow, int nCol)
	{
		*pBuf = (BYTE) m_ppValues[nRow][nCol];
	};

	/// Add two images by pixel to pixel.
	/**
	 * If the two images are of different forms, convert the parameter image to current form
	 */
	virtual TImage& operator+=(const TImage& img)
	{
		if (img.Width() != m_nWidth || img.Height() != m_nHeight)
            throw std::runtime_error("Array indices out of bounds!!");

		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
			{
				// The cast here is very dangerous.  We are assuming that the TImage is a TGrayImage.
				m_ppValues[i][j] = (T) (m_ppValues[i][j] + (*((TGrayImage*) &img))(i, j)); 
			}
		}

		return *this;
	};

    /// Sets the values of the image to be the pixel-by-pixel difference of the two given images
    /**
     *
     */
    virtual TGrayImage<T>& operator-=(const TGrayImage<T> & other)
    {
        if (other.Width() != m_nWidth || other.Height() != m_nHeight)
            throw std::runtime_error("Array indices out of bounds!!");
    
        for (int i = 0; i < m_nHeight; i ++) {
			for (int j = 0; j < m_nWidth; j ++)	{
				m_ppValues[i][j] = m_ppValues[i][j] - other.m_ppValues[i][j];
			}
		}
        return *this;
    }

    /// Compare two images pixel to pixel
    /**
     * Returns true if the images are exactly the same
     */
    virtual bool operator==(const TGrayImage<T> & other) 
    {
        if ( m_nWidth != other.Width() || m_nHeight != other.Height() )
            return false;

        for ( int i = 0; i < m_nHeight; ++i ) {
            for ( int j = 0; j < m_nWidth; ++j ) {
                if ( m_ppValues[i][j] != other.m_ppValues[i][j] ) { return false; }
            }
        }

        return true;
    }

    /// Compare two images pixel to pixel
    /**
     * Returns true if the images are different
     */
    virtual bool operator!=(const TGrayImage<T> & other) 
    {
        if ( m_nWidth != other.Width() || m_nHeight != other.Height() )
            return true;

        for ( int i = 0; i < m_nHeight; ++i ) {
            for ( int j = 0; j < m_nWidth; ++j ) {
                if ( m_ppValues[i][j] != other.m_ppValues[i][j] ) { return true; }
            }
        }

        return false;
    }
	
	/// Multiply two images by pixel to pixel.
	/**
	 * If the two images are of different forms, convert the parameter image to current form
	 */
	virtual TImage& operator*=(const TImage& img)
	{
		if (img.Width() != m_nWidth || img.Height() != m_nHeight)
            throw std::runtime_error("Array indices out of bounds!!");

		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
			{
				// The cast here is very dangerous.  We are assuming that the TImage is a TGrayImage.
				m_ppValues[i][j] = (T) (m_ppValues[i][j] * (*((TGrayImage*) &img))(i, j));
			}
		}

		return *this;
	};

	/// Add a value to current image.
	virtual TImage& operator+=(FLOAT val)
	{
		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
				m_ppValues[i][j] = (T) (m_ppValues[i][j] + val);
		}
		return *this;
	};

	/// Multiply all images pixels by an value
	virtual TImage& operator*=(FLOAT val)
	{
		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
				m_ppValues[i][j] = (T)(m_ppValues[i][j] * val);
		}
		return *this;
	};

	/// Get all absolute pixel values of the current image
	virtual TImage& Absolute()
	{
		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
			{
				if (m_ppValues[i][j] < 0)
					m_ppValues[i][j] = -m_ppValues[i][j];
				else
					m_ppValues[i][j] = m_ppValues[i][j];
			}
		}
		return *this;
	};

	/// Negate the current image
	virtual TImage& Negate()
	{
		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
				m_ppValues[i][j] = - m_ppValues[i][j];
		}
		return *this;
	};

	/// Get all inverse pixel values of the current image
	virtual TImage& Inverse()  //QIN_PROBLEM
	{
		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
				m_ppValues[i][j] = (T) (1.0 / m_ppValues[i][j]);
		}

		return *this;
	};

	/// Normalizes [minimum, maximum] values of image to [targetmin, targetmax].
	/**
	 * If usesourcerange = true, then normalizes [sourcemin, sourcemax] to [targetmin, targetmax] and clips all values outside of [sourcemin, sourcemax].
	 */
	virtual TImage& Normalize(T targetmin, T targetmax, bool usesourcerange = false, T sourcemin = 0, T sourcemax = 0)
	{
		if (usesourcerange)
		{
			Clip(sourcemin, sourcemax);
		}
		else
		{
			// find mininum / maximum.
			sourcemin = m_ppValues[0][0];
			sourcemax = m_ppValues[0][0];

			for (int i = 0; i < m_nHeight; i ++)
				for (int j = 0; j < m_nWidth; j ++)
			{
				if (m_ppValues[i][j] < sourcemin)
					sourcemin = m_ppValues[i][j];
				if (m_ppValues[i][j] > sourcemax)
					sourcemax = m_ppValues[i][j];
			}
		}

		// convert to double to avoid integer divides.
		DOUBLE sourcerange = (DOUBLE) sourcemax -  (DOUBLE) sourcemin;
		DOUBLE targetrange = (DOUBLE) targetmax -  (DOUBLE) targetmin;

		for (int i = 0; i < m_nHeight; i ++)
			for (int j = 0; j < m_nWidth; j ++)
		{
			m_ppValues[i][j] = (T) (((m_ppValues[i][j] - sourcemin) / sourcerange) * targetrange + targetmin);
			
		}
		return *this;
	}

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
	virtual TImage& Shift(int nDist, int nDir, int nType)
	{
		int nShiftLen;
		int r, c;
		T* pBuf;
		switch (nDir)
		{
		case SHIFTLEFT:
			if (nDist >= m_nWidth)
                throw std::runtime_error("Array indices out of bounds!!");
			nShiftLen = m_nWidth - nDist;
			pBuf = new T[nDist];
			for (r = 0; r < m_nHeight; r ++)
			{
				if (nType == PADZERO)
				{
					for (c = 0; c < nShiftLen; c ++)
						pBuf[c] = 0;
				}
				else if (nType == PADLASTLINE)
				{
					for (c = 0; c < nShiftLen; c ++)
						pBuf[c] = m_ppValues[r][m_nWidth - 1];
				}
				else if (nType == PADCYCLE)
					memcpy(pBuf, m_ppValues[r], nDist * sizeof(T));	
				memmove(m_ppValues[r], m_ppValues[r] + nDist, nShiftLen * sizeof(T));
				//memmove_s(m_ppValues[r], nShiftLen * sizeof(T), m_ppValues[r] + nDist, nShiftLen * sizeof(T));
				memcpy(m_ppValues[r] + nShiftLen, pBuf, nDist * sizeof(T)); 
			}
			break;

		case SHIFTRIGHT:
			if (nDist >= m_nWidth)
                throw std::runtime_error("Array indices out of bounds!!");
			nShiftLen = m_nWidth - nDist;
			pBuf = new T[nDist];
			for (r = 0; r < m_nHeight; r ++)
			{
				if (nType == PADZERO)
				{
					for (c = 0; c < nDist; c ++)
						pBuf[c] = 0;
				}
				else if (nType == PADLASTLINE)
				{
					for (c = 0; c < nDist; c ++)
						pBuf[c] = m_ppValues[r][0];
				}
				else if (nType == PADCYCLE)
					memcpy(pBuf, m_ppValues[r] + nShiftLen, nDist * sizeof(T));
				memmove(m_ppValues[r] + nDist, m_ppValues[r], nShiftLen * sizeof(T));
				//memmove_s(m_ppValues[r] + nDist, nShiftLen * sizeof(T), m_ppValues[r], nShiftLen * sizeof(T));
				memcpy(m_ppValues[r], pBuf, nDist * sizeof(T));
			}
			break;

		case SHIFTUP:
			if (nDist >= m_nHeight)
                throw std::runtime_error("Array indices out of bounds!!");
			nShiftLen = m_nHeight - nDist;
			pBuf = new T[m_nWidth * nDist];
			for (r = 0; r < nDist; r ++)
			{
				if (nType == PADZERO)
					memset(pBuf + r * m_nWidth, 0, m_nWidth * sizeof(T));
				else if (nType == PADLASTLINE)
					memcpy(pBuf + r * m_nWidth, m_ppValues[m_nHeight - 1], m_nWidth * sizeof(T));
				else if (nType == PADCYCLE)
					memcpy(pBuf + r * m_nWidth, m_ppValues[r], m_nWidth * sizeof(T));
			}
			//memmove_s(m_ppValues[0], nShiftLen * m_nWidth * sizeof(T), m_ppValues[nDist], nShiftLen * m_nWidth * sizeof(T));
			memmove(m_ppValues[0], m_ppValues[nDist], nShiftLen * m_nWidth * sizeof(T));
			memcpy(m_ppValues[nShiftLen], pBuf, nDist * m_nWidth * sizeof(T));
			break;

		case SHIFTDOWN:
			if (nDist >= m_nHeight)
                throw std::runtime_error("Array indices out of bounds!!");
			nShiftLen = m_nHeight - nDist;
			pBuf = new T[m_nWidth * nDist];
			for (r = 0; r < nDist; r ++)
			{
				if (nType == PADZERO)
					memset(pBuf + r * m_nWidth, 0, m_nWidth * sizeof(T));
				else if (nType == PADLASTLINE)
					memcpy(pBuf + r * m_nWidth, m_ppValues[0], m_nWidth * sizeof(T));
				else if (nType == PADCYCLE)
					memcpy(pBuf + r * m_nWidth, m_ppValues[nShiftLen + r], m_nWidth * sizeof(T));
			}
			//memmove_s(m_ppValues[nDist], nShiftLen * m_nWidth * sizeof(T), m_ppValues[0], nShiftLen * m_nWidth * sizeof(T));
			memmove(m_ppValues[nDist], m_ppValues[0], nShiftLen * m_nWidth * sizeof(T));
			memcpy(m_ppValues[0], pBuf, nDist * m_nWidth * sizeof(T));
			break;
		}

		delete pBuf;
		return *this;
	};

	/// Padding images
	/**
	 * Param in
	 * nPadOnTop: distance of padding for top
	 * nPadOnBottom: distance of padding for bottom
	 * nPadOnLeft: distance of padding for left
	 * nPadOnRight: distance of padding for right
	 * nPadType: padding type
	*/
	virtual TImage& Pad(int nPadOnTop, int nPadOnBottom, int nPadOnLeft, int nPadOnRight, int nPadType)
	{
		int nTotalWidth = m_nWidth + nPadOnLeft + nPadOnRight;
		int nTotalHeight = m_nHeight + nPadOnTop + nPadOnBottom;
		T* pDestMem = new T[nTotalWidth * nTotalHeight];
		T** ppDest = new T*[nTotalHeight];
		int i, j;
		for (i = 0; i < nTotalHeight; i ++)
			ppDest[i] = pDestMem + i * nTotalWidth;

		// Copy data
		for (i = 0; i < m_nHeight; i ++)
			memcpy(ppDest[i + nPadOnTop] + nPadOnLeft, m_ppValues[i], m_nWidth * sizeof(T));

		// Padding with zero
		if (nPadType == PADZERO)
		{
			for (i = 0; i < nPadOnTop; i ++)
			{
				for (j = 0; j < nTotalWidth; j ++)
					ppDest[i][j] = 0;
			}
			for (; i < nPadOnTop + m_nHeight; i ++)
			{
				for (j = 0; j < nPadOnLeft; j ++)
					ppDest[i][j] = 0;
				for (j = nPadOnLeft + m_nWidth; j < nTotalWidth; j ++)
					ppDest[i][j] = 0;
			}
			for (; i < nTotalHeight; i ++)
			{
				for (int j = 0; j < nTotalWidth; j ++)
					ppDest[i][j] = 0;
			}
		}
		// Padding with cycle values.
		else if (nPadType == PADCYCLE)
		{
			// Padding horizontally
			int nPadOnLeftTmp = nPadOnLeft;
			int nPadOnRightTmp = nPadOnRight;
			for (i = nPadOnTop; i < nPadOnTop + m_nHeight; i ++)
			{
				while ( nPadOnLeftTmp - m_nWidth > 0 )
				{
					for ( j = nPadOnLeftTmp - m_nWidth; j < nPadOnLeftTmp; j++ )
						ppDest[i][j] = m_ppValues[i - nPadOnTop][j - nPadOnLeftTmp + m_nWidth];
					nPadOnLeftTmp -= m_nWidth;
				}
				for (j = 0; j < nPadOnLeftTmp; i++)
					ppDest[i][j] = ppDest[i][j + m_nWidth];

				while ( nPadOnRightTmp - m_nWidth > 0)
				{
					for (j = nTotalWidth - nPadOnRightTmp; j < nTotalWidth - nPadOnRightTmp + m_nWidth; j ++)
						ppDest[i][j] = m_ppValues[i - nPadOnTop][j - nTotalWidth + nPadOnRightTmp];
					nPadOnRightTmp -= m_nWidth;
				}

				for (j = nTotalWidth - nPadOnRightTmp; j < nTotalWidth; j ++)
					ppDest[i][j] = ppDest[i][j - m_nWidth];
			}
			// Padding vertically
			int nPadOnTopTmp = nPadOnTop;
			int nPadOnBottomTmp = nPadOnBottom;
			while ( nPadOnTopTmp - m_nHeight > 0 )
			{
				for ( i = nPadOnTopTmp - m_nHeight; i < nPadOnTopTmp; i++ )
					memcpy(ppDest[i], ppDest[i + m_nHeight], nTotalWidth * sizeof(T));
				nPadOnTopTmp -= m_nHeight;
			}
			for ( i = 0; i < nPadOnTopTmp; i++ )
				memcpy(ppDest[i], ppDest[i + m_nHeight], nTotalWidth * sizeof(T));

			while ( nPadOnBottomTmp - m_nHeight > 0 )
			{
				for ( i = nTotalHeight - nPadOnBottomTmp; i < nTotalHeight - nPadOnBottomTmp + m_nHeight; i++ )
					memcpy(ppDest[i], ppDest[i - m_nHeight], nTotalWidth * sizeof(T));
				nPadOnTopTmp -= m_nHeight;
			}
			for ( i = nTotalHeight - nPadOnBottomTmp; i < nTotalHeight; i++ )
				memcpy(ppDest[i], ppDest[i - m_nHeight], nTotalWidth * sizeof(T));
		}
		// Padding with mirror values.
		else if (nPadType == PADMIRROR)
		{
			// Padding horizontally
			int nPadOnLeftTmp = nPadOnLeft;
			int nPadOnRightTmp = nPadOnRight;
			for (i = nPadOnTop; i < nPadOnTop + m_nHeight; i ++)
			{
				while ( nPadOnLeftTmp - m_nWidth > 0)
				{
					for ( j = nPadOnLeftTmp - m_nWidth; j < nPadOnLeftTmp; j++ )
						ppDest[i][j] = m_ppValues[i - nPadOnTop][nPadOnLeftTmp - 1 - j];
					nPadOnLeftTmp -= m_nWidth;
				}
				for (j = 0; j < nPadOnLeftTmp; j ++)
					ppDest[i][j] = ppDest[i][2*nPadOnLeftTmp - 1 - j];

				while ( nPadOnRightTmp - m_nWidth > 0)
				{
					for ( j = nTotalWidth - nPadOnRightTmp; j < nTotalWidth - nPadOnRightTmp + m_nWidth; j++ )
						ppDest[i][j] = m_ppValues[i - nPadOnTop][m_nWidth - 1 + nTotalWidth - nPadOnRightTmp - j];
					nPadOnRightTmp -= m_nWidth;
				}
				for (j = nTotalWidth - nPadOnRightTmp; j < nTotalWidth; j ++)
					ppDest[i][j] = ppDest[i][2*(nTotalWidth - nPadOnRightTmp) - 1 - j];
			}
			// Padding vertically
		    int nPadOnTopTmp = nPadOnTop;
			int nPadOnBottomTmp = nPadOnBottom;
			while (nPadOnTopTmp - m_nHeight > 0)
			{
				for (i = nPadOnTopTmp - m_nHeight; i < nPadOnTopTmp; i ++)
					memcpy(ppDest[i], ppDest[nPadOnTopTmp - 1 - i], nTotalWidth * sizeof(T));
				nPadOnTopTmp -= m_nHeight;
			}
			for (i = 0; i < nPadOnTopTmp; i ++)
				memcpy(ppDest[i], ppDest[2*nPadOnTopTmp - 1 - i], nTotalWidth * sizeof(T));

			while (nPadOnBottomTmp - m_nHeight > 0)
			{
				for (i = nTotalHeight - nPadOnBottomTmp; i < nTotalHeight - nPadOnBottomTmp + m_nHeight; i ++)
					memcpy(ppDest[i], ppDest[m_nHeight + nTotalHeight - nPadOnBottomTmp - 1 - i], nTotalWidth * sizeof(T));
				nPadOnBottomTmp -= m_nHeight;
			}
			for (i = nTotalHeight - nPadOnBottomTmp; i < nTotalHeight; i ++)
				memcpy(ppDest[i], ppDest[2*(nTotalHeight - nPadOnBottomTmp) - 1 - i], nTotalWidth * sizeof(T));
		}

		delete[] m_pMem;
		delete[] m_ppValues;
		m_pMem = pDestMem;
		m_ppValues = ppDest;
		m_nWidth = nTotalWidth;
		m_nHeight = nTotalHeight;
		return *this;
	};

	/// Depadding images
	/** 
	 * Param in
	 * nPadOnTop: distance of padding for top
	 * nPadOnBottom: distance of padding for bottom
	 * nPadOnLeft: distance of padding for left
	 * nPadOnRight: distance of padding for right
	 */
	virtual TImage& Depad(int nPadOnTop, int nPadOnBottom, int nPadOnLeft, int nPadOnRight)
	{
		int nTotalWidth = m_nWidth - nPadOnLeft - nPadOnRight;
		int nTotalHeight = m_nHeight - nPadOnTop - nPadOnBottom;
		if (nTotalWidth <= 0 || nTotalHeight <= 0)
			return *this;  // do nothing

		T* pDestMem = new T[nTotalWidth * nTotalHeight];
		T** ppDest = new T*[nTotalHeight];
		int i;
		for (i = 0; i < nTotalHeight; i ++)
			ppDest[i] = pDestMem + i * nTotalWidth;

		for (i = 0; i < nTotalHeight; i ++)
			memcpy(ppDest[i], m_ppValues[i + nPadOnTop] + nPadOnLeft, nTotalWidth * sizeof(T));

		delete [] m_pMem;
		delete [] m_ppValues;
		m_pMem = pDestMem;
		m_ppValues = ppDest;
		m_nWidth = nTotalWidth;
		m_nHeight = nTotalHeight;
		return *this;
	};

	/// Resize the image by resampling and interpolation
	virtual TImage& Resize(int nWidth, int nHeight)
	{
		T* pDestMem = new T[nWidth * nHeight];
		T** ppDest = new T*[nHeight];
		int i;
		for (i = 0; i < nHeight; i ++)
			ppDest[i] = pDestMem + i * nWidth;

		double dbWidthRatio = ((double) m_nWidth - 1) / nWidth;
		double dbHeightRatio = ((double) m_nHeight - 1) / nHeight;
		int r = 0, c = 0;
		double dbRowOff = 0, dbColOff = 0;
		for (i = 0; i < nHeight; i ++)
		{
			for (int j = 0; j < nWidth; j ++)
			{
				ppDest[i][j] = (T)((1 - dbRowOff) * (1- dbColOff) * m_ppValues[r][c] +
								   (1 - dbRowOff) * dbColOff * m_ppValues[r][c + 1] +
									dbRowOff * (1 - dbColOff) * m_ppValues[r + 1][c] +
									dbRowOff * dbColOff * m_ppValues[r + 1][c + 1]);
				dbColOff += dbWidthRatio;
				while (dbColOff > 1)
				{
					dbColOff -= 1.0;
					c ++;
				}
			}

			dbRowOff += dbHeightRatio;
			dbColOff = 0;
			c = 0;
			while (dbRowOff > 1)
			{
				dbRowOff -= 1.0;
				r ++;
			}
		}

		delete [] m_pMem;
		delete [] m_ppValues;
		m_pMem = pDestMem;
		m_ppValues = ppDest;
		m_nWidth = nWidth;
		m_nHeight = nHeight;
		return *this;
	};

    /// Nearest neighbor interpolation
    TImage& UpsampleNearestNeighbor(int nWidth, int nHeight)
    {
        // don't do anything if the sizes are bad.
        if (nWidth < m_nWidth || nHeight < m_nHeight)
            return *this;

        T* pDestMem = new T[nWidth * nHeight];
        T** ppDest = new T * [nHeight];
        int i, j;
        for (i = 0; i < nHeight; i++)
            ppDest[i] = pDestMem + i * nWidth;

        float xScale = m_nWidth / (float)nWidth;
        float yScale = m_nHeight / (float)nHeight;

        for (i = 0; i < nHeight; i++)
        {
            for (j = 0; j < nWidth; j++)
            {
                int px = (int) floor(j * xScale);
                int py = (int) floor(i * yScale);
                ppDest[i][j] = m_ppValues[py][px];
            }
        }

        delete[] m_pMem;
        delete[] m_ppValues;
        m_pMem = pDestMem;
        m_ppValues = ppDest;
        m_nWidth = nWidth;
        m_nHeight = nHeight;
        return *this;
    }

    /// Upsample by a factor of 2 with EPX ("Eric's Pixel Expansion")
    // see https://en.wikipedia.org/wiki/Pixel-art_scaling_algorithms#EPX/Scale2?AdvMAME2?
    TImage& Upsample2EPX()
    {
        int nWidth = 2 * m_nWidth;
        int nHeight = 2 * m_nHeight;

        T* pDestMem = new T[nWidth * nHeight];
        T** ppDest = new T * [nHeight];
        int i, j;
        for (i = 0; i < nHeight; i++)
            ppDest[i] = pDestMem + i * nWidth;

        for (i = 0; i < m_nHeight; i++)
        {
            for (j = 0; j < m_nWidth; j++)
            {
                ppDest[2 * i][2 * j] = m_ppValues[i][j];
                ppDest[2 * i + 1][2 * j] = m_ppValues[i][j];
                ppDest[2 * i][2 * j + 1] = m_ppValues[i][j];
                ppDest[2 * i + 1][2 * j + 1] = m_ppValues[i][j];

                if ((i > 0) && (i > 0) && (i < m_nHeight-1) && (j < m_nWidth - 1)) {
                    if ((m_ppValues[i][j-1] == m_ppValues[i+1][j]) && (m_ppValues[i][j-1] != m_ppValues[i-1][j]) && (m_ppValues[i+1][j] != m_ppValues[i][j+1]))
                        ppDest[2 * i + 1][2 * j] = m_ppValues[i+1][j];
                    if ((m_ppValues[i+1][j] == m_ppValues[i][j+1]) && (m_ppValues[i+1][j] != m_ppValues[i][j-1]) && (m_ppValues[i][j+1] != m_ppValues[i-1][j]))
                        ppDest[2 * i + 1][2 * j + 1] = m_ppValues[i][j+1];
                    if ((m_ppValues[i-1][j] == m_ppValues[i][j-1]) && (m_ppValues[i-1][j] != m_ppValues[i][j+1]) && (m_ppValues[i][j-1] != m_ppValues[i+1][j]))
                        ppDest[2 * i][2 * j] = m_ppValues[i][j-1];
                    if ((m_ppValues[i][j+1] == m_ppValues[i-1][j]) && (m_ppValues[i][j+1] != m_ppValues[i+1][j]) && (m_ppValues[i-1][j] != m_ppValues[i][j-1]))
                        ppDest[2 * i][2 * j + 1] = m_ppValues[i-1][j];
                }

            }
        }

        delete[] m_pMem;
        delete[] m_ppValues;
        m_pMem = pDestMem;
        m_ppValues = ppDest;
        m_nWidth = nWidth;
        m_nHeight = nHeight;
        return *this;
    }

	/// Block Average
	TImage& BlockAverage(int blocksize)
	{
		if (blocksize < 2)
			return *this;

		// new height and width. Truncates to nearest integer.
		int nWidth = m_nWidth / blocksize;
		int nHeight = m_nHeight / blocksize;

		// don't do anything if the sizes are bad.
		if (nWidth < 1 || nHeight < 1)
			return *this;

		T* pDestMem = new T[nWidth * nHeight];
		T** ppDest = new T*[nHeight];
		int i;
		for (i = 0; i < nHeight; i ++)
			ppDest[i] = pDestMem + i * nWidth;


		double blockpixels = blocksize * blocksize; // number of pixels in block.
		double averagedvalue; // average value of pixels.
		int j, r, c; // loop indices.
		// block average the values, with origin at top left corner
		// i.e. 2x2 block average populates upper left pixel with 
		// average of the 2x2 block at upper left corner of original.
		for (i = 0; i < nHeight; i ++)
		{
			for (j = 0; j < nWidth; j ++)
			{
				averagedvalue = 0.0;
	
				// find the average value for the block.
				// edge pixels without sufficient pixels for block size
				// are simply averaged mirror padding.
				for (r = i*blocksize; r < i*blocksize + blocksize; r++)
				{
					for (c = j*blocksize; c < j*blocksize + blocksize; c++)
					{
						averagedvalue += (*this)(r, c, PADMIRROR);
					}
				}

				ppDest[i][j] = (T) (averagedvalue / blockpixels);
			
			}
		}

		delete [] m_pMem;
		delete [] m_ppValues;
		m_pMem = pDestMem;
		m_ppValues = ppDest;
		m_nWidth = nWidth;
		m_nHeight = nHeight;
		return *this;
	}

	/// Get value range of the image intensity for subarea
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
	                               int nStartRow, int nStartCol, int nWidth, int nHeight)
	{
		int nEndRow = nStartRow + nHeight - 1;
		int nEndCol = nStartCol + nWidth - 1;
		if ((nStartRow >= m_nHeight || nStartRow < 0) ||
			(nStartCol >= m_nWidth || nStartCol < 0) ||
			(nEndRow >= m_nHeight || nEndRow < 0) ||
			(nEndCol >= m_nWidth || nEndCol < 0))
			return false;

		TGrayImage<BYTE>* pMaskFlag = 0;
		if (pMask)
			pMaskFlag = new TGrayImage<BYTE>(pMask);

		T val;
		T tMaxValue = m_ppValues[nStartRow][nStartCol];
		T tMinValue = tMaxValue;
		for (int i = nStartRow; i <= nEndRow; i ++)
		{
			for (int j = nStartCol; j <= nEndCol; j ++)
			{
				if (pMaskFlag && !(*pMaskFlag)(i, j))
					continue;

				if ((val = m_ppValues[i][j]) > tMaxValue)
					tMaxValue = val;
				if (val < tMinValue)
					tMinValue = val;
			}
		}

		dbMax = tMaxValue;
		dbMin = tMinValue;
		if (pMaskFlag)
			delete pMaskFlag;
		return true;
	};

	
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
	virtual int* GetIntensityHisto(double dbMin, double dbMax, int nHistCnt, TImage* pMask)
	{
		TGrayImage<BYTE>* pMaskFlag = 0;
		if (pMask)
			pMaskFlag = new TGrayImage<BYTE>(pMask);

		int* histo = new int[nHistCnt];
		int i;
		for (i = 0; i < nHistCnt; i ++)
			histo[i] = 0;

		double interval = (dbMax - dbMin) / nHistCnt;
		int index;

		for (i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
			{
				if (pMaskFlag && !(*pMaskFlag)(i, j))
					continue;

				index = (int)((m_ppValues[i][j] - dbMin) / interval);
				if (index < 0)
					index = 0;
				if (index >= nHistCnt)
					index = nHistCnt - 1;
				histo[index] ++;
			}
		}

		if (pMaskFlag)
			delete pMaskFlag;

		return histo;
	};

	//Shuhrat otsu
	virtual int* GetRawIntensityHisto(int nHistCnt, TImage* pMask)
	{
		TGrayImage<BYTE>* pMaskFlag = 0;
		if (pMask)
			pMaskFlag = new TGrayImage<BYTE>(pMask);

		int* histo = new int[nHistCnt];
		int i;
		for (i = 0; i < nHistCnt; i ++)
			histo[i] = 0;

		int index;

		for (i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
			{
				if (pMaskFlag && !(*pMaskFlag)(i, j))
					continue;

				index = (int)m_ppValues[i][j];
				histo[index] ++;
			}
		}

		if (pMaskFlag)
			delete pMaskFlag;

		return histo;
	};


	/// Compute the mean value (DC) of the image
	virtual FLOAT GetIntensityMean(TImage* pMask)
	{
		TGrayImage<BYTE>* pMaskFlag = 0;
		if (pMask)
			pMaskFlag = new TGrayImage<BYTE>(pMask);

		double dbMean = 0;
		int nTotalCnt = 0;
		for (int i = 0; i < m_nHeight; i ++)
			for (int j = 0; j < m_nWidth; j ++)
		{
			if (pMaskFlag && !(*pMaskFlag)(i, j))
				continue;

			dbMean += m_ppValues[i][j];
			nTotalCnt ++;
		}
		dbMean = dbMean / nTotalCnt;

		if (pMaskFlag)
			delete pMaskFlag;
		return (FLOAT) dbMean;
	};

	/// Dump the gray scale image pixel values
	/**
	 * Param out
	 * os: output stream where to dump to
	 */
	virtual void dump(ostream& os)
	{
		for (int i = 0; i < m_nHeight; i ++)
		{
			for (int j = 0; j < m_nWidth; j ++)
				os << (*this)(i, j) << " ";

			os << endl;
		}
	};

	void disconnect_regions(){
		int m_width = m_nWidth;
		int m_height = m_nHeight; 
		TGrayImage<int> temp_img = TGrayImage<int>(m_width,m_height);
		temp_img.SetValue(-1);
		int** temp_img_data = temp_img.GetData();
		int** img_data = this->GetData();

		int temp_regions_filled=0;

		//find the first 0 entry in temp_img and start filling from there
		int row,col;
		bool found_new_region;

		while(true){
			found_new_region=false;
			for(row=0;row<m_height && !found_new_region;row++)
				for(col=0;col<m_width && !found_new_region;col++){
					if(temp_img_data[row][col]==-1)
						found_new_region=true;
				}

				row--;
				col--;

				if(!found_new_region)  //true if all regions have been processed.
					break;  

				int img_fill_val = img_data[row][col];

				//Perform stack based flood fill algorithm.  See wikipedia flood fill.
				stack<int> rows;
				stack<int> cols;
				rows.push(row);
				cols.push(col);
				temp_img_data[row][col]=temp_regions_filled;
				while(!rows.empty()){
					row=rows.top();
					rows.pop();
					col=cols.top();
					cols.pop();

					//4 connected
					if(row+1<m_height && temp_img_data[row+1][col]==-1 && img_data[row+1][col] == img_fill_val){  //down
						temp_img_data[row+1][col]=temp_regions_filled;
						rows.push(row+1);
						cols.push(col);
					}
					if(row-1>=0 && temp_img_data[row-1][col]==-1 && img_data[row-1][col] == img_fill_val){  //up
						temp_img_data[row-1][col]=temp_regions_filled;
						rows.push(row-1);
						cols.push(col);
					}
					if(col+1<m_width && temp_img_data[row][col+1]==-1 && img_data[row][col+1] == img_fill_val){  //right
						temp_img_data[row][col+1]=temp_regions_filled;
						rows.push(row);
						cols.push(col+1);
					}
					if(col-1>=0 && temp_img_data[row][col-1]==-1 && img_data[row][col-1] == img_fill_val){  //left
						temp_img_data[row][col-1]=temp_regions_filled;
						rows.push(row);
						cols.push(col-1);
					}
				}
				temp_regions_filled++;
		}

		//copy temp image to original image
		for(row=0;row<m_height && !found_new_region;row++)
			for(col=0;col<m_width && !found_new_region;col++){
				img_data[row][col]=temp_img_data[row][col];
			}
	}

};

#endif

