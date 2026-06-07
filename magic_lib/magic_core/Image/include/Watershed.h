/** \file
 *
 * \brief TWatershed class
 *
 */

#ifndef TWaterShed_Class
#define TWaterShed_Class


#include "DataType.h"
#include "GrayImage.h"
#include "MonoImage.h"
#include "Filter.h"

///Class of edge detection using watershed method. 
class TWaterShed
{
public:
	TWaterShed(int nMode);

// Attributes
public:
	// Constants
	
	/// neighbour mode to use for watershed algorithm.
	enum neighbor_mode {NEIGHBOR4 = 0,  // 4-neighbor system
	                    NEIGHBOR8 = 1};  // 8-neighbor system

	/// values for several cases of pixel / site types.
	enum map_constant {UNSUBMERGED = -1,  // value corresponding to unsubmerged sites
	                   MAPBORDER = -2,  // value indicating the border of regions in the map
	                   UNDETERMINED = -3,  // undetermined sites
	                   POSSIBLEMINIMA = -4}; // possible minima for seeds

protected:
	/// Neighborhood system to use.
	int m_nNeighborCnt;
	POINTex m_Neighbors[8];  // maximum 8 neighbors

	/// Pointer to the begining elements of corresponding water heights
	int* m_pRowPtr[257];  // here make the dimension 1 larger
	/// Pointer to the beginning elements of corresponding water heights
	int* m_pColPtr[257];  // here make the dimension 1 larger
	/// Array storing the positions of the pixels in ascending order of water height
	int* m_pRowData;
	int* m_pColData;

	/// Number of segments
	int m_nSegmentCnt;

// Operations
public:

	TGrayImage<int>* Segment(TGrayImage<int>* pGradientImg, TMonoImage* pMask,
	                         TFilter* pSmoothFilter);
	virtual void Segment(TGrayImage<int>*& pMapImg, TGrayImage<int>* pGradientImg,
	                     TMonoImage* pMask, int nBeginHeight);

	int GetNumSegments();


protected:

	TGrayImage<int>* GetInitialMap(TGrayImage<FLOAT>* GradientImg, TMonoImage* pMask, int numColBlks, int numRowBlks, bool bGlocal);

	void Sort(TGrayImage<int>* pGradientImg, TMonoImage* pMask);

	void ExtendRegion(int** Map, TGrayImage<int>* pGradientImg, TMonoImage* pMask,
	                  int nWaterHeight, int nRow, int nCol,
	                  TGrayImage<int>* pPrevPath, TGrayImage<int>* pCurrPath);
	BOOL IsBorder(int** Map, int Row, int Col, int Width, int Height);

// Implementation
public:
	virtual ~TWaterShed();

};

#endif

/////////////////////////////////////////////////////////////////////////////
