// WaterShed.cpp : implementation of the TWaterShed class
//
// Class of edge detection using watershed method
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Watershed.h"
#include <vector>

//
// Constructor
/////////////////////////////////////////////////////////////////////////////
TWaterShed::TWaterShed(int nMode)
{
	if (nMode == NEIGHBOR4)
	{
		m_nNeighborCnt = 4;
		memcpy(m_Neighbors, neighbor4, 4 * sizeof(POINTex));
	}
	else
	{
		m_nNeighborCnt = 8;
		memcpy(m_Neighbors, neighbor8, 8 * sizeof(POINTex));
	}

	m_pRowData = 0;
	m_pColData = 0;
}


//
// Destructor
/////////////////////////////////////////////////////////////////////////////
TWaterShed::~TWaterShed()
{
	if (m_pRowData)
		delete[] m_pRowData;  //QIN_MOD
	if (m_pColData)
		delete[] m_pColData;  //QIN_MOD
}


//
// Sort the gradient image by the value
// 
// Param in
//  pGradientImg: the gradient image
//  pMask: the mask
// 
void TWaterShed::Sort(TGrayImage<int>* pGradientImg, TMonoImage* pMask)
{
	// Count the occurrence of each value
	int nHist[256];
	memset(nHist, 0, 256 * sizeof(int));
	int i;
	int** imgdata;
	/*
	try
	{
		imgdata = pGradientImg->GetData();
	}
	catch (Exception^)
	{
		return;
	}*/
	imgdata = pGradientImg->GetData();
	for (i = pGradientImg->Height() - 1; i >= 0; i --)
		for (int j = pGradientImg->Width() - 1; j >= 0; j --)
		{
			if (pMask && !(*pMask)(i, j))
				continue;
			nHist[imgdata[i][j]] ++;
		}

		// Make the pointers point to the beginning element of each value
		m_pRowPtr[0] = m_pRowData;
		m_pColPtr[0] = m_pColData;
		for (i = 1; i <= 256; i ++)
		{
			m_pRowPtr[i] = m_pRowPtr[i - 1] + nHist[i - 1];
			m_pColPtr[i] = m_pColPtr[i - 1] + nHist[i - 1];
		}

		// TODO: The value k in the code below sometimes goes > 255.
		// TODO: This causes a crash. Fix it.

		// Fill in the elements
		memset(nHist, 0, 256 * sizeof(int));
		int k;
		for (i = pGradientImg->Height() - 1; i >= 0; i --)
			for (int j = pGradientImg->Width() - 1; j >= 0; j --)
			{
				if (pMask && !(*pMask)(i, j))
					continue;	
				k = imgdata[i][j];
				m_pRowPtr[k][nHist[k]] = i;
				m_pColPtr[k][nHist[k]] = j;
				nHist[k] ++;
			}
}



// Locate the edge using watershed algorithm
// 
// Param in
//  pGradientImg: the gradient image
//  pMask: the mask
//  fThreshold: initial threshold for determining segment seeds
// 
// Return
//  result map
// 
// TODO (rdmartin): 2011-07-11 Commented this function because it appears to be unused; can likely be deleted
//TGrayImage<int>* TWaterShed::Segment(TGrayImage<int>* pGradientImg, TMonoImage* pMask,
//                                     FLOAT fThreshold)
//{
//	// Get initial water segments
//	TGrayImage<int>* pMap = GetInitialMap(pGradientImg, pMask, fThreshold);
//	// Segmentation
//	Segment(pMap, pGradientImg, pMask, (int) (fThreshold + 1));
//	return pMap;
//}

// Locate the edge using watershed algorithm
// 
// Param in
//  pGradientImg: the gradient image
//  pMask: the mask
//  pSeeds: seeds
//  nSeedsLabel: seeds labels
//  nSeedsCnt: number of seeds
//
// Return
//  result map
//
// TODO (rdmartin): 2011-07-11 Commented this function because it appears to be unused; can likely be deleted
//TGrayImage<int>* TWaterShed::Segment(TGrayImage<int>* pGradientImg, TMonoImage* pMask,
//                                     POINTex* pSeeds, int* nSeedsLabel, int nSeedsCnt)
//{
//	// Get initial water segments
//	TGrayImage<int>* pMap = GetInitialMap(pGradientImg, pMask,
//	                                         pSeeds, nSeedsLabel, nSeedsCnt);
//	// Segmentation
//	Segment(pMap, pGradientImg, pMask, 0);
//	return pMap;
//}


//
// Locate the edge using watershed algorithm
// 
// Param in
//  pGradientImg: the gradient image
//  pMask: the mask
//  pFilter: the smoothing filter for locating local minima as seeds
// 
// Return
//  result map
// 

int TWaterShed::GetNumSegments()
{
	return m_nSegmentCnt;
}

TGrayImage<int>* TWaterShed::Segment(TGrayImage<int>* pGradientImg, TMonoImage* pMask,
	TFilter* pFilter)
{
	// Smooth the gradient image  // NOT suggested since gradient magnitude is normalized within mask and this operation is not related to mask!!
	TGrayImage<FLOAT>* pSmoothedImg;
	if (pFilter)
	{
		TGrayImage<FLOAT>* pImg = new TGrayImage<FLOAT>(pGradientImg);
		pSmoothedImg = pFilter->Conv2(pImg);
		delete pImg;
	}
	else
		pSmoothedImg = new TGrayImage<FLOAT>(pGradientImg);

	// Get initial water segments
	TGrayImage<int>* pMap = GetInitialMap(pSmoothedImg, pMask, 0, 0, false);
	// Segmentation
	Segment(pMap, pGradientImg, pMask, 0);

	delete pSmoothedImg;
	return pMap;
}


//
// Locate the edge using watershed algorithm
// 
// Param in
//  pMap: the initial segmentation map
//  pGradientImg: the gradient image
//  pMask: the mask
//  nBeginHeight: the begining water line height
// 
// Param out
//  pMap: the segmentation map
// 
void TWaterShed::Segment(TGrayImage<int>*& pMap, TGrayImage<int>* pGradientImg,
	TMonoImage* pMask, int nBeginHeight)
{
	int** map = pMap->GetData();

	// Initialize the struct
	int nWidth = pGradientImg->Width();
	int nHeight = pGradientImg->Height();
	if (m_pRowData)
		delete[] m_pRowData;
	if (m_pColData)
		delete[] m_pColData;
	m_pRowData = new int[nWidth * nHeight];
	m_pColData = new int[nWidth * nHeight];
	Sort(pGradientImg, pMask);

	TGrayImage<int>* pCurrPath = new TGrayImage<int>(nWidth, nHeight);
	TGrayImage<int>* pPrevPath = new TGrayImage<int>(nWidth, nHeight);

	// Water flooding
	int *pRow, *pCol;
	int r, c;
	for (int i = nBeginHeight; i <= 255; i ++)
	{
		pRow = m_pRowPtr[i];
		pCol = m_pColPtr[i];
		for (; pRow < m_pRowPtr[i + 1]; pRow = pRow + 1, pCol = pCol + 1)
		{
			if (map[*pRow][*pCol] != UNSUBMERGED ||  // already flooded
				(pMask && !(*pMask)(*pRow, *pCol)))
				continue;

			for (int k = m_nNeighborCnt - 1; k >= 0; k --)
			{
				r = *pRow + m_Neighbors[k].Row;
				c = *pCol + m_Neighbors[k].Col;
				if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
					continue;

				if (map[r][c] >= 0)
				{
					map[*pRow][*pCol] = map[r][c];
					break;
				}
			}

			if (map[*pRow][*pCol] >= 0)
			{
				// Check if current pixel should be a boundary pixel
				if (IsBorder(map, *pRow, *pCol, nWidth, nHeight))
					map[*pRow][*pCol] = MAPBORDER;
				else
					// Basin boundary point
					ExtendRegion(map, pGradientImg, pMask, i, *pRow, *pCol,
					pPrevPath, pCurrPath);
			}
		}
	}

	// Free
	delete pPrevPath;
	delete pCurrPath;

}


//
// Check if current flooding has reached the border
// 
// Param in
//  pMap: the current segmentation configuration obtained
//  nRow, nCol: current position

// 
// Return
//  true: border; false: not border
// 
BOOL TWaterShed::IsBorder(int** pMap, int nRow, int nCol, int nWidth, int nHeight)
{
	// Check if pixel on the other side is of another segment, the direction
	// can be deviated by 2 from the coming direction
	int nLabel = pMap[nRow][nCol], m;
	int r, c;
	int k;
	for (k = m_nNeighborCnt - 1; k >= 0; k --)
	{
		r = nRow + m_Neighbors[k].Row;
		c = nCol + m_Neighbors[k].Col;
		if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
			continue;

		m = pMap[r][c];
		if (m >= 0 && m != nLabel)
			break;
	}
	if (k >= 0)
		// belong to final shed line
		return true;
	else
		return false;
}

// Get the initial map (seeds region with others undetermined)
// 
// Param in
//  pGradientImg: the gradient image
//  pMask: the mask
// 
// Return
//  initial configuration map
// 
TGrayImage<int>* TWaterShed::GetInitialMap(TGrayImage<FLOAT>* pGradientImg, TMonoImage* pMask, int numColBlks, int numRowBlks, bool bGlocal)
{
	FLOAT** img = pGradientImg->GetData();
	int nWidth = pGradientImg->Width();
	int nHeight = pGradientImg->Height();

	// Pad the map image for convenience and initialize it
	TGrayImage<int>* pMap = new TGrayImage<int>(nWidth + 2, nHeight + 2);
	int** map = pMap->GetData();

	int i, j;

	//create a raster to store watershed minimas
	vector<vector<int>> minimas;
	for (i = 0; i < nHeight; i ++)
	{		 
		vector<int> row;
		for (j = 0; j < nWidth; j ++)
		{
			row.push_back(0);
		}
		minimas.push_back(row);
	}


	for (i = 0; i < nWidth + 2; i ++)
	{
		map[0][i] = MAPBORDER;  // image border
		map[nHeight + 1][i] = MAPBORDER;  // image border
	}
	for (i = 0; i < nHeight + 2; i ++)
	{
		map[i][0] = MAPBORDER;  // image border
		map[i][nWidth + 1] = MAPBORDER;  // image border
	}
	for (i = 1; i <= nHeight; i ++)
	{
		for (j = 1; j <= nWidth; j ++)
			map[i][j] = UNDETERMINED;  // undetermined
	}

	// Label local minima of gradient as seeds
	int nRow, nCol, r, c;
	int nLabel = 0;
	double dbLocalMinima;  // local minima of gradient
	for (i = 0; i < nHeight; i ++)
	{
		nRow = i + 1;
		for (j = 0; j < nWidth; j ++)
		{
			nCol = j + 1;
			if (pMask && !(*pMask)(i, j))
			{
				map[nRow][nCol]=MAPBORDER;
				continue;
			}
			// Find the neighbor that has minimum gradient
			dbLocalMinima = img[i][j];
			for (int k = m_nNeighborCnt - 1; k >= 0; k --)
			{
				r = nRow + m_Neighbors[k].Row;
				c = nCol + m_Neighbors[k].Col;
				if (map[r][c] == MAPBORDER)
					// Out of boundary
					continue;
				if (pMask && !(*pMask)(r - 1, c - 1))
					continue;
				if (img[r-1][c-1] < dbLocalMinima)
				{
					// Neighbor with lower gradient found, current cannot be minima
					map[nRow][nCol] = UNSUBMERGED;
					break;
				}
				else if(img[r-1][c-1] == dbLocalMinima)
				{
					if (map[r][c] == UNSUBMERGED)  
					{
						// Neighbor with same gradient is not minima, current
						// cannot be minima either
						map[nRow][nCol] = UNSUBMERGED;
						break;
					}
					else if (map[r][c] != POSSIBLEMINIMA)
						// Current may be a minima, however need to further check
						// those neighbors that have the same gradient
						map[nRow][nCol] = POSSIBLEMINIMA;
				}
			}
			if (map[nRow][nCol] == UNDETERMINED)  //Neighbors are either with higher gradient magnitudes or labeled as POSSIBLEMINIMA
			{
				// Local minima
				map[nRow][nCol] = nLabel ++;
				//Recording original minima
				minimas[nRow-1][nCol-1] = -1;
			}


		}
	}



	if(bGlocal)
	{
		bool found = false;
		int nRowBlockSize = nHeight/numRowBlks;
		int nColBlockSize = nWidth/numColBlks;
		int prevRegionCount = nLabel;
		nLabel = 0;
		int minrow, mincol, maxrow, maxcol;/*
		minrow = 1;
		mincol = 1;
		maxrow = nRowBlockSize;
		maxcol = nColBlockSize;*/

		for (i=0;i<numRowBlks;i++)
			for (j=0;j<numColBlks;j++)
			{
				minrow=i*nRowBlockSize;
				maxrow=(i+1)*nRowBlockSize-1;
				mincol=j*nColBlockSize;
				maxcol=(j+1)*nColBlockSize-1;
				if (i==numRowBlks-1)
					maxrow=nHeight-1;
				if (j==numColBlks-1)
					maxcol=nWidth-1;

				vector<FLOAT> pixelGrad;
				vector<int> pixelRealIndex;
				int realIndex=0;
				for(int row = minrow; row<=maxrow; row++)
				{
					for(int col = mincol; col<=maxcol; col++)
					{
						if(minimas[row][col] == -1)
						{
							pixelGrad.push_back(img[row][col]);
							pixelRealIndex.push_back(realIndex);
						}
						realIndex++;
					}
				}

				//FindMin
				if(pixelGrad.size()>0)
				{
					int index= 0;
					int col1;
					float min = pixelGrad[0];

					for (col1 = 0; col1<pixelGrad.size(); col1++)
					{
						if(min > pixelGrad[col1])
						{
							min = pixelGrad[col1];
							index = col1;
						}
					}

					int rowMin, colMin;

					rowMin = minrow+pixelRealIndex[index]/(maxcol-mincol+1);
					colMin = mincol+pixelRealIndex[index]%(maxcol-mincol+1);

					//mark all other pixels as unsubmerged in the window except the minima in that window
					for(int row = minrow; row<=maxrow; row++)
					{
						for(int col = mincol; col<=maxcol; col++)
						{
							if(minimas[row][col] == -1)
							{
								map[row+1][col+1] = UNSUBMERGED;
								minimas[row][col] = 0;
							}
						}										
					}

					map[rowMin+1][colMin+1] = nLabel ++;
					minimas[rowMin][colMin] = -1;
				}
			}
	}

	// Finalize the initial segmentation map
	m_nSegmentCnt = nLabel;
	pMap->Depad(1, 1, 1, 1);
	map = pMap->GetData();
	for (i = 0; i < nHeight; i ++)
	{
		for (j = 0; j < nWidth; j ++)
		{
			if (map[i][j] == POSSIBLEMINIMA)
				map[i][j] = UNSUBMERGED;
		}
	}

	return pMap;
}


//
// Extend the region where the height is lower than current water line height
// 
// Param in
//  map: the current segmentation map
//  pGradientImg: the gradient image
//  nRow, nCol: starting position
//  pPrevPath, pCurrPath: temporary space to preserve path. Pass in here 
//                        to save workload
// 
// Param out
//  map: the new segmentation map
/////////////////////////////////////////////////////////////////////////////
void TWaterShed::ExtendRegion(int** map, TGrayImage<int>* pGradientImg,
	TMonoImage* pMask, int nWaterHeight,
	int nRow, int nCol,
	TGrayImage<int>* pPrevPath, TGrayImage<int>* pCurrPath)
{
	int nWidth = pGradientImg->Width();
	int nHeight = pGradientImg->Height();
	int** img = pGradientImg->GetData();

	// The following labels all pixels connected to and below the water height 
	// of current point
	(*pCurrPath)(nRow, nCol) = 0;
	(*pPrevPath)(nRow, nCol) = -1;  // a number not in [0, 7]
	int nOff = m_nNeighborCnt / 2;  // offset corresponding to oppisite neighbor pair
	int nLabel = map[nRow][nCol], k;

PREVPIXEL:
	while ((k = (*pCurrPath)(nRow, nCol)) < m_nNeighborCnt)
	{
		(*pCurrPath)(nRow, nCol)++;

		// Exclude the coming path
		if (k == (*pPrevPath)(nRow, nCol))
			continue;

		int r = nRow + m_Neighbors[k].Row;
		int c = nCol + m_Neighbors[k].Col;
		if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
			continue;
		if (pMask && !(*pMask)(r, c))
			continue;
		if (map[r][c] != UNSUBMERGED)
			continue;  // already flooded
		if (img[r][c] > nWaterHeight)
			// exclude above flood pixel
			continue;

		// same segment
		map[r][c] = nLabel;
		if (IsBorder(map, r, c, nWidth, nHeight))
		{
			map[r][c] = MAPBORDER;
			continue;
		}

		(*pPrevPath)(r, c) = (k + nOff) % m_nNeighborCnt;
		(*pCurrPath)(r, c) = 0;
		nRow = r;
		nCol = c;
	}

	if ((k = (*pPrevPath)(nRow, nCol)) != -1)
	{
		// go back to previous path pixel
		nRow += m_Neighbors[k].Row;
		nCol += m_Neighbors[k].Col;
		goto PREVPIXEL;
	}
}
