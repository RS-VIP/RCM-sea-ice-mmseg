/** \file
 *
 * \brief TSimpleClustering class.
 *
 * \author Peter Yu
 *
 */

#ifndef TSimpleClustering_Class
#define TSimpleClustering_Class

#include "SegmentationStatistics.h"
#include "GrayImage.h"

/// Base class for some (but not all) clustering methods.
/**
 * See subclasses for which specific clustering methods, there's quite a mix of
 * them.
 *
 */
class TSimpleClustering
{
public:

	TSimpleClustering(int nClusterCnt) : m_nClusterCnt(0), m_nClassPixelCnt(0)
	{
		if (nClusterCnt > 0)
		{
			m_nClusterCnt = nClusterCnt;
			AllocPixelCnts();
		}
	}

	/// Get segmentation statistics object that lists all the class statistics.
	/**
	 * The statistics will not be updated before being returned
	 */
	virtual TSegmentationStatistics* GetSegStats() = 0;

	/// Get segmentation statistics object that lists all the class statistics.
	/**
	 * The statistics will first be updated by the classes in pMap.
	 */
	virtual TSegmentationStatistics* GetSegStats(TGrayImage<int>* pMap) = 0;

protected:
	/// Number of clusters
	int m_nClusterCnt;

	/// number of pixels in each cluster
	int* m_nClassPixelCnt;

	void AllocPixelCnts()
	{
		if (m_nClusterCnt > 0)
		{
			if (m_nClassPixelCnt)
				delete [] m_nClassPixelCnt;

			m_nClassPixelCnt = new int[m_nClusterCnt];
		}
	}

// Implementation
public:
	virtual ~TSimpleClustering()
	{
		if (m_nClassPixelCnt)
			delete [] m_nClassPixelCnt;
	}

};

#endif
