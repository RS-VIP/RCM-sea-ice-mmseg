/** \file
 *
 * \brief TSegmentationStatistics class.
 *
 */

#pragma once

#include <vector>

#include "Matrix.h"
#include "Vector.h"

using namespace std;
//#include "PolarimetricCovarianceMatrix.h"

/// Stores statistics for a segmentation
class TSegmentationStatistics
{
public:
	/// Type of statistic.
	enum StatisticType { GAUSSIAN, POLARIMETRIC_WISHART };

	TSegmentationStatistics() : m_nClusterCnt(0), m_nPolygonNumber(0)
	{
	}

	TSegmentationStatistics(int nClusterCnt) : m_nClusterCnt(0), m_nPolygonNumber(0)
	{
		if (nClusterCnt > 0)
		{
			m_nClusterCnt = nClusterCnt;
			m_vPixelCounts.resize(m_nClusterCnt, 0);
		}
	}

	virtual TSegmentationStatistics* Clone() = 0;

	/// Sets pixel counts from an int array with nClusterCnt elements.
	void SetPixelCounts(long* nPixelCnt)
	{
		for (int i = 0; i < m_nClusterCnt; i++)
			m_vPixelCounts[i] = nPixelCnt[i];
	}

	/// Gets the type of this statistic set.
	virtual StatisticType GetType() = 0;

	/// Number of clusters
	int m_nClusterCnt;

	/// number of pixels in each class.
	vector<long> m_vPixelCounts;

	/// stores polygon number (for GUI's all polygons mode; set by GUI itself)
	int m_nPolygonNumber;

public:
	virtual ~TSegmentationStatistics()
	{
	}
};

/// Stores statistics for a Gaussian statistics segmentation
class TGaussianSegmentationStatistics : public TSegmentationStatistics
{
public:

	TGaussianSegmentationStatistics() : TSegmentationStatistics(0)
	{

	}

	TGaussianSegmentationStatistics(int nClusterCnt) : TSegmentationStatistics(nClusterCnt)
	{
		if (m_nClusterCnt > 0)
		{
			m_vMu.resize(m_nClusterCnt);
			m_vCov.resize(m_nClusterCnt);
		}
	}

	/// Has the effect of a virtual copy constructor.
	virtual TGaussianSegmentationStatistics* Clone()
	{
		return new TGaussianSegmentationStatistics(*this);
	};

	/// Gets the type of this statistic set.
	virtual StatisticType GetType()
	{
		return TSegmentationStatistics::GAUSSIAN;
	}

	/// means of the clusters
	vector<TVector> m_vMu;
	/// covariance matrix of the clusters	
	vector<TMatrix> m_vCov;

public:

	virtual ~TGaussianSegmentationStatistics()
	{
	}	

};

///// Stores statistics for a Polarimetric statistics segmentation
//class TPolarimetricSegmentationStatistics : public TSegmentationStatistics
//{
//public:
//
//	TPolarimetricSegmentationStatistics() : TSegmentationStatistics(0)
//	{
//
//	}
//
//	TPolarimetricSegmentationStatistics(int nClusterCnt) : TSegmentationStatistics(nClusterCnt)
//	{
//		if (m_nClusterCnt > 0)
//		{
//			m_vCov.resize(m_nClusterCnt);
//		}
//	}
//
//	/// Has the effect of a virtual copy constructor.
//	virtual TPolarimetricSegmentationStatistics* Clone()
//	{
//		return new TPolarimetricSegmentationStatistics(*this);
//	};
//
//	/// Gets the type of this statistic set.
//	virtual StatisticType GetType()
//	{
//		return TSegmentationStatistics::POLARIMETRIC_WISHART;
//	}
//
//	/// mean polarimetric covariance of the clusters
//	vector<TPolarimetricCovarianceMatrix> m_vCov;
//
//public:
//	virtual ~TPolarimetricSegmentationStatistics()
//	{
//	}
//};
