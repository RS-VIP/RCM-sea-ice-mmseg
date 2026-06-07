/** \file
 *
 * \brief TSampler class
 *
 */

#ifndef TSampler_Class
#define TSampler_Class

/// Utility class for sampling.
class TSampler
{
public:
	TSampler();

// Attributes
protected:
	/// For Gaussian sampling
	double m_dbGaussian[512];

	/// Number of total random scans to be drawn
	int m_nScanTotalCnt;
	/// Number of remaining random scans to be drawn
	int m_nScanRemainCnt;
	/// Record for the scans
	int* m_pScanRecord;

// Operations
public:
	/// Uniformly distributed sampling
	double UniformSampling(double dbMean, double dbStd);
	/// Gaussian distributed sampling
	double GaussianSampling(double dbMean, double dbStd);

	/// Uniformly distributed scan order
	bool InitScan(int nScanTotalCnt);

	/// Get the next scan
	int GetNextScan();

// Implementation
public:
	virtual ~TSampler();

};

#endif

/////////////////////////////////////////////////////////////////////////////
