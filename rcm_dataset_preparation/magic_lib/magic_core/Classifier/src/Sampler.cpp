// Sampler.cpp : implementation of the TSampler class
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include "Sampler.h"


//
// Constructor for Sampler
// 
TSampler::TSampler()
{
	// Initialize the Gaussian sampling array
	int i;
	for (i = 0; i < 512; i ++)
		m_dbGaussian[i] = exp(-i * i / (160 * 160 * 2.0));
	for (i = 1; i < 512; i ++)
		m_dbGaussian[i] = m_dbGaussian[i] + m_dbGaussian[i - 1];

	m_nScanTotalCnt = 0;
	m_nScanRemainCnt = 0;
	m_pScanRecord = 0;
}


//
// Destructor
// 
TSampler::~TSampler()
{
	if (m_pScanRecord)
		delete[] m_pScanRecord; //QIN_MOD
}


// 
// Uniformly distributed sampling
//
double TSampler::UniformSampling(double dbMean, double dbStd)
{
	double dbLimit = sqrt(3.0) * dbStd; //QIN_WHY
	double dbRandNum = rand() * dbLimit * 2 / (RAND_MAX + 1.0);
	return dbRandNum - dbLimit + dbMean;
}

// 
// Gaussian distributed sampling //QIN_MOD
// 
double TSampler::GaussianSampling(double dbMean, double dbStd)
{
	double dbRandNum = rand() * m_dbGaussian[511] * 2 / (RAND_MAX + 1.0);
	double dbVal = fabs(dbRandNum - m_dbGaussian[511]);
	int l = 0, r = 511, m;
	while (r - l > 1)                     //QIN_MOD
	{
		m = (r + l) / 2;
		if (dbVal <= m_dbGaussian[m])
			r = m;
		else
			l = m;
	}
	if (dbRandNum < m_dbGaussian[511])
		r = - r;

	return r * dbStd / 160 + dbMean;
}

// 
// Uniformly randomly distributed scan
// 
bool TSampler::InitScan(int nScanTotalCnt)
{
	if (m_pScanRecord && (m_nScanTotalCnt != nScanTotalCnt))
	{
		delete[] m_pScanRecord; //QIN_MOD
		m_pScanRecord = 0;
	}

	if ((m_nScanTotalCnt = nScanTotalCnt) <= 0)
		return false;

	m_nScanRemainCnt = m_nScanTotalCnt;
	if (!m_pScanRecord)
		m_pScanRecord = new int[m_nScanTotalCnt];
	for (int i = m_nScanTotalCnt - 1; i >= 0; i --)
		m_pScanRecord[i] = i;

	return true;
}

// 
// Get the next scan
// 
int TSampler::GetNextScan()
{
	if (m_nScanRemainCnt <= 0)
		return -1;

	double dbRandNum = ((double) rand()) / (RAND_MAX + 1.0) * m_nScanRemainCnt;
	int nIndex = (int) dbRandNum;
	int nRet = m_pScanRecord[nIndex];
	m_nScanRemainCnt --;
	m_pScanRecord[nIndex] = m_pScanRecord[m_nScanRemainCnt];
	return nRet;
}

