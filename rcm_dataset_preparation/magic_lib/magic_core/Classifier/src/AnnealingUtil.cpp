// AnnealingUtil.cpp : implementation of the TAnnealingUtil class
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "GrayImage.h"
#include "AnnealingUtil.h"

#define Pi 3.1415926535897


//
// Constructor for AnnealingUtil
//
TAnnealingUtil::TAnnealingUtil(double InitTemparature, double Ratio)
{
	T = InitTemparature;
	DecRatio = Ratio;
//	C = T * log(2);
//	nIter = 1;

	E = 0;  // null array
	m_nClusterCnt = 0;
}


//
// Destructor
//
TAnnealingUtil::~TAnnealingUtil()
{
	if (E)
		delete[] E; //QIN_MOD
}


// 
// Compute the temperature for the next iteration
//
void TAnnealingUtil::NextTemperature()
{
	T = T * DecRatio;
//	nIter ++;
//	T = C / log(1 + nIter);
}


// 
// Gibbs sampler
//
// Param in
//  Energies: the energies
//  ClusterCnt: number of clusters
// 
// Return
//  the index of cluster sampled [0 ... nClusterCnt-1]
//
int TAnnealingUtil::GibbsSampler(double dbEnergies[], int nClusterCnt)
{
	if (E && m_nClusterCnt != nClusterCnt)
	{
		delete[] E; //QIN_MOD
		E = 0;
	}
	if (!E)
	{
		m_nClusterCnt = nClusterCnt;
		E = new double[m_nClusterCnt];
	}

	double dbMin = 100000000;  // a large enough number 10^8
	int k;
	for (k = nClusterCnt - 1; k >= 0; k --)
	{
		if (dbMin > dbEnergies[k])
			dbMin = dbEnergies[k];
	}

	// Sampling
	double dbTotal = 0;
	for (k = nClusterCnt - 1; k >= 0; k --)
	{
		// Subtract a const so that not all numbers after exp() will be near zero
		E[k] = exp(-(dbEnergies[k]-dbMin)/T);  
		dbTotal += E[k];
	}

	double RandNum = rand() * dbTotal / (RAND_MAX + 1.0);
	for (k = 0; k < nClusterCnt; k ++)
	{
		if (RandNum < E[k])
			break;
		E[k + 1] += E[k];
	}
	//Fix for crashing that occurs rarely when RandNum and E[k] are both close to 0.
	if (k == nClusterCnt){
		k --;
	}
	return k;
}

