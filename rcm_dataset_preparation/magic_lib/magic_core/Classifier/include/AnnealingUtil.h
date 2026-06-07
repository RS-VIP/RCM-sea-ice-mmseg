/** \file
 *
 * \brief TAnnealingUtil class
 *
 */

#ifndef TAnnealingUtil_Class
#define TAnnealingUtil_Class

/// Utility class for simulated annealing
class TAnnealingUtil
{
public:
	TAnnealingUtil(double InitTemparature, double Ratio);

// Attributes
protected:
	/// temperature
	double T;  

	/// ratio for decreasing temerature
	double DecRatio;

	//double C;  // const for annealing schedule
	//int nIter;  // interation number

private:
	/// space for sampling use (allocated in class to save time)
	double* E;

	int m_nClusterCnt;

// Operations
public:
	const double GetTemparature() { return T; };
	void NextTemperature();

	/// Gibbs sampling using energies in supplied array, returns index of cluster sampled.
	int GibbsSampler(double dbEnergies[], int nClusterCnt);

// Implementation
public:
	virtual ~TAnnealingUtil();

};

#endif

