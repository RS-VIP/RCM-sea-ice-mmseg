/** \file
 *
 * \brief TVector class
 *
 */

#ifndef TVector_Class
#define TVector_Class

#include "Matrix.h"

/// Class for vector operations
class TVector : public TMatrix
{
public:
	TVector() : TMatrix(), m_bColumn(0), m_nElementCnt(0){}
	TVector(int nElementCnt, bool bColumn);
	TVector(double* pData, int nElementCnt, bool bColumn);
//	TVector(TVector& Other);

protected:

	/// Column matrix or not.
	bool m_bColumn;
	/// number of elements.
	int m_nElementCnt;

// Operations
public:
	/// Access the elements
	double& operator()(int nIndex);

	/// Access to the number of elements
	const int GetElementCnt() { return m_nElementCnt; };

	/// Obtain the transpose of the vector
	TVector* Transpose();

	/// Compute the multiplication of current vector and the specified vector
	double VecMul(TVector* A);

	/// Obtain the L-2 norm
	double Norm2();

	/// Get the maximum / minimum values
	double Max();
	double Min();

    /// Get the indices of the maximum / minimum values
    int Argmin();
    int Argmax();

// Implementation
public:
	virtual ~TVector();

};

#endif

/////////////////////////////////////////////////////////////////////////////
