/** \file
 *
 * \brief TMatrix class.
 *
 */

#ifndef TMatrix_Class
#define TMatrix_Class

#include <string>
#include <sstream>

#include "DataType.h"

/// Class for Matrix operations.
class TMatrix
{
public:
	/// Constructor
	TMatrix();
	/// Constructor
	TMatrix(int nRowCnt, int nColCnt);
	/// Constructor with data.
	TMatrix(double* pData, int nRowCnt, int nColCnt);
	/// Copy constructor.
	TMatrix(const TMatrix& m);
	/// Assignment operator
	TMatrix& operator=(const TMatrix& m);
public:
	/// data array itself
	double* m_pDataArray;
	/// allows row and column access
	double** m_ppData;
	/// number of rows
	int m_nRowCnt;
	/// number of columns.
	int m_nColCnt;

// Operations
public:
	/// Load data
	void LoadData(double* pData);
	/// Copy matrix
	BOOL CopyFrom(TMatrix* pMatrix);
	/// Set all elements to be zero
	void SetZero();

	/// Access the elements
	double& operator()(int nRow, int nCol);
	/// Arithmetic operation
	void operator*=(double dbVal);
    void operator/=(double dbVal);
    void operator-=(double dbVal);
    void operator+=(double dbVal);
	void operator+=(const TMatrix& A);
	void operator-=(const TMatrix& A);

	TMatrix& operator*(const double& dbVal);
    TMatrix& operator/(const double& dbVal);
    TMatrix& operator+(const double& dbVal);
    TMatrix& operator-(const double& dbVal);

    // Sum of all elements
    double Sum();

	/// Obtain the transpose of the matrix
	TMatrix* Transpose();
	BOOL Transpose(TMatrix* T);

	/// Compute the inverse of the matrix
	TMatrix* Inv();
	BOOL Inv(TMatrix* A);

	/// Compute the multiplication of current matrix and the specified matrix
	TMatrix* Mul(TMatrix* A);
	BOOL Mul(TMatrix* Result, TMatrix* A);

	/// Perform the Cholesky decomposition
	TMatrix* CholeskyDecompostion();
	BOOL CholeskyDecompostion(TMatrix* L);

	/// Obtain the determinant of the matrix
	double Determinant();

    	/// Obtain the eigenvalues and eigenvectors of the matrix, return the smallest eigenvalue
	double Eigen(double* EigVals, TMatrix* EigVecs);

	/// Obtain the smallest eigenvalue
	double Eigen();

	/// Exchange row
	void ExchangeRow(int nRow1, int nRow2);

	/// Converts this matrix to a single line string with rows concatenated one after the other.
	std::string ToString();

private:
	/// Allocate the memory for the matrix
	void AllocMemory(int nRowCnt, int nColCnt);

// Implementation
public:
	virtual ~TMatrix();

};

#endif

/////////////////////////////////////////////////////////////////////////////
