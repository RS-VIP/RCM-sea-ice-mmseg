// Matrix.cpp : implementation of the TMatrix class
/////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>
#include <cstring>
#include "Matrix.h"
#include "Eigen.h"

#define ZERO_THRESH 1e-12

// Default constructor
TMatrix::TMatrix() : m_nRowCnt(0), m_nColCnt(0), m_pDataArray(0), m_ppData(0)
{
}

//
// Constructor for Matrix
//
TMatrix::TMatrix(double* d, int r, int c)
{
	AllocMemory(r, c);
	LoadData(d);
}

//
// Constructor for Matrix
//
TMatrix::TMatrix(int r, int c)
{
	AllocMemory(r, c);
}


// Copy constructor.
TMatrix::TMatrix(const TMatrix& m)
{
	AllocMemory(m.m_nRowCnt, m.m_nColCnt);
	for (int r = m_nRowCnt - 1; r >= 0; r --)
		memcpy(m_ppData[r], m.m_ppData[r], m_nColCnt * sizeof(double));
}

// Assignment operator
TMatrix& TMatrix::operator=(const TMatrix& m)
{
	if (this != &m)
	{
		if (m_pDataArray) delete[] m_pDataArray;
		if ( m_ppData ) delete[] m_ppData;

		AllocMemory(m.m_nRowCnt, m.m_nColCnt);
		for (int r = m_nRowCnt - 1; r >= 0; r --)
			memcpy(m_ppData[r], m.m_ppData[r], m_nColCnt * sizeof(double));
	}
	return *this;
}

//
// Destructor
// 
TMatrix::~TMatrix()
{
	if (m_pDataArray) delete[] m_pDataArray;
	if ( m_ppData ) delete[] m_ppData;

}

//
// Allocate the memory for the matrix
// 
void TMatrix::AllocMemory(int r, int c)
{
	m_nRowCnt = r;
	m_nColCnt = c;
	m_pDataArray = new double[m_nRowCnt * m_nColCnt];
	m_ppData = new double*[m_nRowCnt];
	for (int i = 0; i < m_nRowCnt; i ++)
		m_ppData[i] = &m_pDataArray[i * m_nColCnt];
}

// 
// Load data into matrix
// 
void TMatrix::LoadData(double* d)
{
	memcpy(m_pDataArray, d, m_nRowCnt * m_nColCnt * sizeof(double));
}

// 
// Copy matrix
// 
BOOL TMatrix::CopyFrom(TMatrix* pMatrix)
{
	if (pMatrix->m_nRowCnt != m_nRowCnt || 
	    pMatrix->m_nColCnt != m_nColCnt)
		return false;
	for (int r = m_nRowCnt - 1; r >= 0; r --)
		memcpy(m_ppData[r], pMatrix->m_ppData[r], m_nColCnt * sizeof(double));
	return true;
}

// 
// Get the element of the matrix at the specified row and column
// 
double& TMatrix::operator()(int r, int c)
{
	return m_ppData[r][c];
}

// 
// Arithmetic operation
// 
void TMatrix::operator*=(double dbVal)
{
	for (int i = m_nRowCnt * m_nColCnt - 1; i >= 0; i --)
		m_pDataArray[i] *= dbVal;
}

void TMatrix::operator/=(double dbVal)
{
    for (int i = m_nRowCnt * m_nColCnt - 1; i >= 0; i --)
        m_pDataArray[i] /= dbVal;
}

void TMatrix::operator+=(double dbVal)
{
    for (int i = m_nRowCnt * m_nColCnt - 1; i >= 0; i --)
        m_pDataArray[i] += dbVal;
}

void TMatrix::operator-=(double dbVal)
{
    for (int i = m_nRowCnt * m_nColCnt - 1; i >= 0; i --)
        m_pDataArray[i] -= dbVal;
}

void TMatrix::operator+=(const TMatrix& A)
{
	double** ppData = A.m_ppData;
	for (int r = m_nRowCnt - 1; r >= 0; r --)
		for (int c = m_nColCnt - 1; c >= 0; c --)
	{
		m_ppData[r][c] += ppData[r][c];
	}
}

void TMatrix::operator-=(const TMatrix& A)
{
	double** ppData = A.m_ppData;
	for (int r = m_nRowCnt - 1; r >= 0; r --)
		for (int c = m_nColCnt - 1; c >= 0; c --)
	{
		m_ppData[r][c] -= ppData[r][c];
	}
}

TMatrix& TMatrix::operator*(const double& dbVal){
    TMatrix* pMat = new TMatrix(*this);
    *pMat *= dbVal;
    return *pMat;
}

TMatrix& TMatrix::operator/(const double& dbVal){
    TMatrix* pMat = new TMatrix(*this);
    *pMat /= dbVal;
    return *pMat;
}

TMatrix& TMatrix::operator+(const double& dbVal){
    TMatrix* pMat = new TMatrix(*this);
    *pMat += dbVal;
    return *pMat;
}

TMatrix& TMatrix::operator-(const double& dbVal){
    TMatrix* pMat = new TMatrix(*this);
    *pMat -= dbVal;
    return *pMat;
}

double TMatrix::Sum()
{
    double sum=0;
    for (int i = 0; i < m_nColCnt; i ++)
        for (int j = 0; j < m_nRowCnt; j ++)
            sum += m_ppData[j][i];
    return sum;
}

// 
// Set all elements to be zero
// 
void TMatrix::SetZero()
{
	memset(m_pDataArray, 0, m_nRowCnt * m_nColCnt * sizeof(double));
}

// 
// Exchange the rows
// 
void TMatrix::ExchangeRow(int r1, int r2)
{
	double* pSwap = m_ppData[r1];
	m_ppData[r1] = m_ppData[r2];
	m_ppData[r2] = pSwap;
}


// 
// Obtain the transpose of the matrix
// 
TMatrix* TMatrix::Transpose()
{
	TMatrix* T = new TMatrix(m_nColCnt, m_nRowCnt);
	Transpose(T);
	return T;
}

// 
// Obtain the transpose of the matrix
// 
BOOL TMatrix::Transpose(TMatrix* T)
{
	if (T->m_nRowCnt != m_nColCnt || T->m_nColCnt != m_nRowCnt)
		return false;

	for (int i = 0; i < m_nColCnt; i ++)
		for (int j = 0; j < m_nRowCnt; j ++)
			(*T)(i, j) = m_ppData[j][i];
	return true;
}

// 
// Compute the inverse of the matrix
// 
TMatrix* TMatrix::Inv()
{
	// Check if it is a square matrix
	if (m_nRowCnt != m_nColCnt)
		return 0;

	TMatrix* pRet = new TMatrix(m_nRowCnt, m_nColCnt);
	Inv(pRet);
	return pRet;
}

// 
// Compute the inverse of the matrix
// 
// Param out
//  I: the obtained inverse of current matrix
// 
BOOL TMatrix::Inv(TMatrix* I)
{
	// Check if it is a square matrix
	if (m_nRowCnt != m_nColCnt ||
	    I->m_nRowCnt != m_nRowCnt || I->m_nColCnt != m_nColCnt)
		return false;
	int nSize = m_nRowCnt;

	// The following compute the inverse using Gauss elimination
	TMatrix* A = new TMatrix(nSize, nSize);
	A->CopyFrom(this);
	int i, j;
	for (i = 0; i < nSize; i ++)
		for (j = 0; j < nSize; j ++)
	{
		// Initialize with the unit matrix
		if (i == j)
			(*I)(i, j) = 1;
		else
			(*I)(i, j) = 0;
	}

	double dbLargest;
	int nIndex;
	int r, c;
	for (i = 0; i < nSize; i ++)
	{
		// Search for largest (pivot) element
		dbLargest = fabs((*A)(i, i));
		nIndex = i;
		for (r = i + 1; r < nSize; r ++)
		{
			if (fabs((*A)(r, i)) > dbLargest)
			{
				dbLargest = fabs((*A)(r, i));
				nIndex = r;
			}
		}

		// Exchange the rows if not the smallest or zero
		if (nIndex != i)
		{
			A->ExchangeRow(i, nIndex);
			I->ExchangeRow(i, nIndex);
		}

		if (dbLargest < 0.000001)
			(*A)(i, i) = 0.000001;
		for (c = i + 1; c < nSize; c ++)
			(*A)(i, c) /= (*A)(i, i);
		for (c = 0; c < nSize; c ++)
			(*I)(i, c) /= (*A)(i, i);

		// Subtract so that current column in all rows are zero except current row
		for (r = 0; r < nSize; r ++)
		{
			if (r != i)
			{
				for (c = i + 1; c < nSize; c ++)
					(*A)(r, c) -= (*A)(r, i) * (*A)(i, c);
				for (c = 0; c < nSize; c ++)
					(*I)(r, c) -= (*A)(r, i) * (*I)(i, c);
			}
		}
	}

	delete A;
	return true;
}


// 
// Matrix multiplication
// 
TMatrix* TMatrix::Mul(TMatrix* A)
{
	if (m_nColCnt != A->m_nRowCnt)
		return 0;
	TMatrix* pResult = new TMatrix(m_nRowCnt, A->m_nColCnt);
	Mul(pResult, A);
	return pResult;
}


// 
// Matrix multiplication
// 
// Param in:
//  A: the specified matrix to be multiplied with current matrix
// 
// Param out:
//  pResult: obtained result
// 
BOOL TMatrix::Mul(TMatrix* pResult, TMatrix* A)
{
	if (m_nColCnt != A->m_nRowCnt || 
	    (pResult->m_nRowCnt != m_nRowCnt || pResult->m_nColCnt != A->m_nColCnt))
		return false;

	// Multiplication
	for (int i = 0; i < m_nRowCnt; i ++)
	{
		for (int j = 0; j < A->m_nColCnt; j ++)
		{
			(*pResult)(i, j) = 0;
			for (int k = 0; k < m_nColCnt; k ++)
				(*pResult)(i, j) += m_ppData[i][k] * (*A)(k, j);
		}
	}
	return true;
}

// 
// Perform the Cholesky decomposition
// 
// Return
//  the lower triangular L such that L*L'=A  
// 
TMatrix* TMatrix::CholeskyDecompostion()
{
	if (m_nRowCnt != m_nColCnt)
		return 0;

	TMatrix* L = new TMatrix(m_nRowCnt, m_nRowCnt);
	CholeskyDecompostion(L);
	return L;
}

// 
// Perform the Cholesky decomposition
// 
// Param out
//  the lower triangular L such that L*L'=A.
// 
BOOL TMatrix::CholeskyDecompostion(TMatrix* L)
{
	if (m_nRowCnt != m_nColCnt ||
	    L->m_nRowCnt != m_nRowCnt || L->m_nColCnt != m_nColCnt)
		return false;

	int i, j, k;
	double sum;
	int nSize = m_nRowCnt;
	for (i = 0; i < nSize; i ++)
	{
		for (j = i; j < nSize; j ++)
		{
			for (sum = m_ppData[i][j], k = i - 1; k >= 0; k --)
				sum -= (*L)(i, k) * (*L)(j, k);
			
			if (i == j) 
			{
				if (sum <= 0.0)
				{
					// printf("\nA is not positive definite!");
					return false;
				}
				else
					(*L)(i, i) = sqrt(sum);
			}
			else 
			{
				(*L)(j, i) = sum / (*L)(i, i);
			}
		}
	}
	for (i = 0; i < nSize; i ++)
		for (j = i + 1; j < nSize; j ++)
	{
		(*L)(i, j) = 0.0;
	}

	return true;
}

// 
// Obtain the determinant of the matrix. In this implementation,
// we assume the matrix is symmetric positive definite, and use
// Cholesky decomposition to compute the determinant. Therefore,
// it cannot work for all square matrice
// 
double TMatrix::Determinant()
{
	if (m_nRowCnt != m_nColCnt)
		return 0;

	TMatrix* L = new TMatrix(m_nRowCnt, m_nRowCnt);
	if (!CholeskyDecompostion(L))
	{
		delete L;
		return 0;
	}

	double d = (*L)(0, 0);
	for (int i = 1; i < m_nRowCnt; i ++)
		d *= (*L)(i, i);

	delete L;
	return (d * d);
}


double TMatrix::Eigen(double* EigVals, TMatrix* EigVecs)
{
	if ( m_nRowCnt != m_nColCnt)
	{
		printf("Error: Cannot take the eigendecomposition of a non-square matrix.");
		//_getchar_nolock();
		exit(1);
	}

	double dbDiff = 0.0;
	
	int i, j;

	for ( i = 0; i < m_nRowCnt; i++ )
		for ( j = 0; j < m_nColCnt; j++ )
			dbDiff += fabs( m_ppData[i][j] - m_ppData[i][j] );

	if ( ( dbDiff / ( m_nRowCnt * m_nColCnt ) ) > ZERO_THRESH )
	{
		printf("Error: The square real matrix to be eigen-decomposed should be symmetric!");
		//_getchar_nolock();
		exit(1);
	}

	double** A = new double*[m_nRowCnt + 1];
	for ( i = 0; i <= m_nRowCnt; i++ )
		A[i] = new double[m_nRowCnt + 1];
	for ( i = 1; i <= m_nRowCnt; i++ )
		for ( j = 1; j <= m_nRowCnt; j++ )
			A[i][j] = m_ppData[i - 1][j - 1];

	double* d = new double[m_nRowCnt + 1];
	double** V = new double*[m_nRowCnt + 1];
	for ( i = 0; i <= m_nRowCnt; i++ )
		V[i] = new double[m_nRowCnt + 1];
	
	eigen(A, m_nRowCnt, d, V);

	for ( i = 0; i < m_nRowCnt; i++ )
	{
		EigVals[i] = d[i + 1];
		for ( j = 0; j < m_nRowCnt; j++ )
			(*EigVecs)(i, j) = V[i + 1][j + 1];
	}

	for ( i = 0; i <= m_nRowCnt; i++ )
	{
		delete[] A[i];
		delete[] V[i];
	}
	delete[] A;
	delete[] d;
	delete[] V;

	double dbMin = 999999999999;
	for ( i = 0; i < m_nRowCnt; i++ )
		if ( EigVals[i] < dbMin )
			dbMin = EigVals[i];

	return dbMin;
}


double TMatrix::Eigen()
{
	if ( m_nRowCnt != m_nColCnt)
	{
		printf("Error: The matrix to be eigen-decomposed should be square real one!");
		//_getchar_nolock();
		exit(1);
	}

	double dbDiff = 0.0;
	
	int i, j;

	for ( i = 0; i < m_nRowCnt; i++ )
		for ( j = 0; j < m_nColCnt; j++ )
			dbDiff += fabs( m_ppData[i][j] - m_ppData[i][j] );

	if ( ( dbDiff / ( m_nRowCnt * m_nColCnt ) ) > ZERO_THRESH )
	{
		printf("Error: The square real matrix to be eigen-decomposed should be symmetric!");
		//_getchar_nolock();
		exit(1);
	}

	double** A = new double*[m_nRowCnt + 1];
	for ( i = 0; i <= m_nRowCnt; i++ )
		A[i] = new double[m_nRowCnt + 1];
	for ( i = 1; i <= m_nRowCnt; i++ )
		for ( j = 1; j <= m_nRowCnt; j++ )
			A[i][j] = m_ppData[i - 1][j - 1];

	double* d = new double[m_nRowCnt + 1];
	double** V = new double*[m_nRowCnt + 1];
	for ( i = 0; i <= m_nRowCnt; i++ )
		V[i] = new double[m_nRowCnt + 1];
	
	eigen(A, m_nRowCnt, d, V);

	double dbMin = 999999999999;
	for ( i = 1; i <= m_nRowCnt; i++ )
		if ( d[i] < dbMin )
			dbMin = d[i];

	for ( i = 0; i <= m_nRowCnt; i++ )
	{
		delete[] A[i];
		delete[] V[i];
	}
	delete[] A;
	delete[] d;
	delete[] V;

	return dbMin;
}

std::string TMatrix::ToString()
{
	std::string ret;

	// concatenate rows into string.
	for (int i = 0; i < m_nRowCnt * m_nColCnt; i++)
	{
		if (i != 0)
			ret += " ";
		ret += std::to_string(m_pDataArray[i]);
		//ret += StringConvertFunctions::ConvertToString(m_pDataArray[i], 0);
	}
	
	return ret;
}
