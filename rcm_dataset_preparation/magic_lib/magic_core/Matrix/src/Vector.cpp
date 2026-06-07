// Vector.cpp : implementation of the TVector class
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Vector.h"
//#include "ProgressEvent.h"


//
// Constructor for vector
//
TVector::TVector(double* pData, int nElementCnt, bool bColumn)
	: TMatrix(pData, (bColumn)?nElementCnt:1, (bColumn)?1:nElementCnt)
{
	m_bColumn = bColumn;
	m_nElementCnt = nElementCnt;
}

//
// Constructor for vector
//
TVector::TVector(int nElementCnt, bool bColumn)
	: TMatrix((bColumn)?nElementCnt:1, (bColumn)?1:nElementCnt)
{
	m_bColumn = bColumn;
	m_nElementCnt = nElementCnt;
}

////
//// Copy Constructor for vector
////
//TVector::TVector(TVector& Other)
//        : TMatrix(Other)
//{
//    m_bColumn = Other.m_bColumn;
//    m_nElementCnt = Other.m_nElementCnt;
//}

//
// Destructor
// 
TVector::~TVector()
{
}

// 
// Get the element of the vector at the specified row and column
// 
double& TVector::operator()(int nIndex)
{
	if (m_bColumn)
		return m_ppData[nIndex][0];
	else
		return m_ppData[0][nIndex];
}

// 
// Obtain the transpose of the vector
// 
TVector* TVector::Transpose()
{
	TVector* T = new TVector(m_nElementCnt, !m_bColumn);
	TMatrix::Transpose(T);
	return T;
}


// 
// Vector multiplication
// 
double TVector::VecMul(TVector* A)
{
//	if (m_nElementCnt != A->m_nElementCnt)
//		return 0;
	double dbRet = 0;
	for (int k = 0; k < m_nElementCnt; k ++)
		dbRet += m_pDataArray[k] * A->m_pDataArray[k];
	return dbRet;
}

// 
// Obtain the L-2 norm of the vector.
// 
double TVector::Norm2()
{
	return sqrt(VecMul(this));
}

//
// Obtain the minimum value of the vector.
//
double TVector::Min()
{
    double val = m_pDataArray[0];
    for (int k = 0; k < m_nElementCnt; k ++)
        val = m_pDataArray[k] < val ? m_pDataArray[k] : val;
    return val;
}

//
// Obtain the maximum value of the vector.
//
double TVector::Max()
{
    double val = m_pDataArray[0];
    for (int k = 0; k < m_nElementCnt; k ++)
        val = m_pDataArray[k] > val ? m_pDataArray[k] : val;
    return val;
}


//
// Obtain the index of the minimum value of the vector.
//
int TVector::Argmin() {
    double val = m_pDataArray[0];
    int idx = 0;
    for (int k = 0; k < m_nElementCnt; k ++)
        if (m_pDataArray[k] < val){ val=m_pDataArray[k]; idx=k; }
    return idx;
}

//
// Obtain the index of the maximum value of the vector.
//
int TVector::Argmax() {
    double val = m_pDataArray[0];
    int idx = 0;
    for (int k = 0; k < m_nElementCnt; k ++)
        if (m_pDataArray[k] > val){ val=m_pDataArray[k]; idx=k; }
    return idx;
}
