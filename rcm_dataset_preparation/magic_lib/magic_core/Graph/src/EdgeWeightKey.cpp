// EdgeWeightKey.cpp : implementation of the TEdgeWeightKey class
//
// Class for keys used for the sorting in AVL tree. The key is the
// merging edge weight
/////////////////////////////////////////////////////////////////////////////

#include "EdgeWeightKey.h"


// 
// Constructor
// 
TEdgeWeightKey::TEdgeWeightKey(double dbVal)
	: TKey()
{
	m_dbVal = dbVal;
}

// 
// Compare if current key is less or larger than the specified one
// 
// Return
//  0: equal; -1: less; 1: larger
// 
int TEdgeWeightKey::Compare(TKey* pKey)
{
	TEdgeWeightKey* p = (TEdgeWeightKey*) pKey;
	if (m_dbVal < p->m_dbVal)
		return -1;
	else if (m_dbVal > p->m_dbVal)
		return 1;
	else
		return 0;
}

// 
// Clone
// 
TKey* TEdgeWeightKey::Clone()
{
	TEdgeWeightKey* pRet = new TEdgeWeightKey(m_dbVal);
	return pRet;
}
