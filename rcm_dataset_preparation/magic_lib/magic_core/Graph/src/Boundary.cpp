// Boundary.cpp : implementation of the TBoundary class
//
// Class for region boundary struct
/////////////////////////////////////////////////////////////////////////////

#include <cstring>
#include "Boundary.h"

#include<iostream>


const static int NeighborCodeTable[3][3] = {
	{ -1, 1, -1}, { 0, 4, 2}, { -1, 3, -1} 
};

//
// Class of chain code. It shall be noted that we use a 4-neighbor system 
// plus the center site itself, and hence each site in the chain can be 
// represented by 3 bits. We use a byte with the lowest 3 bits for each 
// site, and the highest 3 bits for representing the boundary direction at
// each site
/////////////////////////////////////////////////////////////////////////////

// 
// Constructor
// 
TChainCode::TChainCode(POINTex StartPt, int nPtCnt)
{
	Init(StartPt, 0, nPtCnt);
}

// 
// Constructor
// 
TChainCode::TChainCode(POINTex StartPt, int nStartPtDir, int nPtCnt)
{
	Init(StartPt, nStartPtDir, nPtCnt);
}

// 
// Initialization
// 
void TChainCode::Init(POINTex StartPt, int nStartPtDir, int nPtCnt)
{
	m_StartPt = StartPt;
	m_EndPt = StartPt;
	m_nStartPtDir = nStartPtDir << 5;  // place into the highest 3 bits in the byte

	m_nTotalCnt = nPtCnt;
	m_pChainData = new BYTE[m_nTotalCnt];
	m_nPtCnt = 0;

	m_pNext = 0;
}

// 
// Destructor
// 
TChainCode::~TChainCode()
{
	if (m_pChainData) delete[] m_pChainData;
}

// 
// Compute the reverse chain
// 
void TChainCode::ReverseChain()
{
	// Update the site locations
	BYTE* pNewData = new BYTE[m_nPtCnt];
	int i;
	for (i = 0; i < m_nPtCnt; i ++)
	{
		if (m_pChainData[m_nPtCnt - 1 - i] & 0x0004)
			pNewData[i] = 4;  // center site
		else
			pNewData[i] = (m_pChainData[m_nPtCnt - 1 - i] + 2) & 0x0003;
	}

	// Update the directions
	for (i = 0; i < m_nPtCnt - 1; i ++)
		pNewData[i] |= m_pChainData[m_nPtCnt - 2 - i] & 0x00e0;
	pNewData[i] |= m_nStartPtDir;
	m_nStartPtDir = m_pChainData[m_nPtCnt - 1] & 0x00e0;

	// Preserve the new data
	POINTex Pt = m_EndPt;
	m_EndPt = m_StartPt;
	m_StartPt = Pt;

	delete[] m_pChainData;
	m_pChainData = pNewData;
	m_nTotalCnt = m_nPtCnt;
}

// 
// Get the first site
// 
void TChainCode::GetFirstSite(POINTex& Pt)
{
	int nDir;
	GetFirstSite(Pt, nDir);
}

// 
// Get the first site
// 
void TChainCode::GetFirstSite(POINTex& Pt, int& nDir)
{
	m_CurPt = m_StartPt;
	m_nCurIndex = 0;
	Pt = m_CurPt;
	nDir = m_nStartPtDir >> 5;
}

// 
// Get the next site
// 
BOOL TChainCode::GetNextSite(POINTex& Pt)
{
	int nDir;
	return GetNextSite(Pt, nDir);
}

// 
// Get the next site
// 
BOOL TChainCode::GetNextSite(POINTex& Pt, int& nDir)
{
	if (m_nCurIndex >= m_nPtCnt)
		return false;

	int j = m_pChainData[m_nCurIndex ++];
//	nDir = j >> 5;
//	j &= 0x0007;
    nDir = j; // Hmmmm.....
	m_CurPt.Row += neighbor5[j].Row;
	m_CurPt.Col += neighbor5[j].Col;
	Pt = m_CurPt;
	return true;
}

// 
// Check if two points can be concatenated by chain code representation.
// If so, return the corresponding chain code.
// 
// Return
//  -1: cannot be concatenated; Otherwise: the chain code
// 
int TChainCode::IsConcatenatable(POINTex CurPt, POINTex NextPt)
{
	int dr = NextPt.Row - CurPt.Row + 1;
	int dc = NextPt.Col - CurPt.Col + 1;
	if (dr < 0 || dr >= 3 || dc < 0 || dc >= 3)
		return -1;
	return NeighborCodeTable[dr][dc];
}

// 
// Add a site
// 
// Param in
//  Pt: the location of the site
//  nDir: the boundary direction at the site
// 
BOOL TChainCode::AddSite(POINTex Pt, int nDir)
{
	// Compute the code
	int code = IsConcatenatable(m_EndPt, Pt);
	if (code < 0)
	{
		code = IsConcatenatable(m_StartPt, Pt);
		if (code < 0)
			return false;
		ReverseChain();
	}

	// Check if need to expand memory
	if (m_nPtCnt >= m_nTotalCnt)
	{
		m_nTotalCnt = (int) (m_nTotalCnt * 1.1 + 4);
		BYTE* pNewData = new BYTE[m_nTotalCnt];
		memcpy(pNewData, m_pChainData, m_nPtCnt * sizeof(BYTE));
		delete[] m_pChainData;
		m_pChainData = pNewData;
	}

	m_pChainData[m_nPtCnt ++] = code | (nDir << 5);
	m_EndPt = Pt;
	return true;
}

// 
// Delete a site. The site has to be either the starting or ending
// site of the chain. Otherwise, return false
// 
BOOL TChainCode::DelSite(POINTex Pt)
{
	if (m_nPtCnt <= 0)  // the whole chain has less than two sites
		return false;

	int j;
	if (Pt.Row == m_StartPt.Row && Pt.Col == m_StartPt.Col)
	{
		m_nPtCnt --;
		j = m_pChainData[0];
		m_nStartPtDir = j >> 5;
		j &= 0x0007;
		m_StartPt.Row += neighbor5[j].Row;
		m_StartPt.Col += neighbor5[j].Col;
		memcpy(m_pChainData, &m_pChainData[1], m_nPtCnt * sizeof(BYTE));
		return true;
	}
	else if (Pt.Row == m_EndPt.Row && Pt.Col == m_EndPt.Col)
	{
		m_nPtCnt --;
		j = m_pChainData[m_nPtCnt] & 0x0007;
		m_EndPt.Row -= neighbor5[j].Row;
		m_EndPt.Col -= neighbor5[j].Col;
		return true;
	}
	return false;
}

//
// Merge data of the specifed chain with current
// 
BOOL TChainCode::MergeData(TChainCode* pChain)
{
	int code = IsConcatenatable(m_EndPt, pChain->m_StartPt);
	if (code < 0)
		return false;

	// Check if need to expand memory
	int nTotalCnt = m_nPtCnt + 1 + pChain->m_nPtCnt;
	if (nTotalCnt > m_nTotalCnt)
	{
		m_nTotalCnt = nTotalCnt;
		BYTE* pNewData = new BYTE[m_nTotalCnt];
		memcpy(pNewData, m_pChainData, m_nPtCnt * sizeof(BYTE));
		delete[] m_pChainData;  //QIN_MOD
		m_pChainData = pNewData;
	}

	// Merge data
	m_pChainData[m_nPtCnt ++] = code | pChain->m_nStartPtDir;
	memcpy(&m_pChainData[m_nPtCnt], pChain->m_pChainData,
	       pChain->m_nPtCnt * sizeof(BYTE));
	m_nPtCnt += pChain->m_nPtCnt;
	m_EndPt = pChain->m_EndPt;

	return true;
}


//
// Class for keys used for the sorting in AVL tree
/////////////////////////////////////////////////////////////////////////////

// 
// Constructor
// 
TBoundaryKey::TBoundaryKey(int nLabel1, int nLabel2)
	: TKey()
{
	if (nLabel1 > nLabel2)
	{
		m_nLabel1 = nLabel1;
		m_nLabel2 = nLabel2;
	}
	else
	{
		m_nLabel1 = nLabel2;
		m_nLabel2 = nLabel1;
	}
}

// 
// Compare if current key is less or larger than the specified one
// 
// Return
//  0: equal; -1: less; 1: larger
// 
int TBoundaryKey::Compare(TKey* pKey)
{
	TBoundaryKey* p = (TBoundaryKey*) pKey;
	if (m_nLabel1 < p->m_nLabel1)
		return -1;
	else if (m_nLabel1 > p->m_nLabel1)
		return 1;
	else
	{
		if (m_nLabel2 < p->m_nLabel2)
			return -1;
		else if (m_nLabel2 > p->m_nLabel2)
			return 1;
		else
			return 0;
	}
}

// 
// Clone
// 
TKey* TBoundaryKey::Clone()
{
	TBoundaryKey* pRet = new TBoundaryKey(m_nLabel1, m_nLabel2);
	return pRet;
}


//
// Class for wrapping each boundary segment
/////////////////////////////////////////////////////////////////////////////

// 
// Constructor
//
TBoundary::TBoundary(TGraph* pGraph, int nLabel1, int nLabel2, TChainCode* pChain, BOOL bOnGrid)
	: TEdge(pGraph), TAVLTreeNode(new TBoundaryKey(nLabel1, nLabel2))
{
	m_nLabel1 = nLabel1;
	m_nLabel2 = nLabel2;
	m_pChain = pChain;
	m_bActive = true;
	m_bOnGrid = bOnGrid;

	// Compute the length of the boundary
	m_nLength = 0;
	TChainCode* pCur = pChain;
	while (pCur)
	{
		m_nLength += pCur->m_nPtCnt + 1;
		pCur = pCur->m_pNext;
	}
}

// 
// Destructor
// 
TBoundary::~TBoundary()
{
	TChainCode* pCur = m_pChain;
	while (pCur)
	{
		m_pChain = m_pChain->m_pNext;
		delete pCur;
		pCur = m_pChain;
	}
}

// 
// Attach a chain code
// 
void TBoundary::AttachChain(TChainCode* pChain)
{
	// First check if the specified chain code can be merged
	int code;
	TChainCode* pCur = m_pChain;
	while (pCur)
	{
		if ((code = pCur->IsConcatenatable(pCur->m_EndPt, pChain->m_StartPt))
			>= 0)
		{
			pCur->MergeData(pChain);
			m_nLength += pChain->m_nPtCnt + 1;
			delete pChain;
			return;
		}
		if ((code = pCur->IsConcatenatable(pCur->m_EndPt, pChain->m_EndPt))
			>= 0)
		{
			pChain->ReverseChain();
			pCur->MergeData(pChain);
			m_nLength += pChain->m_nPtCnt + 1;
			delete pChain;
			return;
		}
		if ((code = pCur->IsConcatenatable(pCur->m_StartPt, pChain->m_StartPt))
			>= 0)
		{
			pCur->ReverseChain();
			pCur->MergeData(pChain);
			m_nLength += pChain->m_nPtCnt + 1;
			delete pChain;
			return;
		}
		if ((code = pCur->IsConcatenatable(pCur->m_StartPt, pChain->m_EndPt))
			>= 0)
		{
			pCur->ReverseChain();
			pChain->ReverseChain();
			pCur->MergeData(pChain);
			m_nLength += pChain->m_nPtCnt + 1;
			delete pChain;
			return;
		}

		pCur = pCur->m_pNext;
	}

	// Insert at the begining of the linked list
	pChain->m_pNext = m_pChain;
	m_pChain = pChain;
	m_nLength += pChain->m_nPtCnt + 1;
}

// 
// Detach a chain code from the linked list. Please note that the
// corresponding memory is not freed in this function
// 
void TBoundary::DetachChain(TChainCode* pChain)
{
	if (m_pChain == pChain)
	{
		m_pChain = m_pChain->m_pNext;
		m_nLength -= pChain->m_nPtCnt + 1;
	}
	else
	{
		TChainCode* pCur = m_pChain;
		while (pCur->m_pNext)
		{
			if (pCur->m_pNext == pChain)
			{
				pCur->m_pNext = pChain->m_pNext;
				m_nLength -= pChain->m_nPtCnt + 1;
				break;
			}
			pCur = pCur->m_pNext;
		}
	}
}

// 
// Add a site
// 
// Param in
//  Pt: the location of the site
//  nDir: the boundary direction at the site
// 
void TBoundary::AddSite(POINTex Pt, int nDir)
{
	TChainCode* pChain = m_pChain;
	while (pChain)
	{
		if (pChain->AddSite(Pt, nDir))
			break;
		pChain = pChain->m_pNext;
	}

	if (!pChain)  // cannot be attached to any existing chain, make a new chain
	{
		TChainCode* pNewChain = new TChainCode(Pt, nDir, 10);
		AttachChain(pNewChain);
	}
	else
		m_nLength ++;
}

// 
// Merge with a specified boundary
// 
void TBoundary::MergeWith(TEdge* pEdge)  //QIN_MOD_POINTER
{
	TBoundary* pBoundary = (TBoundary*) pEdge;
	TChainCode* pCur = pBoundary->m_pChain;
	TChainCode* pNext;
	while (pCur)
	{
		pNext = pCur->m_pNext;  // Note that this is a neccessary step: since pCur will be freed in AttachChain, its next must be first extracted
		// Merge the chain code
		pBoundary->DetachChain(pCur);
		AttachChain(pCur);

		pCur = pNext;
	}
}
