/** \file
 *
 * \brief TBoundary, TChainCode and TBoundaryKey classes
 *
 */

#ifndef TBoundary_Class
#define TBoundary_Class

#include "Graph.h"
#include "AVLTree.h"
#include "DataType.h"


/// Chain code representation of region boundaries.
class TChainCode
{
public:
	TChainCode(POINTex StartPt, int nPtCnt);
	TChainCode(POINTex StartPt, int nStartPtDir, int nPtCnt);

// Attributes
public:
	/// starting position of the chain
	POINTex m_StartPt;  

	// ending position of the chain
	POINTex m_EndPt; 

	// boundary direction symbols for the starting position.
	/**
	 * The directions are defined by:
	 * -------
	 * |1|2|3|
	 * -------
     * |0| |4|
	 * -------
	 * |7|6|5|
	 * ------- 
	 */
	BYTE m_nStartPtDir;

	/// chain data
	BYTE* m_pChainData;

	/// number of points in the chain, excluding the starting point
	int m_nPtCnt;

	/// number of total elements that the chain can contain
	int m_nTotalCnt;

	/// current point
	POINTex m_CurPt;

	/// index of current point
	int m_nCurIndex;

	/// for constructing linked list
	TChainCode* m_pNext;

// Operations
private:
	/// initialization
	void Init(POINTex StartPt, int nStartPtDir, int nPtCnt);  

public:

	/// Add a site
	BOOL AddSite(POINTex Pt, int nDir);

	/// Delete a site. 
	/**
	 * The site has to be either the starting or ending site of the chain.
	 * Otherwise, return false
	 */
	BOOL DelSite(POINTex Pt);

	/// Check if two points can be concatenated by chain code representation, if so, return the corresponding chain code.
	/** 
	 * returns -1 if cannot be concatenated.
	 */
	int IsConcatenatable(POINTex CurPt, POINTex NextPt);

	void GetFirstSite(POINTex& Pt);
	void GetFirstSite(POINTex& Pt, int& nDir);
	BOOL GetNextSite(POINTex& Pt);
	BOOL GetNextSite(POINTex& Pt, int& nDir);

	/// Compute the reverse chain
	void ReverseChain();

	/// Merge data of the specifed chain with current
	BOOL MergeData(TChainCode* pChain);

// Implementation
public:
	virtual ~TChainCode();

};


/// class for keys used for the sorting in AVL tree
class TBoundaryKey : public TKey
{
public:
	TBoundaryKey(int nLabel1, int nLabel2);

// Attribute
private:
	int m_nLabel1, m_nLabel2;

// Operations
public:
	virtual int Compare(TKey* pKey);
	virtual TKey* Clone();

};


/// Class for wrapping each boundary segment
class TBoundary : public TEdge, public TAVLTreeNode
{
public:
	TBoundary(TGraph* pGraph, int nLabel1, int nLabel2, TChainCode* pChain, BOOL bOnGrid);

// Attributes
public:
	/// boundary sites in chain code
	TChainCode* m_pChain;

	/// length of boundary
	int m_nLength;

	/// flag indicating if the regions connected by the boundary are active (not broken)
	BOOL m_bActive;

	/// flag indicating if the boundary sites are on grid or in between grid sites
	BOOL m_bOnGrid;

// Operations
public:
	/// attach a chain
	virtual void AttachChain(TChainCode* pChain);

	/// detach a chain
	virtual void DetachChain(TChainCode* pChain);

	/// add a site
	virtual void AddSite(POINTex Pt, int nDir);

	/// Merge with a specified boundary
	virtual void MergeWith(TEdge* pEdge);

// Implementation
public:
	virtual ~TBoundary();

};

 
#endif

/////////////////////////////////////////////////////////////////////////////
