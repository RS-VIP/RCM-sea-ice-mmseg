// TAVLTree.cpp : implementation of the TAVLTree class
//
// Class for AVL tree
/////////////////////////////////////////////////////////////////////////////

#include "AVLTree.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))


//
// Class for AVL tree node
/////////////////////////////////////////////////////////////////////////////

// 
// Constructor
// 
TAVLTreeNode::TAVLTreeNode(TKey* pKey)
{
	m_pKey = 0;
	Reset(pKey);
}

// 
// Destructor
// 
TAVLTreeNode::~TAVLTreeNode()
{
	if (m_pKey)
		delete m_pKey;
}

// 
// Reset
// 
void TAVLTreeNode::Reset(TKey* pKey)
{
	m_pLeft = 0;
	m_pRight = 0;
	m_pParent = 0;
	if (m_pKey)
		delete m_pKey;
	m_pKey = pKey;
}

//
// Class for AVL tree
/////////////////////////////////////////////////////////////////////////////

// 
// Constructor
// 
TAVLTree::TAVLTree()
{
	m_pRoot = 0;
	m_pUpperKey = 0;
	m_nNodeCnt = 0;
}

//
// Destructor
// 
TAVLTree::~TAVLTree()
{
	DelTree(m_pRoot);
	if (m_pUpperKey)
		delete m_pUpperKey;
}

// 
// Reset
// 
void TAVLTree::Reset()
{
	m_pRoot = 0;
	m_nNodeCnt = 0;
}

//
// Delete the subtree with the specified node as the tree root
// 
void TAVLTree::DelTree(TAVLTreeNode* pTreeRoot)
{
	if (pTreeRoot)  // not null node
	{
		DelTree(pTreeRoot->m_pLeft);
		DelTree(pTreeRoot->m_pRight);
		delete(pTreeRoot);
	}
}

// 
// Return the height of the specified node
// 
int TAVLTree::Height(TAVLTreeNode* pNode)
{
	if (!pNode)
		return -1;
	else
		return pNode->m_nHeight;
}

//
// Insert a node into the AVL tree
// 
void TAVLTree::Insert(TAVLTreeNode* pNode)
{
	m_pRoot = Insert(pNode, m_pRoot);
	m_pRoot->m_pParent = 0;
	m_nNodeCnt ++;
}

// 
// Insert a node into the AVL subtree with the specified node as the tree root
// 
TAVLTreeNode* TAVLTree::Insert(TAVLTreeNode* pNode, TAVLTreeNode* pTreeRoot)
{
	if (!pTreeRoot)  // empty subtree
	{
		pTreeRoot = pNode;
		pTreeRoot->m_nHeight = 0;
	}
	else
	{
		if (pTreeRoot->m_pKey->Compare(pNode->m_pKey) > 0)  // the node to be inserted is smaller, insert left
		{
			pTreeRoot->m_pLeft = Insert(pNode, pTreeRoot->m_pLeft);
			if (Height(pTreeRoot->m_pLeft) - Height(pTreeRoot->m_pRight) == 2)
			{
				if (pTreeRoot->m_pLeft->m_pKey->Compare(pNode->m_pKey) > 0)
					pTreeRoot = SingleRotateWithLeft(pTreeRoot);
				else
					pTreeRoot = DoubleRotateWithLeft(pTreeRoot);
			}
			else
				pTreeRoot->m_pLeft->m_pParent = pTreeRoot;
		}
		else
		{
			pTreeRoot->m_pRight = Insert(pNode, pTreeRoot->m_pRight);
			if (Height(pTreeRoot->m_pRight) - Height(pTreeRoot->m_pLeft) == 2)
			{
				if (pTreeRoot->m_pRight->m_pKey->Compare(pNode->m_pKey) <= 0)
					pTreeRoot = SingleRotateWithRight(pTreeRoot);
				else
					pTreeRoot = DoubleRotateWithRight(pTreeRoot);
			}
			else
				pTreeRoot->m_pRight->m_pParent = pTreeRoot;
		}
	}

	pTreeRoot->m_nHeight = 1 + MAX(Height(pTreeRoot->m_pLeft), Height(pTreeRoot->m_pRight));
	return pTreeRoot;
}

// 
// Perform a rotate between a node and its left child, update heights,
// then return new root
// This function can be called only if the specified node has a left child
// 
TAVLTreeNode* TAVLTree::SingleRotateWithLeft(TAVLTreeNode* pNode)
{
	TAVLTreeNode* pLeft = pNode->m_pLeft;
	pNode->m_pLeft = pLeft->m_pRight;
	if (pLeft->m_pRight)
		pLeft->m_pRight->m_pParent = pNode;
	pLeft->m_pRight = pNode;
	pNode->m_pParent = pLeft;

	pNode->m_nHeight = MAX(Height(pNode->m_pLeft), Height(pNode->m_pRight)) + 1;
	pLeft->m_nHeight = MAX(Height(pLeft->m_pLeft), pNode->m_nHeight) + 1;

	return pLeft;  // new root
}

// 
// Perform a rotate between a node and its right child, update heights,
// then return new root
// This function can be called only if the specified node has a right child
// 
TAVLTreeNode* TAVLTree::SingleRotateWithRight(TAVLTreeNode* pNode)
{
	TAVLTreeNode* pRight = pNode->m_pRight;
	pNode->m_pRight = pRight->m_pLeft;
	if (pRight->m_pLeft)
		pRight->m_pLeft->m_pParent = pNode;
	pRight->m_pLeft = pNode;
	pNode->m_pParent = pRight;

	pNode->m_nHeight = MAX(Height(pNode->m_pLeft), Height(pNode->m_pRight)) + 1;
	pRight->m_nHeight = MAX(Height(pRight->m_pRight), pNode->m_nHeight) + 1;

	return pRight;  // new root
}

// 
// Perform the right-left double rotation
// 
TAVLTreeNode* TAVLTree::DoubleRotateWithLeft(TAVLTreeNode* pNode)
{
	pNode->m_pLeft = SingleRotateWithRight(pNode->m_pLeft);
	return SingleRotateWithLeft(pNode);
}

// 
// Perform the left-right double rotation
// 
TAVLTreeNode* TAVLTree::DoubleRotateWithRight(TAVLTreeNode* pNode)
{
	pNode->m_pRight = SingleRotateWithLeft(pNode->m_pRight);
	return SingleRotateWithRight(pNode);
}

// 
// Find the node with the specified key value
// 
TAVLTreeNode* TAVLTree::Find(TKey& Key)
{
	int nCompare;
	TAVLTreeNode* pNode = m_pRoot;
	while (pNode)
	{
		nCompare = pNode->m_pKey->Compare(&Key);
		if (nCompare > 0)  // search left
			pNode = pNode->m_pLeft;
		else if (nCompare < 0)  // search right
			pNode = pNode->m_pRight;
		else
			break;
	}

	return pNode;
}

// 
// Find the first node in the specified range
// 
TAVLTreeNode* TAVLTree::FindFirstInRange(TKey& lKey, TKey& hKey)
{
	if (!m_pRoot)  // empty tree
		return 0;

	if (m_pUpperKey)
		delete m_pUpperKey;  // free previous upper limit key
	m_pUpperKey = hKey.Clone();

	// Find the first node that is equal or greater than lKey
	int nCompare;
	TAVLTreeNode* pNode = m_pRoot;
	m_pCur = 0;
	do {
		nCompare = pNode->m_pKey->Compare(&lKey);
		if (nCompare > 0)  // current greater than lkey, preserve current and continue to search left
		{
			m_pCur = pNode;
			pNode = pNode->m_pLeft;
		}
		else if (nCompare < 0)  // current smaller, continue to search right
			pNode = pNode->m_pRight;
		else  // equal
		{
			m_pCur = pNode;
			return m_pCur;
		}
	} while (pNode);

	if (m_pCur && m_pCur->m_pKey->Compare(&hKey) <= 0)
		return m_pCur;
	else
		return 0;
}

// 
// Find the next node in the specified range
// 
TAVLTreeNode* TAVLTree::FindNextInRange()
{
	if (m_pCur = GetNextNode(m_pCur))
	{
		if (m_pCur->m_pKey->Compare(m_pUpperKey) <= 0)
			return m_pCur;
	}
	return 0;
}

// 
// Get the most left child
// 
TAVLTreeNode* TAVLTree::GetMostLeftChild(TAVLTreeNode* pNode)
{
	while (pNode->m_pLeft)
		pNode = pNode->m_pLeft;
	return pNode;
}

// 
// Get the next node in the key order
// 
TAVLTreeNode* TAVLTree::GetNextNode(TAVLTreeNode* pNode)
{
	// The node has right child, select the deepest left child of its right child
	if (pNode->m_pRight)
		return GetMostLeftChild(pNode->m_pRight);

	// The node has no right child, go to the position where the path begins
	// as a left child of the corresponding parent
	while (pNode->m_pParent && pNode == pNode->m_pParent->m_pRight)
		pNode = pNode->m_pParent;

	// The node has no right child, but is left child of its parent, select its parent
	return pNode->m_pParent;
}

/// FUNCTIONS ADDED BY MAX ///

/// Get the size of the subtree whose root is pNode
int TAVLTree::GetSize(TAVLTreeNode *pNode) {
    if (!pNode)
        return 0;
    else
        return GetSize(pNode->m_pLeft) + GetSize(pNode->m_pRight) + 1;
}

/// Get the rank of node pNode within the AVL tree (the number of nodes with smaller keys)
int TAVLTree::GetRank(TAVLTreeNode* pNode) {
    int rank = 0, nCompare;
    TAVLTreeNode* pCur = m_pRoot;

    while(pCur){

        if (pNode->m_pKey->Compare(pCur->m_pKey) == 0){ // matching key found.
            break;
        }

        nCompare = pCur->m_pKey->Compare(pNode->m_pKey);
        if (nCompare < 0){
            pCur = pCur->m_pRight;
        } else {
            rank += GetSize(pCur->m_pRight) + 1;
            pCur = pCur->m_pLeft;
        }
    }

    if (!pCur){
        return -1; // Key was not found
    } else {
        rank += GetSize(pCur->m_pRight);
        return rank;
    }
}

/// Get the rank of node pNode within the AVL tree (the number of nodes with smaller keys)
int TAVLTree::GetRank(TKey* pKey) {
    int rank = 0, nCompare;
    TAVLTreeNode* pCur = m_pRoot;

    while(pCur){

        if (pKey->Compare(pCur->m_pKey) == 0){ // matching key found.
            break;
        }

        nCompare = pCur->m_pKey->Compare(pKey);
        if (nCompare < 0){
            pCur = pCur->m_pRight;
        } else {
            rank += GetSize(pCur->m_pRight) + 1;
            pCur = pCur->m_pLeft;
        }
    }

    if (!pCur){
        return -1; // Key was not found
    } else {
        rank += GetSize(pCur->m_pRight);
        return rank;
    }
}