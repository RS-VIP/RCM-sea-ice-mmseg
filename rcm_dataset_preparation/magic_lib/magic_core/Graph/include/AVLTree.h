/** \file
 *
 * \brief TAVLTree, TKey and TAVLTreeNode classes
 *
 */

#ifndef TAVLTree_Class
#define TAVLTree_Class


/// Key class for helping to sort AVL tree.
class TKey
{
public:
	TKey() {};

// Operations
public:
	virtual int Compare(TKey* pKey) = 0;
	virtual TKey* Clone() = 0;

};

/// AVL tree node class
class TAVLTreeNode
{
	friend class TAVLTree;

public:
	TAVLTreeNode(TKey* pKey);
	void Reset(TKey* pKey);
	TKey* GetKey(){return m_pKey->Clone();}

// Attributes
protected:
	TAVLTreeNode* m_pLeft;
	TAVLTreeNode* m_pRight;
	TAVLTreeNode* m_pParent;
	int m_nHeight;
	TKey* m_pKey;

// Implementation
public:
	~TAVLTreeNode();
};


/// AVL Tree class
class TAVLTree
{
public:
	TAVLTree();

public:
	/// root of the whole tree
	TAVLTreeNode* m_pRoot;

	/// upper limit of searching range
	TKey* m_pUpperKey;

	/// current node in searching
	TAVLTreeNode* m_pCur;

	/// number of nodes
	int m_nNodeCnt;

// Operations
public:
	/// insert a node into the AVL tree
	void Insert(TAVLTreeNode* pNode);

	/// find the node with the specified key value
	TAVLTreeNode* Find(TKey& Key);

	/// find the first node in the specified range
	TAVLTreeNode* FindFirstInRange(TKey& lKey, TKey& hKey);

	/// find the next node in the specified range
	TAVLTreeNode* FindNextInRange();

	/// get the most left child
	TAVLTreeNode* GetMostLeftChild(TAVLTreeNode* pNode);

	/// get the next node in the key order
	TAVLTreeNode* GetNextNode(TAVLTreeNode* pNode);

    /// get the size or rank of a node
	int GetSize(TAVLTreeNode* pNode);
	int GetRank(TAVLTreeNode* pNode);
	int GetRank(TKey* pKey);

	/// Resets root node and node count to zero
	void Reset();

protected:
	/// delete the subtree with the specified node as the tree root
	void DelTree(TAVLTreeNode* pTreeRoot); 

	/// insert into the AVL subtree with the specified node as the tree root
	TAVLTreeNode* Insert(TAVLTreeNode* pNode, TAVLTreeNode* pTreeRoot);

	/// perform a rotate between a node and its left child
	TAVLTreeNode* SingleRotateWithLeft(TAVLTreeNode* pNode);

	/// perform a rotate between a node and its right child
	TAVLTreeNode* SingleRotateWithRight(TAVLTreeNode* pNode);

	/// perform the left-right double rotation
	TAVLTreeNode* DoubleRotateWithLeft(TAVLTreeNode* pNode);

	/// perform the right-left double rotation
	TAVLTreeNode* DoubleRotateWithRight(TAVLTreeNode* pNode);

	/// return the height of the specified node
	int Height(TAVLTreeNode* pNode);  

// Implementation
public:
	~TAVLTree();
};

//// Example Key implementation
//class TTestKey: public TKey
//{
//public:
//    int m_nValue;
//    TTestKey(int value):TKey(){m_nValue = value;}
//
//// Operations
//public:
//    virtual int Compare(TKey* pKey){
//        if (m_nValue < (( TTestKey*) pKey)->m_nValue)
//            return 1;
//        else if (m_nValue > (( TTestKey*) pKey)->m_nValue)
//            return -1;
//        else return 0;
//    }
//    virtual TKey* Clone(){return new TTestKey(m_nValue);}
//
//};

#endif
