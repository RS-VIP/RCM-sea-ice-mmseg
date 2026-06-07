/** \file
 *
 * \brief TEdgeWeightKey class
 *
 */

#ifndef TEdgeWeightKey_Class
#define TEdgeWeightKey_Class

#include "AVLTree.h"

/// Class for keys used for the sorting in AVL tree
/**
 * The key is the merging edge weight.
 */
class TEdgeWeightKey : public TKey
{
public:
	TEdgeWeightKey(double dbVal);

// Attribute
private:
	double m_dbVal;

// Operations
public:
	virtual int Compare(TKey* pKey);
	virtual TKey* Clone();

};

 
#endif
