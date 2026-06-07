// Graph.cpp : implementation of the TGraph class
/////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <time.h>
#include "Graph.h"
#include "Boundary.h"

// initial max. number of neighbours for each node. This is not a hard constraint, the
// size of the neighbour list is expanded if necessary.
#define MAX_NEIGHBORS 4


//
// Base class for edges between graph vertices
/////////////////////////////////////////////////////////////////////////////

//
// Constructor for graph edge
// 
TEdge::TEdge(TGraph* pGraph)
{
	Init(0, 0);
	m_pGraph = pGraph;
}

//
// Constructor for graph edge
// 
TEdge::TEdge(TVertex* pNeighbor1, TVertex* pNeighbor2)
{
	Init(pNeighbor1, pNeighbor2);
	m_pGraph = pNeighbor1->m_pGraph;
}

// 
// Initialization
// 
void TEdge::Init(TVertex* pNeighbor1, TVertex* pNeighbor2)
{
	m_nLabel1 = -1;
	m_nLabel2 = -1;
	m_nIndex = 0;
	SetVertexNeighbor(pNeighbor1, pNeighbor2);

	m_pPrev = 0;
	m_pNext = 0;
}

// 
// Set the neighboring vertice
// 
void TEdge::SetVertexNeighbor(TVertex* pNeighbor1, TVertex* pNeighbor2)
{
	m_pNeighbor1 = pNeighbor1;
	m_pNeighbor2 = pNeighbor2;
	if (m_pNeighbor1)
		m_nLabel1 = m_pNeighbor1->m_nLabel;
	if (m_pNeighbor2)
		m_nLabel2 = m_pNeighbor2->m_nLabel;
}

//
// Destructor for edges
// 
TEdge::~TEdge()
{
}

//
// Base class for graph nodes
/////////////////////////////////////////////////////////////////////////////

//
// Constructor for node
// 
TVertex::TVertex(TGraph* pGraph)
{
	Init(pGraph, -1);
}

//
// Constructor for node
// 
TVertex::TVertex(TGraph* pGraph, int nLabel)
{
	Init(pGraph, nLabel);
}

//
// Initialization
// 
void TVertex::Init(TGraph* pGraph, int nLabel)
{
	m_pGraph = pGraph;
	m_nLabel = nLabel;

	m_nNeighborCnt = 0;
	m_nMaxNeighborCnt = MAX_NEIGHBORS;
	m_pEdges = new TEdge*[m_nMaxNeighborCnt];

	m_pPrev = 0;
	m_pNext = 0;
}

//
// Destructor for node
// 
TVertex::~TVertex()
{
	if (m_pEdges) delete[] m_pEdges; //QIN_MOD
}


//
// Connect the specified edge to current vertex. 
// 
void TVertex::Connect(TEdge* pEdge)
{
	// Check if there is still storage space
	if (m_nNeighborCnt >= m_nMaxNeighborCnt)
	{
		// Need to expand the array
		m_nMaxNeighborCnt = (int) (m_nMaxNeighborCnt * 1.5);
		TEdge** pNewEdges = new TEdge*[m_nMaxNeighborCnt];
		memcpy(pNewEdges, m_pEdges, m_nNeighborCnt * sizeof(TEdge*));  // copy to the expanded storage space
		delete[] m_pEdges;
		m_pEdges = pNewEdges;
	}

	// Connect
	m_pEdges[m_nNeighborCnt ++] = pEdge;
}


//
// Disconnect the specified edge from current vertex.
//
void TVertex::Disconnect(TEdge* pEdge)
{
	for (int i = 0; i < m_nNeighborCnt; i ++)
	{
		if (m_pEdges[i] == pEdge)
		{
			m_nNeighborCnt --;
			memcpy(&m_pEdges[i], &m_pEdges[i+1], (m_nNeighborCnt - i) * sizeof(TEdge*));
			break;
		}
	}
}

// 
// Get the neighboring vertex on the specified edge
// 
// Return
//  0: if invalid link; otherwise the obtained neighbor
// 
TVertex* TVertex::GetNeighborVertex(TEdge* pEdge)
{
	if (pEdge->m_pNeighbor1 == this)
		return pEdge->m_pNeighbor2;
	else if (pEdge->m_pNeighbor2 == this)
		return pEdge->m_pNeighbor1;
	else
		return 0;
}

// 
// Get the edge between current and the specified neighboring
// vertex
// 
// Return
//  0: if the specified vertex not a neighbor; otherwise the obtained edge 
// 
TEdge* TVertex::GetEdge(TVertex* pVertex)
{
	for (int i = 0; i < m_nNeighborCnt; i ++)
	{
		if (GetNeighborVertex(m_pEdges[i]) == pVertex)
			return m_pEdges[i];
	}
	return 0;
}

//
// Base class for Graph
/////////////////////////////////////////////////////////////////////////////

//
// Constructor for graph
// 
TGraph::TGraph(bool bDirected)
{
	m_bDirected = bDirected;
	m_nVertexCnt = 0;
	m_nEdgeCnt = 0;
	m_nMaxDegree = 0;
	Clear();
}


//
// Destructor
// 
TGraph::~TGraph()
{
	Clear();
}


//
// Clear
// 
void TGraph::Clear()
{
	// Delete all vertice in the linked list
	TVertex* pVertex;
	while (m_nVertexCnt > 0)
	{
		pVertex = m_pFirstVertex;
		m_pFirstVertex = m_pFirstVertex->m_pNext;
		delete pVertex;
		m_nVertexCnt --;
	}
	m_pFirstVertex = 0;

	// Delete all edges in the linked list
	TEdge* pEdge;
	while (m_nEdgeCnt > 0)
	{
		pEdge = m_pFirstEdge;
		m_pFirstEdge = m_pFirstEdge->m_pNext;
		delete pEdge;
		m_nEdgeCnt --;
	}
	m_pFirstEdge = 0;
	m_pLastEdge = 0;
}


//
// Add a vertex to graph and the corresponding edges between the vertex
// and its neighbors. For constructing edges of neighboring pair, in
// directed graph, it adds two edges for each direction, while in undirected 
// graph only one edge is added. It should be noted that in directed graph
// there are always two edges for each neighboring pair, even if you do not
// wish to place a link for one of the direction. Brokeness on that direction
// can be implemented by setting some flags on the corresponding edges
// 
// Param in:
//  pVertex: pointer to the vertex to be added
// 
void TGraph::AddVertex(TVertex* pVertex)
{
	// Insert the specified vertex into the link list at the beginning position
	pVertex->m_pGraph = this;
	pVertex->m_pNext = m_pFirstVertex;
	if (m_pFirstVertex)
		m_pFirstVertex->m_pPrev = pVertex;
	m_pFirstVertex = pVertex;
	m_nVertexCnt ++;


	// Add corresponding edges
	TVertex* pNeighbor = FirstNeighbor(pVertex);  // get the first neighboring vertex that should be connected to
	if (pNeighbor)
	{
		do {

			AddEdge(pNeighbor, pVertex);
			if (m_bDirected)
				AddEdge(pVertex, pNeighbor);
			else
				pVertex->Connect(m_pFirstEdge);  // the first edge in the linked list is the common edge just generated
		} while (pNeighbor = NextNeighbor(pVertex));  // get the next neighboring vertex that should be connected to
	}
}

//
// Add directed edge to the graph.
// 
void TGraph::AddEdge(TVertex* pVertex, TVertex* pNeighbor)
{
	TEdge* pEdge = ConstructEdge(pVertex, pNeighbor);
	pEdge->m_pGraph = this;
	pVertex->Connect(pEdge);

	// max 2021: this line needs to be commented or else the number of neighbours of each vertex is incorrect
	//pNeighbor->Connect(pEdge);

	// Insert the specified edge into the linked list at the begining position
	pEdge->m_pNext = m_pFirstEdge;
	if (m_pFirstEdge)
		m_pFirstEdge->m_pPrev = pEdge;
	else
		m_pLastEdge = pEdge;
	m_pFirstEdge = pEdge;
	m_nEdgeCnt ++;
}

//
// Delete a vertex from graph
// 
// Param in:
//  pVertex: pointer to the vertex to be removed
/////////////////////////////////////////////////////////////////////////////
void TGraph::DelVertex(TVertex* pVertex)
{
	//printf("Node[%d] is being deleted!\n", pVertex->m_nLabel);
	// Delete all connected edges
	TEdge* pEdge;
	TVertex* pNeighbor;
	while (pVertex->m_nNeighborCnt > 0)
	{
		pEdge = pVertex->m_pEdges[0];
		pNeighbor = pVertex->GetNeighborVertex(pEdge);
		if (m_bDirected)
			DelEdge(pNeighbor, pVertex);
		else
			pNeighbor->Disconnect(pEdge);  // remove connection at the neighbor side
		DelEdge(pVertex, pNeighbor);
	}

	// Delete the vertex from the linked list
	m_nVertexCnt --;
	if (pVertex->m_pPrev)
		pVertex->m_pPrev->m_pNext = pVertex->m_pNext;
	else
	{
		if (m_pFirstVertex == pVertex)  //QIN_MOD_POINTER
			m_pFirstVertex = pVertex->m_pNext;
	}
	if (pVertex->m_pNext)
		pVertex->m_pNext->m_pPrev = pVertex->m_pPrev;
	delete pVertex;
}


//
/// Delete directed edge between two vertices
// 
void TGraph::DelEdge(TVertex* pVertex, TVertex* pNeighbor)
{
	//printf("Edge [%d:%d] is being deleted!\n", pVertex->m_nLabel, pNeighbor->m_nLabel);


	TEdge* pEdge;
	for (int i = 0; i < pVertex->m_nNeighborCnt; i ++)
	{
		pEdge = pVertex->m_pEdges[i];

		if (pNeighbor == pEdge->m_pNeighbor1 ||
		    pNeighbor == pEdge->m_pNeighbor2)
		{
			// Disconnect
			pVertex->Disconnect(pEdge);

            // Same deal as line 288
			// pNeighbor->Disconnect(pEdge);

			// Delete from the linked list
			m_nEdgeCnt --;

			if (pEdge->m_pPrev)
				pEdge->m_pPrev->m_pNext = pEdge->m_pNext;
			else
			{
				if (pEdge == m_pFirstEdge)  //QIN_MOD
					m_pFirstEdge = pEdge->m_pNext;
			}
			
			if (pEdge->m_pNext)
				pEdge->m_pNext->m_pPrev = pEdge->m_pPrev;
			else
			{
				if (pEdge == m_pLastEdge)
				{
					m_pLastEdge = pEdge->m_pPrev;
				}
			}
			delete pEdge;
			return;
		}
	}
}


//
// Merge two vertices in the graph. We merge the second vertex into the first
// 
// Param in:
//  pVertex1, pVertex2: the two vertices to be merged
// 
void TGraph::Merge(double& dbT1, double& dbT2, TVertex* pVertex1, TVertex* pVertex2)
{
	time_t start, end;

	time(&start);

	// Merge the vertex data
	pVertex1->MergeWith(pVertex2);

	// Connect to all neighbors of the second vertex being merged. If some of them are
	// already the neighbors of the first vertex, merge the corresponding links
	TEdge *pEdge1;  // edge of the first vertex
	TEdge *pEdge2;  // edge of the second vertex
	TVertex *pNeighbor;
	while (pVertex2->m_nNeighborCnt > 0)
	{
		// Find a neighbor
		pEdge2 = pVertex2->m_pEdges[0];
		pNeighbor = pVertex2->GetNeighborVertex(pEdge2);
		if (pNeighbor == pVertex1)  // the other vertex being merged
		{
			if (m_bDirected)
				DelEdge(pNeighbor, pVertex2);
			else
				pNeighbor->Disconnect(pEdge2);

			DelEdge(pVertex2, pNeighbor);
		}
		else
		{
			// Check if it is already a neighbor of the first vertex
			pEdge1 = pVertex1->GetEdge(pNeighbor);
			if (pEdge1)
			{
				// Already a neighbor of the first vertex, merge with current
                pEdge1->MergeWith(pEdge2);
				if (m_bDirected)
				{
					pEdge1 = pNeighbor->GetEdge(pVertex1);
					pEdge2 = pNeighbor->GetEdge(pVertex2);
					pEdge1->MergeWith(pEdge2);
					DelEdge(pNeighbor, pVertex2);
				}
				else 
					pNeighbor->Disconnect(pEdge2);
				DelEdge(pVertex2, pNeighbor);
			}
			else
			{
				// Not a neighbor yet, add it
				pVertex2->Disconnect(pEdge2);
				pVertex1->Connect(pEdge2);
				pEdge2->SetVertexNeighbor(pVertex1, pNeighbor);

				if (m_bDirected)
				{
					pEdge2 = pNeighbor->GetEdge(pVertex2);
					pEdge2->SetVertexNeighbor(pNeighbor, pVertex1);
				}
			}
		}

	}
	// Delete the nodes
	DelVertex(pVertex2);
	
	time(&end);
	dbT1 += difftime(end, start);

	time(&start);
	// Update the parameters associated with the merged vertex. Although some operations have
	// been completed in MergeWith(), there are some other operations that need
	// to be performed after the merging is over (e.g. shape feature)
	Update(pVertex1);
	time(&end);
	dbT2 += difftime(end, start);
}

//
// Merge two vertices in the graph. We merge the second vertex into the first
// 
// Param in:
//  pVertex1, pVertex2: the two vertices to be merged
// 
void TGraph::Merge(TVertex* pVertex1, TVertex* pVertex2)
{
	// Merge the vertex data
	pVertex1->MergeWith(pVertex2);

	// Connect to all neighbors of the second vertex being merged. If some of them are
	// already the neighbors of the first vertex, merge the corresponding links
	TEdge *pEdge1;  // edge of the first vertex
	TEdge *pEdge2;  // edge of the second vertex
	TVertex *pNeighbor;
	while (pVertex2->m_nNeighborCnt > 0)
	{
		// Find a neighbor
		pEdge2 = pVertex2->m_pEdges[0];
		pNeighbor = pVertex2->GetNeighborVertex(pEdge2);
		if (pNeighbor == pVertex1)  // the other vertex being merged
		{
			if (m_bDirected)
				DelEdge(pNeighbor, pVertex2);
			else
				pNeighbor->Disconnect(pEdge2);

			DelEdge(pVertex2, pNeighbor);
		}
		else
		{
			// Check if it is already a neighbor of the first vertex
			pEdge1 = pVertex1->GetEdge(pNeighbor);
			if (pEdge1)
			{
				// Already a neighbor of the first vertex, merge with current
				pEdge1->MergeWith(pEdge2);
//				pEdge1->MergeWith((TStrengthPenaltyBoundary*)pEdge2);
				if (m_bDirected)
				{
					pEdge1 = pNeighbor->GetEdge(pVertex1);
					pEdge2 = pNeighbor->GetEdge(pVertex2);
					pEdge1->MergeWith(pEdge2);
					DelEdge(pNeighbor, pVertex2);
				}
				else 
					pNeighbor->Disconnect(pEdge2);
				DelEdge(pVertex2, pNeighbor);
			}
			else
			{
				// Not a neighbor yet, add it
				pVertex2->Disconnect(pEdge2);
				pVertex1->Connect(pEdge2);
				pEdge2->SetVertexNeighbor(pVertex1, pNeighbor);

				if (m_bDirected)
				{
					pEdge2 = pNeighbor->GetEdge(pVertex2);
					pEdge2->SetVertexNeighbor(pNeighbor, pVertex1);
				}
			}
		}

	}
	DelVertex(pVertex2);
}

// If you want the node degrees to be correct, you need to call this method after graph construction.
void TGraph::UpdateMaxDegree() {
    int nMaxDegree = 0;
    TVertex* pCurVertex = this->m_pFirstVertex;
    while(pCurVertex){
        int nNeighb = pCurVertex->m_nNeighborCnt;
        if (nNeighb > nMaxDegree) nMaxDegree = nNeighb;
        pCurVertex = pCurVertex->m_pNext;
    }
    m_nMaxDegree = nMaxDegree;
}
