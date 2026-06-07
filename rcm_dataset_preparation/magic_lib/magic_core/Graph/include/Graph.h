/** \file
 *
 * \brief TGraph class.
 *
 */

#ifndef TGraph_Class
#define TGraph_Class

#define INIT_MAXDIST 100000000  // 10^8 a large enough number

class TVertex;
class TGraph;

/// Class for link between graph nodes
class TEdge
{
public:
	TEdge(TGraph* pGraph);
	TEdge(TVertex* pNeighbor1, TVertex* pNeighbor2);

public:
	/// graph this edge belongs to
	TGraph* m_pGraph;

	TVertex* m_pNeighbor1;
	TVertex* m_pNeighbor2;
	int m_nLabel1, m_nLabel2;  // the two labels of the vertice linked by the edge

	int m_nIndex; // edge index, only sometimes used, you can probably ignore it

	// For constructing a doubly linked list of all edges in graph
	TEdge* m_pPrev;
	TEdge* m_pNext;

private:
	void Init(TVertex* pNeighbor1, TVertex* pNeighbor2);  // initialization

public:
	void SetVertexNeighbor(TVertex* pNeighbor1, TVertex* pNeighbor2);

	virtual double GetEdgeWeight() { return 0; }
	virtual void MergeWith(TEdge* pEdge) { }
public:
	virtual ~TEdge();
};


/// Class for graph nodes
class TVertex
{
public:
	TVertex(TGraph* pGraph);
	TVertex(TGraph* pGraph, int nLabel);

public:
	/// graph this vertex belongs to
	TGraph* m_pGraph;

	/// region identifier
	int m_nLabel;

	/// pointer to neighbor links
	TEdge** m_pEdges;

	/// number of neighbors
	int m_nNeighborCnt;

	/// maximum number of neighbors, for allocating memory
	int m_nMaxNeighborCnt;

	/// For constructing a doubly linked list of all vertice in graph
	TVertex* m_pPrev;
	/// For constructing a doubly linked list of all vertice in graph
	TVertex* m_pNext;

private:
	/// initialization
	void Init(TGraph* pGraph, int nLabel);  

public:
	void Connect(TEdge* pEdge);
	void Disconnect(TEdge* pEdge);

	TVertex* GetNeighborVertex(TEdge* pEdge);
	TEdge* GetEdge(TVertex* pVertex);

	virtual double GetVertexWeight() { return 0; };
	virtual void MergeWith(TVertex* pNeighbor) {};

public:
	virtual ~TVertex();
};

/// Base class for graph representation
class TGraph
{
public:
	TGraph(bool bDirected);

// Attributes
protected:
	/// flag indicating if the graph is directed or not
	bool m_bDirected; 

	/// the first graph vertex in the link list
	TVertex* m_pFirstVertex;
	
	/// number of vertice
	int m_nVertexCnt;

	/// the first graph edge in the link list
	TEdge* m_pFirstEdge;

	/// the last graph edge in the link list
	TEdge* m_pLastEdge;

	/// number of edges
	int m_nEdgeCnt;

	/// maximum vertex degree
	int m_nMaxDegree;

// Operations
public:
	void AddVertex(TVertex* pVertex);
	void AddEdge(TVertex* pVertex, TVertex* pNeighbor);
	void DelVertex(TVertex* pVertex);
	void DelEdge(TVertex* pVertex, TVertex* pNeighbor);
	void Merge(TVertex* pVertex1, TVertex* pVertex2);
	void Merge(double& dbT1, double& dbT2, TVertex* pVertex1, TVertex* pVertex2);
	void UpdateMaxDegree();

	TVertex* GetVertices() { return m_pFirstVertex; };
	TEdge* GetEdges() { return m_pFirstEdge; };

	int GetVertexCnt() { return m_nVertexCnt; };
	int GetEdgeCnt() { return m_nEdgeCnt; };
	int GetMaxDegree() { return m_nMaxDegree; };

protected:
	void Clear();

	/// Find the first neighbor for adding the current vertex
	virtual TVertex* FirstNeighbor(TVertex* pVertex){return (TVertex*)0;}

	// Find the next neighbor for adding the current vertex
	virtual TVertex* NextNeighbor(TVertex* pVertex){return (TVertex*)0;}

	/// Construct an edge for adding the current vertex 
	virtual TEdge* ConstructEdge(TVertex* pNeighbor1, TVertex* pNeighbor2)
	{ return (new TEdge(pNeighbor1, pNeighbor2)); }

	/// Update the parameters associated with the specified vertex
	virtual void Update(TVertex* pVertex) {};

// Implementation
public:
	virtual ~TGraph();

};

#endif
