/** \file
 *
 * \brief TRegionAdjacencyGraph class.
 *
 */
#ifndef TRegionAdjacencyGraph_Class
#define TRegionAdjacencyGraph_Class

#include <vector>
#include "Graph.h"
#include "MonoImage.h"
#include "GrayImage.h"
#include "Boundary.h"
#include "Region.h"
#include "GaussianRegion.h"


/// Region Adjacency Graph
/**
 *
 * This class is the base RAG class. This can be constructed to handle
 * polarimetric or gaussian statistics. The initial nodes that it starts off
 * with can be individual pixels (very computationally expensive and not used)
 * or watershed regions.
 *
 */
class TRegionAdjacencyGraph : public TGraph
{

// Attributes
public:

    /** Construction mode.
    * PIXEL: each graph node is a pixel in the original image; boundary is
    *         defined as located between pixels.
    *  REGION: each graph node is a region obtained from a preliminary segmentation
    *          result; boundary is defined on pixel sites.
    */
    enum construction_mode {PIXEL = 0,  // construction of graph is pixel based
        REGION = 1};  // construction of graph is region based

	/// statistics mode: Gaussian or Polarimetric Wishart)
	// (Polarimetric Wishart not implemented after a major code overhaul but should be straighforward to bring back)
	enum statistics_mode {  GAUSSIAN = 0,               // Gaussian statistics.
	                        POLARIMETRIC_WISHART = 1};  // Polarimetric Wishart

    enum boundary_mode {FULL = 0, NO_INTERSECTIONS=1};

    /// boundary segment list in the form of AVL tree, for fast neighbour lookup.
    TAVLTree m_BoundaryList;


public:

	/// Constructor for region mode; requires a watershed map.
	TRegionAdjacencyGraph(TGrayImage<int>* pMap, TMonoImage* pMask,
	                      TGrayImage<FLOAT>** ppSrcImg,
                          int nFeatureCnt,
                          TRegionAdjacencyGraph::statistics_mode enStatMode = GAUSSIAN,
                          TRegionAdjacencyGraph::boundary_mode enBoundaryMode = FULL);


protected:

    /// construction mode: pixel or region
    construction_mode m_enConstructMode;

	/// statistics mode: gaussian or polarimetric
	statistics_mode m_enStatisticsMode;

    /// boundary mode: full (with intersection points) or without boundary intersection points
    boundary_mode m_enBoundaryMode;
	
	/// the map with labels corresponding to id of nodes in graph
	TGrayImage<int>* m_pMap;

	/// the mask of pixels to use and to ignore.
	TMonoImage* m_pMask;

	/// the original feature image(s)
	TGrayImage<FLOAT>** m_ppSrcImg;

	/// number of feature images
	int m_nFeatureCnt;

    /// vertex lookup table
    TVertex** m_pVertexLUT;

    /// edge lookup table
    TEdge** m_pEdgeLUT;

// Operations
public:

	/// construct the region adjacency graph
	virtual void ConstructGraph();

	TGrayImage<int> * GetMap() const { return m_pMap; };
	TMonoImage* GetMask() { return m_pMask; };
	TGrayImage<FLOAT>** GetImageSource() { return m_ppSrcImg; };
	int GetFeatureCnt() { return m_nFeatureCnt; };
	const statistics_mode GetStatMode() { return m_enStatisticsMode; };
    const boundary_mode GetBoundaryMode() { return m_enBoundaryMode; };


    // These are dangerous, try not to use them since m_pVertexLUT and m_pEdgeLUT might change
    TRegion* GetRegion(int nLabel){ return (TRegion*) m_pVertexLUT[nLabel];}
    TEdge* GetEdge(int nLabel) {return m_pEdgeLUT[nLabel];}


protected:

	/// Creates correct subclass of TVertex based on the statistics mode (TGaussianRegion or TPolarimetricRegion)
	virtual TVertex* CreateNewRegion(int nLabel);

	/// generate a struct for the specified boundary
	virtual TBoundary* CreateNewBoundary(int nLabel1, int nLabel2, TChainCode* pChain, BOOL bOnGrid) { return new TBoundary(this, nLabel1, nLabel2, pChain, bOnGrid); };

	/// get the first neighbor for constructing graph
	virtual TVertex* FirstNeighbor(TVertex* pVertex);

	/// get the next neighbor for constructing graph
	virtual TVertex* NextNeighbor(TVertex* pVertex);

	/// construct edge between vertex neighbors
	virtual TEdge* ConstructEdge(TVertex* pNeighbor1, TVertex* pNeighbor2);

	/// construct the graph in REGION mode
	virtual void ConstructGraphRegion();

	/// load the boundaries into a list
	virtual void LoadBoundarySegments();

    /// load the boundaries, excluding the boundary intersection points
    virtual void LoadInnerBoundarySegments();

	/// add a boundary segment into the list
	virtual void AddBoundarySegment(TGrayImage<BYTE>* pVisitRecord, POINTex StartPt);

	/// compute the corresponding chain code
	virtual TChainCode* ComputeChainCode(POINTex StartPt, int nLabel1, int nLabel2,
	                             TGrayImage<BYTE>* pVisitRecord);

	/// get labels of the two neighboring segments on the specifed site
	virtual BOOL GetNeighboringRegionLabels(int& nLabel1, int& nLabel2,
	                                POINTex Pt);

	/// check if boundary segment with the same labels already exists in the list
	TBoundary* IsBoundarySegmentInList(int nLabel1, int nLabel2);

	/// check if the site is boundary of segments specified by the two labels
	virtual BOOL IsBoundaryOfRegions(int nLabel1, int nLabel2, POINTex Pt);

    /// check if the site is boundary of segments, excluding the boundary intersection points
    virtual BOOL IsInnerBoundaryOfRegions(int nLabel1, int nLabel2, POINTex Pt);

	/// Check if is is on the border (no matter grid or mask)
	BOOL IsBorderSegment(int nRow, int nCol);

    /// Set the label mapping for the specified vertex
    void SetLabelMapping(int* nLabels, TRegion* pCurrent, int nNewLabel);

// Implementation
public:

    /// get the map in which each label (color) represents a connected region of the same class
//    TGrayImage<int>* GetConnectedClassRegionMap();

    BOOL SetClassLabelByIndex(int nClassLabel, int nIndex);

	virtual ~TRegionAdjacencyGraph();

};

#endif

/////////////////////////////////////////////////////////////////////////////
