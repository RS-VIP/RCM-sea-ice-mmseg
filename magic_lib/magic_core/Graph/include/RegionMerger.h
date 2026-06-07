/** \file
 *
 * \brief TRegionMerger class.
 *
 */

//#pragma warning( disable : 4250 )

#ifndef TRegionMerger_Class
#define TRegionMerger_Class

#include "GrayImage.h"
#include "Region.h"
#include "GIEPGraph.h"
#include<tuple>


/// A GIEPGraph that merges similar regions based on IRGS energy function.
class TRegionMerger : public TGIEPGraph
{
public:
	/// Region based constructor - nodes are watershed regions.
	TRegionMerger(TRegionClassifier* pClassifier,
	              TGrayImage<int>* pMap, TMonoImage* pMask, 
				  TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt, TRegionAdjacencyGraph::statistics_mode enStatMode,
				  TGrayImage<FLOAT>* pGrad,
                  TRegionAdjacencyGraph::boundary_mode enBoundMode=FULL);

//	/// Pixel based constructor - nodes are pixels.
//	TRegionMerger(TRegionClassifier* pClassifier,
//				  TMonoImage* pMask,
//	              TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt);

// Attributes
public:
	/// Update mode of Beta parameter.
	enum enBetaMode { FIXED = 0x00, MLFIXED = 0x01, MLFISH = 0x02 };

protected:

	/// Which edge connects the two most similar regions.
	TEdge* m_pNearestEdge;

// Operations
public:
	//TImage* Merge(double& dbInitTime, double& dbFinalBeta, TRegionClassifier*& pRetClassifier, int nBetaMode, bool bOutIntermResult, int nIterations, double dbK, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2);  // region merge
	//TImage* Merge(int nTmp, double& dbFinalBeta, TRegionClassifier*& pRetClassifier, int nBetaMode, bool bOutIntermResult, int nIterations, double dbK, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2);  // region merge
	
	/// Region merge - drives the entire region labeling and merging iterations.
	TImage* Merge(double& dbFinalBeta, TRegionClassifier*& pRetClassifier, int nBetaMode, bool bOutIntermResult, int nIterations, double dbK, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2);  // region merge

    std::tuple<double, double, bool> MergingStep(double dbK, double dbBeta1, double dbBeta2, int nIt);

	/// Sort the edge by its weight (merging criterion) //QIN_MOD_20090810
	void SortEdge();

	/// get a new graph with each vertex representing a connected region of the same class
	TGraph* GetConnectedClassRegions();

	/// get the map in which each label (color) represents a connected region of the same class
	TGrayImage<int>* GetConnectedClassRegionMap();

	/// Save Intermediate Result - completely commented out for some reason.
	void SaveIntermediateResult(int iIter, int iVrtxCnt );

	/// merge nearest vertices
	BOOL Merge();

    /// merge nearest vertices and return their indices
    std::vector<std::pair<int,int>> MergeAndReturnIndex();

protected:

	/// update the parameters associated with the specified vertex
	virtual void Update(TVertex* pVertex);

	/// find the edge that corresponds to the two nearest vertice
	TEdge* FindNearest(); 
	
	/// determining if merging iteration should stop	
	BOOL ShouldStop();

	/// Set the label mapping for the specified vertex
	void SetLabelMapping(int* nLabels, TRegion* pCurrent, int nNewLabel);


// Implementation
public:
	virtual ~TRegionMerger();

};

#endif

/////////////////////////////////////////////////////////////////////////////
