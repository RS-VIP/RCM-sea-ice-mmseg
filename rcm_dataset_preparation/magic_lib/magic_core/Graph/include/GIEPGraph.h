/** \file
 *
 * \brief TGIEPGraph class.
 *
 * \author Qiyao Yu
 *
 * \date 2006
 *
 */

// GIEPGraph.h : interface of the TGIEPGraph class
//
// Region adjancency graph (RAG) for graduated increased edge penalty
//
// Qiyao Yu, 2006
/////////////////////////////////////////////////////////////////////////////


#ifndef TGIEPGraph_Class
#define TGIEPGraph_Class

#include "RegionAdjacencyGraph.h"
#include "RegionClassifier.h"

/// Region adjancency graph (RAG) for graduated increased edge penalty (GIEP)
/**
 * Both Gaussian and Polarimetric Statistics are supported; this class
 * generates the proper TStrengthPenaltyBoundary subclass.
 */
class TGIEPGraph : public TRegionAdjacencyGraph
{
public:
	TGIEPGraph(TRegionClassifier* pClassifier,
	           TGrayImage<int>* pMap, TMonoImage* pMask, TGrayImage<FLOAT>** ppSrcImg,
	           int nFeatureCnt, TRegionAdjacencyGraph::statistics_mode enStatMode, TGrayImage<FLOAT>* pGrad,
               TRegionAdjacencyGraph::boundary_mode enBoundMode=FULL);

//	TGIEPGraph(TRegionClassifier* pClassifier,
//	          TMonoImage* pMask,  TGrayImage<FLOAT>** ppSrcImg, int nFeatureCnt); //QIN_MOD

// Attributes
protected:
	/// classfier for the regions
	TRegionClassifier* m_pClassifier;
	/// flag indicating if no classifier is specified
	BOOL m_bNullClassifier;

	/// Gradient image (for edge strength).
	TGrayImage<FLOAT>* m_pGradient;

	/// feature model weight
	double m_dbAlpha;
	/// spatial context + edge penalty weight
	double m_dbBeta;

	/// edge strength histogram for computing the penalty function parameter
	double m_pEdgeHist[256];
	/// total number of edge sites
	int m_nEdgeSiteTotal;  
	/// maximum of edge strength
	double m_dbEdgeStrengthMax;

	/// penalty function parameter	
	double m_dbK;

public:
	/// number of iterations
	int m_nIter;
	/// update mode of Beta on each iteration.
	enum enBetaMode { FIXED = 0x00, MLFIXED = 0x01, MLFISH = 0x02 };

// Operations
public:
	/// classification
	TImage* Classify(TRegionClassifier*& pRetClassifier, int nBetaMode, bool bOutIntermResult, int nIterations, double dbK, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2);

	/// update the penalty for binary relation of regions
	void UpdatePenalty(); 

	TRegionClassifier* GetClassifier() { return m_pClassifier; };	
	TGrayImage<FLOAT>* GetGradient() { return m_pGradient; }; 
	double& GetCurBeta() { return m_dbBeta; };
	double& GetCurAlpha() { return m_dbAlpha; };
	double& GetCurK() { return m_dbK; };

	void UpdateBeta(double dbBeta1, int mode);
	void UpdateBetaFisher(double dbK, double dbBeta1, double dbBeta2, double dbFisher, int it);

    void SetClassifier(TRegionClassifier* pClassifier);
    /// estimate the boundary pixel population over the overall
    double GetBoundaryPopulationRatio();
    /// get weight for edge penalty
    double UpdateBeta(double dbEdgeRatio);
    /// update the penalty function parameter
    void UpdateK(double dbK, int nIter);

private:
	/// initialization
	void Init(TRegionClassifier* pClassifier); 

protected:
	/// generate a struct for the specified boundary
	/**
	 * actual subclass of TBoundary created depends on
	 * TRegionAdjacencyGraph::m_enStatisticsMode of this instance)
	 */
	virtual TBoundary* CreateNewBoundary(int nLabel1, int nLabel2, TChainCode* pChain, BOOL bOnGrid);

	/// compute the maximum of edge strength
	double ComputeMaxEdgeStrength();
	/// compute the edge strength histogram
	void ComputeEdgeStrengthHist();



	/// estimate the boundary pixel population over the overall	
	double GetBoundaryPopulationRatio(TGrayImage<int>* pMap);
	/// compute the number of pairs of neighboring sites belonging to different classes
	int ComputeEdgePairCnt(TGrayImage<int>* pMap);  

	/// get weight for feature model
	double UpdateAlpha(int nIter, double dbAlpha1, double dbAlpha2);

public:

	void SetClassLabels(std::vector<int> vClassLabels);



// Implementation
public:
	virtual ~TGIEPGraph();

};

#endif

/////////////////////////////////////////////////////////////////////////////
