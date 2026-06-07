/** \file
 *
 * \brief TRegionClassifier class
 *
 */

#ifndef TRegionClassifier_Class
#define TRegionClassifier_Class

//#include "ProcObject.h"
#include "RegionAdjacencyGraph.h"
#include "Region.h"
#include "Vector.h"
#include "SegmentationStatistics.h"

/// Base class for classifying the regions
class TRegionClassifier
{
public:
	TRegionClassifier();

// Attributes
protected:
	/// The region adjacency graph
	TRegionAdjacencyGraph* m_pGraph;

// Operations
public:
	/// Set the graph to the specified one.
	void SetGraph(TRegionAdjacencyGraph* pGraph) { m_pGraph = pGraph; };

	/// initialize the classifier
	virtual bool Init() { return true; }; 

	/// intermediate classification
	virtual bool IntermediateClassify(); 

	/// final classification
	virtual bool FinalClassify(); 

	/// update class stats based on current assignment of labels
	virtual bool UpdateParam() { return true; };

	/// final update class stats based on current assignment of labels
	virtual bool FinalUpdateParam() {return true; };

	/// return overall energy
	virtual double OverallEnergy();

	/// compute the minimum Fisher criterion
	virtual double ComputeMinFisher() { return 1; };

	/// class label of specified region.
	virtual int GetClassLabel(TRegion* pRegion) { return pRegion->m_nClassLabel; };

	/// Gets the classified map.
	virtual TImage* GetClassMap() { return 0; };


	/// Get segmentation statistics object that lists all the class statistics.
	virtual TSegmentationStatistics* GetSegStats(TGrayImage<int>* pMap)
	{
		return 0;
	}

protected:
	/// classifies vecfeature vector into one of the two class labels
	virtual int BoundaryPixelClassify(TVector* VecFeature, int nLabel1, int nLabel2) { return 0; }; //QIN_ADD

// Implementation
public:
	virtual ~TRegionClassifier();

};

#endif

