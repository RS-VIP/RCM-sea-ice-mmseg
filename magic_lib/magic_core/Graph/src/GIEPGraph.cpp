// GIEPGraph.cpp : implementation of the TGIEPGraph class
//
// Region adjancency graph (RAG) for graduated increased edge penalty
//
// Qiyao Yu, 2006
/////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "GIEPGraph.h"
#include "EstimateUtil.h"
//#include "PolarimetricStrengthPenaltyBoundary.h"
#include "GaussianStrengthPenaltyBoundary.h"

#define SCHEDULE_STEP 0.02

#define BETA_CNT 23
const static double betaLUT[BETA_CNT][2] = {
	{0.6410440, 0.05},
	{0.5941100, 0.15},
	{0.5617630, 0.20}, 
	{0.5195670, 0.25}, 
	{0.4613010, 0.30}, 
	{0.3728720, 0.35}, 
	{0.2246080, 0.40}, 
	{0.1375060, 0.45}, 
	{0.0928173, 0.50}, 
	{0.0663682, 0.55}, 
	{0.0503680, 0.60}, 
	{0.0397756, 0.65}, 
	{0.0327342, 0.70}, 
	{0.0280129, 0.75}, 
	{0.0248621, 0.80}, 
	{0.0228758, 0.85}, 
	{0.0212809, 0.90}, 
	{0.0188605, 1.00}, 
	{0.0175266, 1.15}, 
	{0.0163295, 1.40}, 
	{0.0157941, 1.70}, 
	{0.0153300, 2.0}, 
	{0.0147133, 3.0}, 
//	{0.0135030, 5.0},
};

//
// Constructor (region based)
// 
TGIEPGraph::TGIEPGraph(TRegionClassifier* pClassifier,
                       TGrayImage<int>* pMap,
					   TMonoImage* pMask, 
                       TGrayImage<FLOAT>** ppSrcImg,
					   int nFeatureCnt, TRegionAdjacencyGraph::statistics_mode enStatMode,
                       TGrayImage<FLOAT>* pGrad, TRegionAdjacencyGraph::boundary_mode enBoundMode)
	: TRegionAdjacencyGraph(pMap, pMask, ppSrcImg, nFeatureCnt, enStatMode, enBoundMode)
{
	/*Init(pClassifier, &pGrad); //QIN_MOD*/
	Init(pClassifier);
	m_pGradient = pGrad;
}

////
//// Constructor (pixel based)
////
//TGIEPGraph::TGIEPGraph(TRegionClassifier* pClassifier,
//					   TMonoImage* pMask,
//                       TGrayImage<FLOAT>** ppSrcImg,
//					   int nFeatureCnt)
//	: TRegionAdjacencyGraph(ppSrcImg[0]->Width(), ppSrcImg[0]->Height(), pMask, ppSrcImg, nFeatureCnt) //QIN_MOD, assuming each image component has the same size
//{
//	/*Init(pClassifier, ppSrcImg); //QIN_MOD*/
//	Init(pClassifier);
//	m_pGradient = NULL;
//}

//
// Initialization
// 
void TGIEPGraph::Init(TRegionClassifier* pClassifier) //QIN_MOD
{
	m_pClassifier = pClassifier;
	if (!m_pClassifier)
	{
		m_bNullClassifier = true;
		m_pClassifier = new TRegionClassifier;
	}
	else
		m_bNullClassifier = false;

	m_pClassifier->SetGraph(this);
	/*m_pEdgeStrength = pEdgeStrength; //QIN_MOD*/

	m_nEdgeSiteTotal = -1;
}

//
// Destructor
// 
TGIEPGraph::~TGIEPGraph()
{
	if (m_bNullClassifier)
		delete m_pClassifier;  // a classifier pre-allocated, need to be deleted
}

//
// Estimate the boundary pixel population over the overall
// 
double TGIEPGraph::GetBoundaryPopulationRatio()
{
	// Compute the number of edge pixels
	TBoundary* pEdges = (TBoundary*) m_pFirstEdge;
	int nEdgeCnt = 0;
	while (pEdges)
	{
		if (m_bNullClassifier ||  // no specified classifier, each region is a class
		    (((TRegion*) pEdges->m_pNeighbor1)->m_nClassLabel !=
		     ((TRegion*) pEdges->m_pNeighbor2)->m_nClassLabel))
			nEdgeCnt += pEdges->m_nLength;
		pEdges = (TBoundary*) pEdges->m_pNext;
	}

	// Compute the total number of pixels
	int nTotalCnt = 0;
	if (!m_pMask)
		nTotalCnt = m_pMap->Width() * m_pMap->Height();
	else
	{
		for (int i = m_pMap->Height() - 1; i >= 0; i --)
			for (int j = m_pMap->Width() - 1; j >= 0; j --)
		{
			if ((*m_pMask)(i, j))
				nTotalCnt ++;
		}
	}

	double dbRatio = ((double) nEdgeCnt) / nTotalCnt;
	//if (m_enConstructMode == PIXEL)  // if boundary is not defined on the original grid, the actual ratio should be around 1/2 (due to scan order)
	//	dbRatio = dbRatio / 2; //QIN_MOD

	return dbRatio;
}

//
// Estimate the boundary pixel population over the overall
// 
double TGIEPGraph::GetBoundaryPopulationRatio(TGrayImage<int>* pMap)  //Here pMap is class map not region label map
{
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();
	// Compute the total number of pixels
	int nTotalCnt = 0, i, j;
	if (!m_pMask)
		nTotalCnt = nWidth  * nHeight;
	else
	{
		for (i = nHeight - 1; i >= 0; i --)
			for (j = nWidth - 1; j >= 0; j --)
		{
			if ((*m_pMask)(i, j))
				nTotalCnt ++;
		}
	}

	double dbRatio;
	int r, c;
	int nLabel, l;
	//if (m_enConstructMode == REGION)
    if(1)
	{
		int nEdgeCnt = 0;
		for (i = pMap->Height() - 1; i >= 0; i --)
			for (j = pMap->Width() - 1; j >= 0; j --)
		{
			if (m_pMask && !(*m_pMask)(i, j))
				continue;
			if ((*pMap)(i, j) < 0)  // boundary pixel
			{
				nLabel = -1;
				for (int k = 0; k < 8; k ++)
				{
					r = i + neighbor8[k].Row;
					c = j + neighbor8[k].Col;
					if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
						continue;
					if (m_pMask && !(*m_pMask)(r, c))
						continue;

					if ((l = (*pMap)(r, c)) >= 0)  // class label not boundary label
					{
						if (nLabel == -1)  // the first class label found in neighborhood
							nLabel = l;
						else if (nLabel != l)  // having different classes in neighborhood
						{
							nEdgeCnt ++;
							break;
						}
					}
				}
			}
		}

		dbRatio = ((double) nEdgeCnt) / nTotalCnt;
	}
	else
		dbRatio = ComputeEdgePairCnt(pMap) / 6.0 / nTotalCnt;

	return dbRatio;
}

//
// Compute the number of pairs of neighboring sites belonging to different
// classes
// 
// Param in
//  pMap: the class map 
// 
int TGIEPGraph::ComputeEdgePairCnt(TGrayImage<int>* pMap)
{
	int** data = pMap->GetData();
	int nWidth = pMap->Width();
	int nHeight = pMap->Height();

	int r, c, nLabel, nEdgeCnt = 0;
	for (int nRow = 0; nRow < nHeight; nRow ++)
		for (int nCol = 0; nCol < nWidth; nCol ++)
	{
		if (m_pMask && !(*m_pMask)(nRow, nCol))
			continue;
		nLabel = data[nRow][nCol];
		for (int j = 0; j < 8; j ++)
		{
			r = nRow + neighbor8[j].Row;
			c = nCol + neighbor8[j].Col;
			if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
				continue;
			if (m_pMask && !(*m_pMask)(r, c))
				continue;
			if (data[r][c] != nLabel)
				nEdgeCnt ++;
		}
	}

	return nEdgeCnt;
}

//
// Get weight for edge penalty
// 
// Return
//  The weight
// 
double TGIEPGraph::UpdateBeta(double dbEdgeRatio)
{
	double dbBeta;
	if (dbEdgeRatio >= betaLUT[0][0])  // larger than all recorded, select the smallest beta
		dbBeta = betaLUT[0][1];
	else
	{
		int i;
		for (i = 1; i < BETA_CNT; i ++)
		{
			if (dbEdgeRatio > betaLUT[i][0])
				break;
		}

		if (i >= BETA_CNT)  // smaller than all recorded, select the largest beta
			dbBeta = betaLUT[BETA_CNT - 1][1];
		else
			dbBeta = betaLUT[i - 1][1] + 
			         (betaLUT[i][1] - betaLUT[i - 1][1]) /
			         (betaLUT[i][0] - betaLUT[i - 1][0]) *
			         (dbEdgeRatio - betaLUT[i - 1][0]);  // interpolate
	}

//	cout << "The ratio of edge pixels is " << dbEdgeRatio << endl;
//	cout << "The beta is " << dbBeta << endl;
	return dbBeta;
}

// update beta according to ML+Fisher criterion
void TGIEPGraph::UpdateBetaFisher(double dbK, double dbBeta1, double dbBeta2, double dbFisher, int it)
{
    UpdateK(dbK, it);
    m_dbBeta = UpdateBeta(GetBoundaryPopulationRatio());
    m_dbBeta *= ( dbBeta1 * dbFisher / dbBeta2 / (1 + dbFisher / dbBeta2) );
    UpdatePenalty();
    // cout << "Fisher: " << dbFisher << ", Adjusted Beta: " << m_dbBeta << endl;
}

void TGIEPGraph::UpdateBeta(double dbBeta1, int mode)
{
    double dbEdgeRatio = GetBoundaryPopulationRatio();
    double dbBeta;
    if (dbEdgeRatio >= betaLUT[0][0])  // larger than all recorded, select the smallest beta
        dbBeta = betaLUT[0][1];
    else
    {
        int i;
        for (i = 1; i < BETA_CNT; i ++)
        {
            if (dbEdgeRatio > betaLUT[i][0])
                break;
        }

        if (i >= BETA_CNT)  // smaller than all recorded, select the largest beta
            dbBeta = betaLUT[BETA_CNT - 1][1];
        else
            dbBeta = betaLUT[i - 1][1] +
                     (betaLUT[i][1] - betaLUT[i - 1][1]) /
                     (betaLUT[i][0] - betaLUT[i - 1][0]) *
                     (dbEdgeRatio - betaLUT[i - 1][0]);  // interpolate
    }

//	cout << "The ratio of edge pixels is " << dbEdgeRatio << endl;
//	cout << "The beta is " << dbBeta << endl;
    m_dbBeta = dbBeta*dbBeta1;
}


//
// Get weight for feature model
// 
// Return
//  The weight
// 
double TGIEPGraph::UpdateAlpha(int nIter, double dbAlpha1, double dbAlpha2)
{
	double dbAlpha;

	dbAlpha = dbAlpha1 * pow(dbAlpha2, nIter) + 1.0;
//	cout << "The Alpha is " << dbAlpha<< endl;
	return dbAlpha;
}


void TGIEPGraph::SetClassifier(TRegionClassifier *pClassifier)
{
    if (m_pClassifier) delete m_pClassifier;
    m_pClassifier = pClassifier;
    m_bNullClassifier = false;

    // this should probably belong in a different function but I am lazy so I am leaving it in here because it fits my
    // use case. But some later students may be confused about why a set operation takes so long lmao
    m_pClassifier->Init();
}

// 
// Classification
// 
// Param in
//  nIterations: number of iterations
// 
// return
//  the segmentation map
// 
TImage* TGIEPGraph::Classify(TRegionClassifier*& pRetClassifier, int nBetaMode, bool bOutIntermResult, int nIterations, double dbK, double dbBeta1, double dbBeta2, double dbAlpha1, double dbAlpha2)
{
	TImage *pRet;

	// Construct the graph
	ConstructGraph();

	// Output intermediate results for each iteration
	char szOut[64];
	szOut[0] = 'o';
	szOut[1] = 'u';
	szOut[2] = 't';
	szOut[6] = '.';
	szOut[7] = 'b';
	szOut[8] = 'm';
	szOut[9] = 'p';
	szOut[10] = 0;

	// Initialize the classifier
	if (!m_bNullClassifier)
	    std::cout<< "Initializing Classifier!!!" << std::endl;
		m_pClassifier->Init();

	if ( nBetaMode == FIXED )
	{
		m_dbBeta = dbBeta1;
		//cout << "Fixed Beta is: " << m_dbBeta << endl;
	}
	// used for updating gui /// Added by SHARIF
	//string smsg = "Iteration # ";
	//char cnum[111] = {0};
	int iTotalMergedRegion = 0;
	double dbFisher;

	for ( m_nIter = 0; m_nIter < nIterations; m_nIter++ )
	{
		// Update GUI to display iteration number
		//smsg = "Iteration # " + StringConvertFunctions::ConvertToString(m_nIter + 1, 0);
		//smsg += itoa(m_nIter + 1,cnum,10);
		//strcpy(&cnum[0], smsg.c_str());
		//if(CanShowStatus()) DisplayStatus(smsg,2);


		// Update the parameters
		UpdateK(dbK, m_nIter);
		cout << "The coeff of edge energy penalty function: " << m_dbK << endl;



		//M2021 -  commented to remove dependency on regionclustering
//		if ( ((TRegionClustering*) m_pClassifier)->GetType() & TRegionClustering::VAFM )
//			m_dbAlpha = UpdateAlpha(m_nIter, dbAlpha1, dbAlpha2);

		if ( nBetaMode == MLFIXED )
		{
			m_dbBeta = UpdateBeta(GetBoundaryPopulationRatio());
			m_dbBeta *= dbBeta1;
			cout << "Adjusted Beta is: " << m_dbBeta << endl;
		}

		if ( nBetaMode == MLFISH )
		{
			m_dbBeta = UpdateBeta(GetBoundaryPopulationRatio());
	 		dbFisher = m_pClassifier->ComputeMinFisher();
			m_dbBeta *= ( dbBeta1 * dbFisher / dbBeta2 / (1 + dbFisher / dbBeta2) );
			cout << "Fisher: " << dbFisher << ", Adjusted Beta: " << m_dbBeta << endl;
		}

		// Compute the energy corresponding to edge information
		UpdatePenalty();

		// Classification
		if (!m_bNullClassifier)
			m_pClassifier->IntermediateClassify();

		if ( !m_bNullClassifier )
		{
			if (!m_pClassifier->UpdateParam())
				goto END;
		}

		//cout << endl << "Iteration " << m_nIter << " finished" << endl;
		//cout << "Energy: " << m_pClassifier->OverallEnergy() << endl;

		if ( bOutIntermResult )
		{
			// Output intermediate results for each iteration
			int k = m_nIter, l;
			for ( int j = 2; j >= 0; j-- )  // 4 digits
			{
				l = k % 10;
				k = k / 10;
				szOut[3 + j] = '0' + l;
			}
			pRet = m_pClassifier->GetClassMap();
			TEstimateUtil* pEstUtil = new TEstimateUtil;
			TGrayImage<int>* pRet1 = new TGrayImage<int>((TGrayImage<int>*) pRet);
			pEstUtil->SortMap(pRet1, m_pMask, m_ppSrcImg, m_nFeatureCnt);
			TGrayImage<int>* pDispMap = pEstUtil->DispMap(pRet1, m_pMask, false);
//			pDispMap->Save(szOut);
			delete pEstUtil;
			delete pRet1;
			delete pDispMap;
		}

	}

END:
	// Final classification
	if (!m_bNullClassifier)
		m_pClassifier->FinalUpdateParam();

	pRetClassifier = m_pClassifier;

	if (!(pRet = m_pClassifier->GetClassMap()))
	//if (!(pRet = GetClassMap()))
	{
		TEstimateUtil* pEstUtil = new TEstimateUtil;
		pRet = new TGrayImage<int>(m_pMap);
		pEstUtil->RelabelMap( (TGrayImage<int>*) pRet, m_pMask );
		delete pEstUtil;
	}

	return pRet;
}

//// 
//// Update the penalty function parameter
//// 
//void TGIEPGraph::UpdateK(double dbK, int nIter)
//{
//	ofstream outfile;
//	if (nIter == 0)
//		outfile.open("f1.txt",ios::out);
//	else
//		outfile.open("f1.txt",ios::app);
//
//	m_dbK = (1.0 / 255.0) * pow(dbK, (double) nIter); 
//
//	outfile << m_dbK << endl;
//	outfile.close();
//}

// 
// Update the penalty function parameter
// 
void TGIEPGraph::UpdateK(double dbK, int nIter)
{
	//ofstream outfile;
	//if (nIter == 0)
	//	outfile.open("f1.txt",ios::out);
	//else
	//	outfile.open("f1.txt",ios::app);
	
	// Compute the edge strength histogram
	if (m_nEdgeSiteTotal < 0)
		ComputeEdgeStrengthHist();

	double dbRatio = nIter * SCHEDULE_STEP;// * m_nEdgeSiteTotal;  //QIN_MAJOR_DIFF
	if (dbRatio > 1)
		m_dbK = m_dbK * dbK;
	else
	{
		int i;
		for (i = 0; i < 256; i ++)
		{
			if (m_pEdgeHist[i] > dbRatio)
				break;
			dbRatio -= m_pEdgeHist[i];
		}
		if ( i == 256)
			m_dbK = m_dbEdgeStrengthMax;
		else
			m_dbK = (i + 0.5) * m_dbEdgeStrengthMax / 256;
	}
	//outfile << m_dbK << endl;
	//outfile.close();
}

//
// generate a struct for the specified boundary
//
TBoundary* TGIEPGraph::CreateNewBoundary(int nLabel1, int nLabel2, TChainCode* pChain, BOOL bOnGrid) { 
	if (m_enStatisticsMode == TRegionAdjacencyGraph::GAUSSIAN)
	{
		return new TGaussianStrengthPenaltyBoundary(this, nLabel1, nLabel2, pChain, bOnGrid); 
	}
//	else if (m_enStatisticsMode == TRegionAdjacencyGraph::POLARIMETRIC_WISHART)
//	{
//		return new TPolarimetricStrengthPenaltyBoundary(this, nLabel1, nLabel2, pChain, bOnGrid);
//	}
	else
	{
		throw std::runtime_error("Unknown statistics mode!@!@!@!@! You should feel bad for causing this.");
	}
}

// 
// Compute the maximum of edge strength
// 
double TGIEPGraph::ComputeMaxEdgeStrength()
{
	double dbVal, dbMax = 0;
	if (m_enConstructMode == PIXEL)
	{
		int nWidth = m_ppSrcImg[0]->Width();
		int nHeight = m_ppSrcImg[0]->Height();
		int r, c, i;
		for (int row = nHeight - 1; row >= 0; row --)
			for (int col = nWidth - 1; col >= 0; col --)
		{
			if (m_pMask && !(*m_pMask)(row, col))
				continue;

			for (i = 0; i < 4; i ++)  // only those neighbors before current wrt the scan order, so that an edge cannot be counted twice
			{
				r = row + neighbor8[i].Row;
				c = col + neighbor8[i].Col;
				if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
					continue;
				if (m_pMask && !(*m_pMask)(r, c))
					continue;
				dbVal = 0;
				for (int k = 0; k < m_nFeatureCnt; k++)
					dbVal += ((*m_ppSrcImg[k])(row, col) - (*m_ppSrcImg[k])(r, c))*((*m_ppSrcImg[k])(row, col) - (*m_ppSrcImg[k])(r, c));
				dbVal = sqrt(dbVal);
				if (dbVal > dbMax)
					dbMax = dbVal;
			}
		}
	}
	else
	{
		TChainCode *pChain;
		POINTex Pt;
		int nDir;
		TBoundary* pCurBoundary = (TBoundary*) GetEdges();
		while (pCurBoundary)
		{
			pChain = pCurBoundary->m_pChain;
			while (pChain)
			{
				pChain->GetFirstSite(Pt, nDir);
				do {
					dbVal = (*m_pGradient)(Pt.Row, Pt.Col);
					if (dbVal > dbMax)
						dbMax = dbVal;
				} while (pChain->GetNextSite(Pt, nDir));
				pChain = pChain->m_pNext;
			}
			pCurBoundary = (TBoundary*) pCurBoundary->m_pNext;
		}
	}

	return dbMax;
}

// 
// Compute the edge strength histogram
// 
void TGIEPGraph::ComputeEdgeStrengthHist(){
	m_nEdgeSiteTotal = 0;
	m_dbEdgeStrengthMax = ComputeMaxEdgeStrength();
	double dbHistBin = m_dbEdgeStrengthMax / 256;
	memset(m_pEdgeHist, 0, 256 * sizeof(double));
	int index;
	if (m_enConstructMode == PIXEL)
	{
		int nWidth = m_ppSrcImg[0]->Width();
		int nHeight = m_ppSrcImg[0]->Height();
		int r, c, i;
		double dbVal;
		for (int row = nHeight - 1; row >= 0; row --)
			for (int col = nWidth - 1; col >= 0; col --)
		{
			if (m_pMask && !(*m_pMask)(row, col))
				continue;

			for (i = 0; i < 4; i ++)  // only those neighbors before current wrt the scan order, so that an edge cannot be counted twice
			{
				r = row + neighbor8[i].Row;
				c = col + neighbor8[i].Col;
				if (r < 0 || r >= nHeight || c < 0 || c >= nWidth)
					continue;
				if (m_pMask && !(*m_pMask)(r, c))
					continue;
				dbVal = 0;
				for (int k = 0; k < m_nFeatureCnt; k++)
					dbVal += ((*m_ppSrcImg[k])(row, col) - (*m_ppSrcImg[k])(r, c))*((*m_ppSrcImg[k])(row, col) - (*m_ppSrcImg[k])(r, c));
				index = (int) (sqrt(dbVal) / dbHistBin);
				if (index >= 256)
					index = 255;
				m_pEdgeHist[index] ++;
				m_nEdgeSiteTotal ++;
			}
		}
	}
	else
	{
		TChainCode *pChain;
		POINTex Pt;
		int nDir;
		TBoundary* pCurBoundary = (TBoundary*) GetEdges();
		while (pCurBoundary)
		{
			pChain = pCurBoundary->m_pChain;
			while (pChain)
			{
				pChain->GetFirstSite(Pt, nDir);
				do {
					index = (int) ((*m_pGradient)(Pt.Row, Pt.Col) / dbHistBin);
					if (index >= 256)
						index = 255;
					m_pEdgeHist[index] ++;
				} while (pChain->GetNextSite(Pt, nDir));
				pChain = pChain->m_pNext;
			}

			m_nEdgeSiteTotal += pCurBoundary->m_nLength;
			pCurBoundary = (TBoundary*) pCurBoundary->m_pNext;
		}
	}
	
	for (index = 0; index < 256; index ++)
	{
		if ( m_pEdgeHist[index] )
			m_pEdgeHist[index] /= ((double) m_nEdgeSiteTotal);
	}
}

// 
// Update the penalty for binary relation of regions
//
void TGIEPGraph::UpdatePenalty()
{
	TStrengthPenaltyBoundary* pCurBoundary = (TStrengthPenaltyBoundary*) GetEdges();
	while (pCurBoundary)
	{
	    pCurBoundary->ComputePenaltyEnergy();
		pCurBoundary = (TStrengthPenaltyBoundary*) pCurBoundary->m_pNext;
	}
}

void TGIEPGraph::SetClassLabels(std::vector<int> vClassLabels) {
    if (vClassLabels.size() != this->GetVertexCnt()) {
        std::cout << "wrong number of inputs!!" << std::endl;
        return; // OH NO
    }

    int i=0;
    for(TRegion* pR = (TRegion*)this->GetVertices(); pR; pR = (TRegion*)pR->m_pNext){
        pR->m_nClassLabel = vClassLabels[i++];
    }
}