//
// Created by max on 2021-12-23.
//

#include "PyRAG.h"


PyRAG::PyRAG(py::array_t<float> srcImage, py::array_t<float> Gradient, py::array_t<int> Watershed, py::array_t<int> Landmask, int nClassCnt, int nBoundaryMode)
{

    if (srcImage.ndim() == 3)
        m_nFeatureCnt = srcImage.shape(2);
    else
        m_nFeatureCnt = 1;

    m_nClassCnt = nClassCnt;
    m_nBoundaryMode = nBoundaryMode;

    m_ppSrcImg = PyArray2GrayImageArray(srcImage);
    m_pGradient = PyArray2GrayImage(Gradient);

    TGrayImage<int>* pMap = PyArray2GrayImage(Watershed);
    m_pWatershedMap = new TGrayImage<int>(pMap);

    TGrayImage<int>* pLandmask = PyArray2GrayImage(Landmask);
    m_pMask = new TMonoImage(pLandmask);
    delete pLandmask;

    if (!m_nBoundaryMode){
        // construct the region adjacency graph
        m_pGraph = new TRegionMerger(0,
                                     m_pWatershedMap,
                                     m_pMask,
                                     m_ppSrcImg,
                                     m_nFeatureCnt,
                                     TRegionAdjacencyGraph::GAUSSIAN,
                                     m_pGradient,
                                     TRegionAdjacencyGraph::FULL);
    } else {
        // construct the region adjacency graph
        // PreprocessWatershed();
        m_pGraph = new TRegionMerger(0,
                                     m_pWatershedMap,
                                     m_pMask,
                                     m_ppSrcImg,
                                     m_nFeatureCnt,
                                     TRegionAdjacencyGraph::GAUSSIAN,
                                     m_pGradient,
                                     TRegionAdjacencyGraph::NO_INTERSECTIONS);
    }


    m_pGraph->ConstructGraph();

    /// create the boundary label map
    m_pBoundaryLabelMap = new TGrayImage<int>(pMap->Width(), pMap->Height());
    m_pBoundaryLabelMap->SetValue(-2);
    int edgeId=0;
    TBoundary* pCurBoundary = (TBoundary*) m_pGraph->GetEdges();
    while (pCurBoundary)
    {
        TChainCode *pChain = pCurBoundary->m_pChain;
        POINTex Pt;
        int nDir;

        while (pChain)
        {
            pChain->GetFirstSite(Pt, nDir);
            do {
                (*m_pBoundaryLabelMap)(Pt.Row, Pt.Col) = edgeId;

            } while (pChain->GetNextSite(Pt, nDir));

            pChain = pChain->m_pNext;
        }
        edgeId++;
        pCurBoundary = (TBoundary*) pCurBoundary->m_pNext;
    }
    delete pMap;

    // create the initial consecutive region label mapping
    for(int i=0; i< m_pGraph->GetVertexCnt(); i++){
        m_vConsecutiveMapping.push_back(i);
    }

    // Instantiate the GMM and allocate space for its parameters
    gmm = new GMM(m_nClassCnt, m_nFeatureCnt);

    // leave all the IRGS stuff uninitialized for now. It can be instantiated later using the
    // SetupIRGS() method.
    m_pClassifier=nullptr;
    m_bIRGSsetup=false;
    m_nIRGSIterElapsed=0;
}

PyRAG::~PyRAG(){
    delete m_pWatershedMap;
    delete m_pBoundaryLabelMap;
    delete m_pGraph;
    for (int i = 0; i < m_nFeatureCnt; i++)
        delete m_ppSrcImg[i];
    delete[] m_ppSrcImg;
    delete m_pGradient;
    delete gmm;
    delete m_pClassifier;
}

py::array_t<int> PyRAG::GetOriginalWatershed() const{
    py::array_t<int> aWatershed = GrayImage2PyArray(m_pWatershedMap);
    return aWatershed;
}

// Get the graph region map. The graph region map is initialized to the watershed map, but changes as soon as a merging
// iteration occurs. In this version the map is relabeled to make the labels consecutive (no gaps), and the consecutive
// label mapping is returned as well
std::tuple<py::array_t<int>,py::array_t<int>> PyRAG::GetWatershed() const {
    TGrayImage<int>* pMap = new TGrayImage<int>(m_pGraph->GetMap());
    TEstimateUtil pEstUtil;
    int nRegions;
    int* pMapping = pEstUtil.GenerateConsecutiveMapping(nRegions, pMap, m_pMask);
    pEstUtil.RelabelMap(pMap, m_pMask);
    py::array_t<int> aWatershed = GrayImage2PyArray(pMap);

    py::array_t<int> LabelMapping(nRegions);

    auto ra = LabelMapping.mutable_unchecked<1>();
    for(int i=0; i<nRegions;i++){
        ra(i) = pMapping[i];
    }

    delete pMap;
    delete[]pMapping;
    return std::make_tuple(aWatershed, LabelMapping);
}

py::array_t<int> PyRAG::GetBoundaryMap() {
    py::array_t<int> aBoundaryMap = GrayImage2PyArray(m_pBoundaryLabelMap);
    return aBoundaryMap;
}

std::tuple<py::array_t<int>, py::array_t<int>, py::array_t<int>, py::array_t<int>> PyRAG::GetChainCodes() const
{
    std::vector<int> startpos, endpos, codes, edgeidx;
    int edgeId=0;

    TBoundary* pCurBoundary = (TBoundary*) m_pGraph->GetEdges();
    while (pCurBoundary)
    {
        TChainCode *pChain = pCurBoundary->m_pChain;
        POINTex Pt;
        int nDir;

        startpos.push_back(pChain->m_StartPt.Col);
        startpos.push_back(pChain->m_StartPt.Row);
        endpos.push_back(pChain->m_EndPt.Col);
        endpos.push_back(pChain->m_EndPt.Row);

        pChain->GetFirstSite(Pt, nDir);
        while (pChain)
        {
            pChain->GetFirstSite(Pt, nDir);
            while (pChain->GetNextSite(Pt, nDir))
            {
                edgeidx.push_back(edgeId);
                codes.push_back(nDir);
            }

            pChain = pChain->m_pNext;
        }
        edgeId++;
        pCurBoundary = (TBoundary*) pCurBoundary->m_pNext;
    }

    // this incurs a copy but it is worth it for the convenience
    py::array_t<int> astartposx = py::array_t<int>(startpos.size(), startpos.data());
    py::array_t<int> astartposy = py::array_t<int>(endpos.size(), endpos.data());
    py::array_t<int> acodes     = py::array_t<int>(codes.size(),     codes.data());
    py::array_t<int> aedgeidx   = py::array_t<int>(edgeidx.size(),   edgeidx.data());

    return std::make_tuple(astartposx, astartposy, aedgeidx, acodes);
}

py::array_t<int> PyRAG::GetBoundaryLengths(){

    int n_edge = m_pGraph->GetEdgeCnt();
    py::array_t<int, py::array::c_style> pyArray({2*n_edge});

    auto ra = pyArray.mutable_unchecked<1>();
    int i=0;
    for(TBoundary* pB = (TBoundary*) m_pGraph->GetEdges(); pB; pB=(TBoundary*)pB->m_pNext){
        int length = pB->m_nLength;
        ra(i++) = length;
        ra(i++) = length; // double count because edge length is same for both directions
    }

    return pyArray;
}

py::array_t<int> PyRAG::GetRegionClassLabels() const{
    int n_node = m_pGraph->GetVertexCnt();
    py::array_t<float, py::array::c_style> pyArray({n_node});

    auto ra = pyArray.mutable_unchecked<1>();
    for(TRegion* pR = (TRegion*) m_pGraph->GetVertices(); pR; pR=(TRegion*)pR->m_pNext){
        ra(m_vConsecutiveMapping[pR->m_nLabel]) = pR->m_nClassLabel;

    }
    return pyArray;
}

py::array_t<float> PyRAG::GetRegionCentroids() const
{
    int n_node = m_pGraph->GetVertexCnt();
    py::array_t<float, py::array::c_style> pyArray({2, n_node});

    auto ra = pyArray.mutable_unchecked<2>();
    for(TRegion* pR = (TRegion*) m_pGraph->GetVertices(); pR; pR=(TRegion*)pR->m_pNext){
        ra(0, m_vConsecutiveMapping[pR->m_nLabel]) = (float) pR->m_dbXMean;
        ra(1, m_vConsecutiveMapping[pR->m_nLabel]) = (float) pR->m_dbYMean;
    }
    return pyArray;
}

py::array_t<float> PyRAG::GetRegionMeans() const
{
    int n_node = m_pGraph->GetVertexCnt();
    py::array_t<float, py::array::c_style> pyArray({m_nFeatureCnt, n_node});

    auto ra = pyArray.mutable_unchecked<2>();
    for(TGaussianRegion* pR = (TGaussianRegion*) m_pGraph->GetVertices(); pR; pR=(TGaussianRegion*)pR->m_pNext){
        for(int i=0;i<m_nFeatureCnt;i++){
            ra(i, m_vConsecutiveMapping[pR->m_nLabel]) = (float) (*pR->m_pMean)(i);
        }
    }
    return pyArray;
}

py::array_t<float> PyRAG::GetRegionCoMeans() const
{
    int n_node = m_pGraph->GetVertexCnt();
    py::array_t<float, py::array::c_style> pyArray({m_nFeatureCnt, m_nFeatureCnt, n_node});

    auto ra = pyArray.mutable_unchecked<3>();
    for(TGaussianRegion* pR = (TGaussianRegion*) m_pGraph->GetVertices(); pR; pR=(TGaussianRegion*)pR->m_pNext){
        for(int i=0;i<m_nFeatureCnt;i++){
            ra(i, i, m_vConsecutiveMapping[pR->m_nLabel]) = (float) (*pR->m_pCoMean)(i,i);
            for(int j=i+1; j<m_nFeatureCnt;j++){
                ra(i, j, m_vConsecutiveMapping[pR->m_nLabel]) = (float) (*pR->m_pCoMean)(i,j);
                ra(j, i, m_vConsecutiveMapping[pR->m_nLabel]) = (float) (*pR->m_pCoMean)(i,j);
            }
        }
    }
    return pyArray;
}

py::array_t<int> PyRAG::GetRegionPixelCounts() const
{
    int n_node = m_pGraph->GetVertexCnt();
    py::array_t<int, py::array::c_style> pyArray({n_node});

    auto ra = pyArray.mutable_unchecked<1>();
    for(TGaussianRegion* pR = (TGaussianRegion*) m_pGraph->GetVertices(); pR; pR=(TGaussianRegion*)pR->m_pNext){
        ra(m_vConsecutiveMapping[pR->m_nLabel]) = pR->m_nPixelCnt;
    }
    return pyArray;
}

py::array_t<float> PyRAG::GetGMMUnary() const
{
    int n_node = m_pGraph->GetVertexCnt();
    int n_class = m_nClassCnt;
    py::array_t<float, py::array::c_style> pyArray({n_class, n_node});

    auto ra = pyArray.mutable_unchecked<2>();

    for(TRegion* pR = (TRegion*) m_pGraph->GetVertices(); pR; pR=(TRegion*)pR->m_pNext){
        int l = m_vConsecutiveMapping[pR->m_nLabel];
        for(int i=0; i<n_class; i++){
            ra(i,l) = (float) GetUnaryEnergy(pR, i);
        }
    }

    return pyArray;
}

py::array_t<float> PyRAG::GetEdgePenalty()const
{
    int n_edge = m_pGraph->GetEdgeCnt();
    py::array_t<float, py::array::c_style> pyArray({2*n_edge});

    auto ra = pyArray.mutable_unchecked<1>();
    int i=0;
    for(TBoundary* pB = (TBoundary*) m_pGraph->GetEdges(); pB; pB=(TBoundary*)pB->m_pNext){
        double strength = ((TStrengthPenaltyBoundary*) pB)->m_dbPenaltyEnergy;
        ra(i++) = strength;
        ra(i++) = strength; // double count because edge strength is same for both directions
    }

    return pyArray;
}


bool PyRAG::UpdateGMMParam()
{
    bool bSucFlag = true;

    int *nSampleCnt = new int[m_nClassCnt];
    int i, j, k;
    for ( k = 0; k < m_nClassCnt; k++ )
    {
        gmm->m_ppMu[k]->SetZero();
        gmm->m_ppCovInv[k]->SetZero();
        nSampleCnt[k] = 0;
    }

    // Mean and Variance
    int nFeatureCnt = gmm->m_ppMu[0]->GetElementCnt();
    double *dbMu, **dbCoMu;
    TGaussianRegion* pCurrent;
    TVertex* pVertex = m_pGraph->GetVertices();
    while (pVertex)
    {
        pCurrent = (TGaussianRegion*) pVertex;
        k = pCurrent->m_nClassLabel;

        dbMu = gmm->m_ppMu[k]->m_ppData[0];
        dbCoMu = gmm->m_ppCovInv[k]->m_ppData;


        for ( i = 0; i < nFeatureCnt; i++ )
        {
            dbMu[i] += (*pCurrent->m_pMean)(i) * pCurrent->m_nPixelCnt;
            for ( j = i; j < nFeatureCnt; j++ )
                dbCoMu[i][j] += (*pCurrent->m_pCoMean)(i, j) * pCurrent->m_nPixelCnt;
        }
        nSampleCnt[k] += pCurrent->m_nPixelCnt;

        pVertex = pVertex->m_pNext;
    }

    int nTotalCnt = 0;
    for ( k = 0; k < m_nClassCnt; k++ )
    {

        if (nSampleCnt[k])
        {
            dbMu = gmm->m_ppMu[k]->m_ppData[0];
            dbCoMu = gmm->m_ppCovInv[k]->m_ppData;
            for ( i = 0; i < nFeatureCnt; i++ )
            {
                dbMu[i] /= nSampleCnt[k];
                for ( j = i; j < nFeatureCnt; j++ )
                    dbCoMu[i][j] /= nSampleCnt[k];
            }
            for ( i = 0; i < nFeatureCnt; i++ )
            {
                dbCoMu[i][i] -= dbMu[i] * dbMu[i];
                for ( j = i + 1; j < nFeatureCnt; j++ )
                {
                    dbCoMu[i][j] -= dbMu[i] * dbMu[j];
                    dbCoMu[j][i] = dbCoMu[i][j];
                }
            }
        }
        else
        {
            bSucFlag = false;
            nSampleCnt[k]++;
        }

        if ( gmm->m_ppCovInv[k]->Eigen() <= 1e-12 )
            for ( i = 0; i < nFeatureCnt; i++ )
                (*gmm->m_ppCovInv[k])(i, i) += 0.001;

        gmm->m_dbCovDet[k] = sqrt(gmm->m_ppCovInv[k]->Determinant());
        gmm->m_dbLnSigma[k] = log(gmm->m_dbCovDet[k]);
        gmm->m_ppCovInv[k]->Inv(gmm->m_ppCovInv[k]);

        nTotalCnt += nSampleCnt[k];
    }

    // Priors
    for ( k = 0; k < m_nClassCnt; k++ )
    {
        gmm->m_dbPrior[k] = ((double) nSampleCnt[k]) / ((double) nTotalCnt);
        gmm->m_dbLnPrior[k] = log(gmm->m_dbPrior[k]); //QIN_ADD
    }

    delete[] nSampleCnt; //QIN_MOD

    return bSucFlag;
}


bool PyRAG::UpdateBoundaryParamFisher(double dbK, double dbBeta1, double dbBeta2, int it) {

    double dbFisher = ComputeMinFisher(gmm);
    m_pGraph->UpdateBetaFisher(dbK, dbBeta1, dbBeta2, dbFisher, it);
    m_pGraph->UpdatePenalty();
    return true;
}

bool PyRAG::UpdateBoundaryParam(double dbK, int it) {

    m_pGraph->UpdateK(dbK, it);
    m_pGraph->UpdateBeta(m_pGraph->GetBoundaryPopulationRatio(), 1);
    m_pGraph->UpdatePenalty();
    return true;
}

PyGMM PyRAG::GetGMMParams() const{
    PyGMM* pygmm = new PyGMM(gmm);
    return *pygmm;
}

// gets the TRegionMerger object ready to do IRGS steps by instatiating its TRegionClassifier object
void PyRAG::SetupIRGS() {

    // default parameters
    double dbInitRatio = 0.000001;
    int nInit1 = 25;
    int nInit2 = 20;
    int nInitIterFinal = 500;
    int	nSAMode = TMarkovNetClassifier::SINGLE_STEP | TMarkovNetClassifier::NO_COOLING; // equal to 3
    double	dbSAInitT = 0.01;
    int	nSAIter = 1;

    if (m_pClassifier) delete m_pClassifier;
    m_pClassifier = new TGaussianRegionClustering(TRegionClustering::GIEP,
                                                  TRegionClustering::KM_REGION,
                                                  TRegionClustering::RAND2,
                                                  dbInitRatio, nInit1, nInit2,
                                                  nInitIterFinal,
                                                  nSAMode,
                                                  dbSAInitT,
                                                  nSAIter,
                                                  m_nClassCnt,
                                                  false);

    m_pClassifier->SetGraph(m_pGraph);
    m_pGraph->SetClassifier(m_pClassifier);

    m_bIRGSsetup = true;
    m_nIRGSIterElapsed=0;
}

void PyRAG::IRGSStep(double dbK, double dbBeta1, double dbBeta2, bool bVerbose)
{
    if (!m_bIRGSsetup) SetupIRGS();

//    double	dbK = 1.1;
//    double	dbBeta1 = 3;
//    double	dbBeta2 = 0.4;

    if (bVerbose){
        py::print("------------------------------------------------------");
        py::print("Iteration: "+ std::to_string(m_nIRGSIterElapsed) );
        py::print("Num. regions before merging: " + std::to_string(m_pGraph->GetVertexCnt()));
    }

    // do the merging step
    std::tuple<double, double, bool> result = m_pGraph->MergingStep(dbK, dbBeta1, dbBeta2, m_nIRGSIterElapsed++);

    if (bVerbose){
        py::print("Adjusted Beta: " + std::to_string(get<0>(result)));
        py::print("Fisher: " + std::to_string(get<1>(result)));
        py::print("Num. regions after merging: " + std::to_string(m_pGraph->GetVertexCnt()));
    }

    // copy the GMM parameters from pClassifier to the gmm object in PyRAG
    TGaussianRegionClustering* pClassifier = (TGaussianRegionClustering*) m_pGraph->GetClassifier();
    gmm->CopyFromGaussClassifier(pClassifier);

    bool changed = get<2>(result);

    if(changed){
        // update the consecutive mapping
        TEstimateUtil pEstUtil;
        int nRegions;
        int* pMapping = pEstUtil.GenerateConsecutiveMapping(nRegions, m_pGraph->GetMap(), m_pMask);
        for(int i=0; i<nRegions; i++){
            m_vConsecutiveMapping[i] = pMapping[i];
        }

        // update the boundary map
        m_pBoundaryLabelMap->SetValue(-2);
        int edgeId=0;
        TStrengthPenaltyBoundary* pCurBoundary = (TStrengthPenaltyBoundary*) m_pGraph->GetEdges();
        while (pCurBoundary)
        {

            TChainCode *pChain = pCurBoundary->m_pChain;
            POINTex Pt;
            int nDir;

            while (pChain)
            {
                pChain->GetFirstSite(Pt, nDir);
                do {
                    (*m_pBoundaryLabelMap)(Pt.Row, Pt.Col) = edgeId;

                } while (pChain->GetNextSite(Pt, nDir));

                pChain = pChain->m_pNext;
            }
            edgeId++;
            pCurBoundary = (TStrengthPenaltyBoundary*) pCurBoundary->m_pNext;
        }
    }
}

// returns an array of size (2, 2*num_edges) which defines the graph connectivity by giving the indices of the source
// and destination node for each edge
py::array_t<int> PyRAG::GetEdges() const
{
    int n_edge = 2*m_pGraph->GetEdgeCnt();
    py::array_t<int, py::array::c_style> pyArray({2, n_edge});

    auto ra = pyArray.mutable_unchecked<2>();
    int i=0;
    for(TEdge* pEdge = m_pGraph->GetEdges(); pEdge; pEdge=pEdge->m_pNext){
        ra(0,i)   = pEdge->m_nLabel1;
        ra(1,i++) = pEdge->m_nLabel2;
        ra(0,i)   = pEdge->m_nLabel2;
        ra(1,i++) = pEdge->m_nLabel1;
    }

    return pyArray;
}

// Generates a LUT for the mapping of old regions to new regions for all region merging steps so far
py::array_t<int> PyRAG::GetMergingLUT()
{
    TEstimateUtil* pEstUtil;
    int N;
    int* MergeLUT = pEstUtil->GenerateMergingLUT(N, m_pWatershedMap, m_pGraph->GetMap(),m_pMask);
    py::array_t<int> aMerged ({N});
    auto ra = aMerged.mutable_unchecked<1>();
    for(int i=0; i<N; i++){
        ra(i) = MergeLUT[i];
    }
    delete[] MergeLUT;

    return aMerged;
}

// Generates a LUT for the mapping of old regions to new regions during a region merging step
// aMapOld is the region map before merging, aMapNew is the region map after merging.
py::array_t<int> PyRAG::GetMergingLUTFromMaps(py::array_t<int> aMapOld, py::array_t<int> aMapNew)
{
    TGrayImage<int>* pMapOld = PyArray2GrayImage(aMapOld);
    TGrayImage<int>* pMapNew = PyArray2GrayImage(aMapNew);

    TEstimateUtil* pEstUtil;
    int N;
    int* MergeLUT = pEstUtil->GenerateMergingLUT(N, pMapOld, pMapNew, m_pMask);
    py::array_t<int> aMerged ({N});
    auto ra = aMerged.mutable_unchecked<1>();
    for(int i=0; i<N; i++){
        ra(i) = MergeLUT[i];
    }
    delete[] MergeLUT;
    delete pMapOld;
    delete pMapNew;

    return aMerged;
}

void PyRAG::PreprocessWatershed() const
{
    int nWidth = m_pWatershedMap->Width();
    int nHeight = m_pWatershedMap->Height();
    int rr, cc;

    for (int r = 1; r < nHeight-1; r ++)
        for (int c = 1; c < nWidth-1; c ++)
        {
            if (m_pMask && !(*m_pMask)(r, c)) // ignore masked pixels
                continue;

            if ((*m_pWatershedMap)(r,c) > 0) // ignore non-boundary pixels
                continue;

            int nAdjacentLabels = 0;
            int nAdjacentBoundaries = 0;
            int nLabel=-2, val;

            for(int i=0; i<8; i++){

                rr = r + neighbor8[i].Row;
                cc = c + neighbor8[i].Col;

                if ((rr < 0) || (cc < 0) || (rr >= nHeight) || (cc >= nWidth)) continue;

                val = (*m_pWatershedMap)(rr,cc);

                if ( i%2==0 && val < 0 )
                    nAdjacentBoundaries++; // count how many 4-connected boundary pixels there are

                if ((val >= 0) && (val != nLabel))
                {
                    nAdjacentLabels++;
                    nLabel = val;
                }
            }

            if ((nAdjacentBoundaries == 2) && (nAdjacentLabels == 1))
                (*m_pWatershedMap)(r,c) = nLabel;
        }
}

void PyRAG::SetClassLabelsByIndex(py::array_t<int> ClassLabels, py::array_t<int> NodeIndices) {

    int L = ClassLabels.shape(0);
    if (NodeIndices.shape(0) != L) throw std::invalid_argument("SetClassLabelsByIndex: incompatible shapes");

    auto ra = ClassLabels.unchecked<1>();
    auto rb = NodeIndices.unchecked<1>();

    for(int i=0;i<L;i++){
        if (!m_pGraph->SetClassLabelByIndex(ra(i), rb(i))){
            throw std::invalid_argument("SetClassLabelsByIndex: Tried to access a region that does not exist!!");
        }

    }

    // Break the links between regions belonging different classes, so
    // that they cannot be merged afterwards
    TEdge* pCur = m_pGraph->GetEdges();
    while (pCur)
    {
        if (((TRegion*) pCur->m_pNeighbor1)->m_nClassLabel !=
            ((TRegion*) pCur->m_pNeighbor2)->m_nClassLabel)
            ((TBoundary*) pCur)->m_bActive = false;  // not active
        else
            ((TBoundary*) pCur)->m_bActive = true;  // active
        pCur = pCur->m_pNext;
    }
}

py::array_t<int> PyRAG::GetConnectedClassMapBoundary()
{
    TGrayImage<int>* tpMap = m_pGraph->GetMap();
    TMonoImage* pMask = m_pGraph->GetMask();

    // Get maximum number of remaining segments
    TEstimateUtil* pEstUtil = new TEstimateUtil();
    int nMaxCnt = pEstUtil->GetMaxLabel(tpMap, pMask) + 1;
    delete pEstUtil;

    // Generate the label mapping
    int* nLabels = new int[nMaxCnt];
    memset(nLabels, 0, nMaxCnt * sizeof(int));
    TRegion* pVertex = (TRegion*) m_pGraph->GetVertices();
    while (pVertex)
    {
        nLabels[pVertex->m_nLabel] = pVertex->m_nClassLabel;
        pVertex = (TRegion*) pVertex->m_pNext;
    }

    // Obtain the classification map
    int nWidth = tpMap->Width();
    int nHeight = tpMap->Height();
    TGrayImage<int>* pMap = new TGrayImage<int>(tpMap);
    int l, r, c;
    for ( r = 0; r < nHeight; r++ )
        for ( c = 0; c < nWidth; c++ )
        {
            if ( pMask && !(*pMask)(r, c) )
                continue;
            if ( ( l = (*tpMap)(r, c) ) >= 0 )  // class label
                (*pMap)(r, c) = nLabels[l];
        }

    delete[] nLabels;

    int nRow, nCol, l2, i;
    for ( r = 0; r < nHeight; r++ )
        for ( c = 0; c < nWidth; c++ )
        {
            if ( pMask && !(*pMask)(r, c) )
                continue;
            if ( (*pMap)(r, c) < 0 )
            {
                // Checks pixels around current boundary pixel.
                // If all non-boundary pixels that touch the current one
                // belong to the same class, then the current one also
                // belongs to that class and is no longer a boundary point.
                l2 = -1;
                for ( i = 0; i < 8; i++ )
                {
                    nRow = r + neighbor8[i].Row;
                    nCol = c + neighbor8[i].Col;

                    if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
                        continue;
                    if ( pMask && !(*pMask)(nRow, nCol) )  // not valid
                        continue;
                    if ( ( l = (*pMap)(nRow, nCol)) < 0 )  // also boundary
                        continue;
                    if ( l != l2 && l2 >= 0 )
                        break;

                    l2 = l;
                }

                if ( i >= 8 && l2 >= 0 )  // inside the same class
                    (*pMap)(r, c) = l2;
            }
        }

    TGrayImage<FLOAT>** ppFeatures = m_pGraph->GetImageSource();
    int nFeatureCnt = m_pGraph->GetFeatureCnt();

    TVector* pVecFeature = new TVector(nFeatureCnt, false);  // feature vector
    double* dbVecProb = new double[m_nClassCnt];
    TVector* pVecTemp1 = new TVector(nFeatureCnt, false);  // a temperoray vector
    TVector* pVecTemp2 = new TVector(nFeatureCnt, false);  // a temperoray vector

    double dbMin;
    int nRefineCnt = 1, j, k, nMinIdx;

    int* nClassLabel = new int[m_nClassCnt]; // indicates existence of class in image.
    int* nPixelClassLabel = new int[m_nClassCnt]; // indicates existence of class among neighbours of pixel.
    int* nClassLabelUsed; // pointer to which of the above will be used when classifying boundary pixels.
    bool bIndepBoundary; // whether a pixel is an independent boundary pixel with no other pixels of any class beside it.
    // In these cases, the pixel must be classified on its own from amongst all the classes in image.

    TVector** ppMu = new TVector*[m_nClassCnt];
    TMatrix** ppCovInv = new TMatrix*[m_nClassCnt];
    for ( i = 0; i < m_nClassCnt; i++ )
    {
        ppMu[i] = new TVector(nFeatureCnt, false);
        ppCovInv[i] = new TMatrix(nFeatureCnt, nFeatureCnt);
    }

    double* dbPrior = new double[m_nClassCnt];
    double* dbCovDet = new double[m_nClassCnt];

    for ( int m = 0; m < nRefineCnt; m++ )  // several iterations to refine
    {
        for ( i = 0; i < m_nClassCnt; i++ )
        {
            ppMu[i]->SetZero();
            ppCovInv[i]->SetZero();
        }
        memset(dbPrior, 0, sizeof(double) * m_nClassCnt);

        for ( r = 0; r < nHeight; r++ )
            for ( c = 0; c < nWidth; c++ )
            {
                if ( pMask && !(*pMask)(r, c) )
                    continue;
                if ( ( k = (*pMap)(r, c) ) < 0 )
                    continue;

                for ( i = 0; i < nFeatureCnt; i++ )
                {
                    (*ppMu[k])(i) += (*ppFeatures[i])(r, c);
                    (*ppCovInv[k])(i, i) += (*ppFeatures[i])(r, c) * (*ppFeatures[i])(r, c);
                    for ( j = i + 1; j < nFeatureCnt; j++ )
                        (*ppCovInv[k])(i, j) += (*ppFeatures[i])(r, c) * (*ppFeatures[j])(r, c);
                }

                dbPrior[k]++;
            }

        double dbPriTotal = 0.0;
        for ( k = 0; k < m_nClassCnt; k++ )
        {
            if ( dbPrior[k] > 0 )
            {
                for ( i = 0; i < nFeatureCnt; i++ )
                    (*ppMu[k])(i) /= dbPrior[k];  // mean

                for ( i = 0; i < nFeatureCnt; i++ )  // covariance
                {
                    (*ppCovInv[k])(i, i) = (*ppCovInv[k])(i, i) / dbPrior[k] - (*ppMu[k])(i) * (*ppMu[k])(i);
                    for ( j = i + 1; j < nFeatureCnt; j++ )
                    {
                        (*ppCovInv[k])(i, j) = (*ppCovInv[k])(i, j) / dbPrior[k] - (*ppMu[k])(i) * (*ppMu[k])(j);
                        (*ppCovInv[k])(j, i) = (*ppCovInv[k])(i, j);
                    }
                }
                nClassLabel[k] = 1;
            }
            else
            {
                nClassLabel[k] = 0;
                dbPrior[k]++;
            }

            if ( ppCovInv[k]->Eigen() <= 1e-9 )
                for ( i = 0; i < nFeatureCnt; i++ )
                    (*ppCovInv[k])(i, i) += 1e-6;
            dbCovDet[k] = sqrt(ppCovInv[k]->Determinant());  // determinant
            ppCovInv[k]->Inv(ppCovInv[k]);  // inverse of covariance

            dbPriTotal += dbPrior[k];
        }

        for ( k = 0; k < m_nClassCnt; k++ )
            dbPrior[k] = dbPrior[k] / dbPriTotal;  // prior

        for ( r = 0; r < nHeight; r++ )
            for ( c = 0; c < nWidth; c++ )
            {
                if ( pMask && !(*pMask)(r, c) )
                    continue;

                // Classify the boundary pixels using a pixel based
                // MRF model. First, check to see if it's a boundary pixel
                // and then if it is, determine what classes the neighbouring
                // pixels are. Choose from amongst these neighbouring classes
                // to classify the boundary pixel.

                memset(nPixelClassLabel, 0, sizeof(int) * m_nClassCnt); // reset the neighbouring class indicator list.

                if (IsBoundary(pMap, pMask, r, c, nPixelClassLabel, bIndepBoundary))
                {
                    if (bIndepBoundary) // pixel has no neighbours with valid class, so classify it using all classes as candidates.
                    {
                        nClassLabelUsed = nClassLabel;
                    }
                    else // pixel has neighbours with valid class, so classify it using only those neighbouring classes.
                    {
                        nClassLabelUsed = nPixelClassLabel;
                    }

                    // Compute the probabilities
                    for ( i = 0; i < nFeatureCnt; i++ )
                        (*pVecFeature)(i) = (*ppFeatures[i])(r, c);

                    for ( i = 0; i < m_nClassCnt; i++ )
                    {
                        if ( nClassLabelUsed [i] )
                        {
                            for ( j = 0; j < nFeatureCnt; j++ )
                                (*pVecTemp1)(j) = (*pVecFeature)(j) - (*ppMu[i])(j);
                            pVecTemp1->Mul(pVecTemp2, ppCovInv[i]);
                            dbVecProb[i] = ( ( pVecTemp2->VecMul(pVecTemp1) / 2 + log(dbCovDet[i]) ) / nFeatureCnt );  //QIN_nFeatureCnt
                        }
                        else
                            dbVecProb[i] = 999999999999;
                    }

                    // Account spatial context constraint
                    for ( j = 0; j < 8; j++ )
                    {
                        nRow = r + neighbor8[j].Row;
                        nCol = c + neighbor8[j].Col;

                        if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
                            continue;
                        if (pMask && !(*pMask)(nRow, nCol))  // not valid
                            continue;
                        if ( ( l = (*pMap)(nRow, nCol) ) < 0 )  // also boundary
                            continue;

                        dbVecProb[l] -= ((TGIEPGraph*) m_pGraph)->GetCurBeta(); // QIN_nFeatureCnt

                    }

                    // Determine the label
                    dbMin = 999999999999;  // a large enough number
                    for ( i = 0; i < m_nClassCnt; i++ )
                    {
                        if ( dbVecProb[i]< dbMin )
                        {
                            dbMin = dbVecProb[i];
                            nMinIdx = i;
                        }
                    }

                    (*pMap)(r, c) = nMinIdx;
                }
            }
    }

    for ( i = 0; i < m_nClassCnt; i++ )
    {
        delete ppMu[i];
        delete ppCovInv[i];
    }
    delete[] ppMu;
    delete[] ppCovInv;
    delete[] dbPrior;
    delete[] dbCovDet;

    delete[] nClassLabel;
    delete[] nPixelClassLabel;

    delete pVecFeature;
    delete[] dbVecProb;
    delete pVecTemp1;
    delete pVecTemp2;

    py::array_t<int> aMap = GrayImage2PyArray(pMap);
    delete pMap;
    return aMap;
}


bool PyRAG::IsBoundary(TGrayImage<int>* pMap, TMonoImage* pMask, int nRow, int nCol, int*& nClassLabel, bool& bIndepBoundary)
{
    int nWidth = pMap->Width();
    int nHeight = pMap->Height();
    int i, r, c;
    int l = (*pMap)(nRow, nCol);
    int l2;

    if (l < 0)
    {
        bIndepBoundary = true;

        for ( i = 0; i < 8; i++ )
        {
            r = nRow + neighbor8[i].Row;
            c = nCol + neighbor8[i].Col;

            if ( r < 0 || r >= nHeight || c < 0 || c >= nWidth )
                continue;
            if ( pMask && !(*pMask)(r, c) )
                continue;

            l2 = (*pMap)(r, c);

            if (l2 < 0)
                continue;

            nClassLabel[l2] = 1; // this class exists in the neighbourhood.
            bIndepBoundary = false; // at least one valid class beside this pixel.

        }
        return true;
    }

    bIndepBoundary = false;
    return false;
}

py::array_t<int> PyRAG::GetConnectedClassMapGreedy()
{

    TGrayImage<int>* tpMap = m_pGraph->GetMap();
    TMonoImage* pMask = m_pGraph->GetMask();

    // Get maximum number of remaining segments
    TEstimateUtil* pEstUtil = new TEstimateUtil();
    int nMaxCnt = pEstUtil->GetMaxLabel(tpMap, pMask) + 1;
    delete pEstUtil;

    // Generate the label mapping
    int* nLabels = new int[nMaxCnt];
    memset(nLabels, 0, nMaxCnt * sizeof(int));
    TRegion* pVertex = (TRegion*) m_pGraph->GetVertices();
    while (pVertex)
    {
        nLabels[pVertex->m_nLabel] = pVertex->m_nClassLabel;
        pVertex = (TRegion*) pVertex->m_pNext;
    }

    // Obtain the classification map
    int nWidth = tpMap->Width();
    int nHeight = tpMap->Height();
    TGrayImage<int>* pMap = new TGrayImage<int>(tpMap);
    int l, r, c;
    for ( r = 0; r < nHeight; r++ )
        for ( c = 0; c < nWidth; c++ )
        {
            if ( pMask && !(*pMask)(r, c) )
                continue;
            if ( ( l = (*tpMap)(r, c) ) >= 0 )  // class label
                (*pMap)(r, c) = nLabels[l];
        }

    delete[] nLabels;

    int nRow, nCol, l2, i;
    for ( r = 0; r < nHeight; r++ )
        for ( c = 0; c < nWidth; c++ )
        {
            if ( pMask && !(*pMask)(r, c) )
                continue;
            if ( (*pMap)(r, c) < 0 )
            {
                l2 = -1;
                for ( i = 0; i < 8; i++ )
                {
                    nRow = r + neighbor8[i].Row;
                    nCol = c + neighbor8[i].Col;

                    if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
                        continue;
                    if ( pMask && !(*pMask)(nRow, nCol) )  // not valid
                        continue;
                    if ( ( l = (*pMap)(nRow, nCol)) < 0 )  // also boundary
                        continue;
                    if ( l != l2 && l2 >= 0 )
                        break;

                    l2 = l;
                }

                if ( i >= 8 && l2 >= 0 )  // inside the same class
                    (*pMap)(r, c) = l2;
                else {
                    // greedily assign to the first valid class around the boundary
                    // note that this is a totally random and unprincipled way to do things
                    for ( i = 0; i < 8; i++ )
                    {
                        nRow = r + neighbor8[i].Row;
                        nCol = c + neighbor8[i].Col;

                        if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
                            continue;
                        if ( pMask && !(*pMask)(nRow, nCol) )  // not valid
                            continue;
                        if ( ( l = (*pMap)(nRow, nCol)) >= 0 )
                        {
                            (*pMap)(r, c) = l;
                            break;
                        }

                    }
                }
            }
        }
    py::array_t<int> aMap = GrayImage2PyArray(pMap);
    delete pMap;
    return aMap;
}


py::array_t<int> PyRAG::GetConnectedClassMap()
{

    TGrayImage<int>* tpMap = m_pGraph->GetMap();
    TMonoImage* pMask = m_pGraph->GetMask();

    // Get maximum number of remaining segments
    TEstimateUtil* pEstUtil = new TEstimateUtil();
    int nMaxCnt = pEstUtil->GetMaxLabel(tpMap, pMask) + 1;
    delete pEstUtil;

    // Generate the label mapping
    int* nLabels = new int[nMaxCnt];
    memset(nLabels, 0, nMaxCnt * sizeof(int));
    TRegion* pVertex = (TRegion*) m_pGraph->GetVertices();
    while (pVertex)
    {
        nLabels[pVertex->m_nLabel] = pVertex->m_nClassLabel;
        pVertex = (TRegion*) pVertex->m_pNext;
    }

    // Obtain the classification map
    int nWidth = tpMap->Width();
    int nHeight = tpMap->Height();
    TGrayImage<int>* pMap = new TGrayImage<int>(tpMap);
    int l, r, c;
    for ( r = 0; r < nHeight; r++ )
        for ( c = 0; c < nWidth; c++ )
        {
            if ( pMask && !(*pMask)(r, c) )
                continue;
            if ( ( l = (*tpMap)(r, c) ) >= 0 )  // class label
                (*pMap)(r, c) = nLabels[l];
        }

    delete[] nLabels;

    int nRow, nCol, l2, i;
    for ( r = 0; r < nHeight; r++ )
        for ( c = 0; c < nWidth; c++ )
        {
            if ( pMask && !(*pMask)(r, c) )
                continue;
            if ( (*pMap)(r, c) < 0 )
            {
                // Checks pixels around current boundary pixel.
                // If all non-boundary pixels that touch the current one
                // belong to the same class, then the current one also
                // belongs to that class and is no longer a boundary point.
                l2 = -1;
                for ( i = 0; i < 8; i++ )
                {
                    nRow = r + neighbor8[i].Row;
                    nCol = c + neighbor8[i].Col;

                    if (nRow < 0 || nRow >= nHeight || nCol < 0 || nCol >= nWidth)
                        continue;
                    if ( pMask && !(*pMask)(nRow, nCol) )  // not valid
                        continue;
                    if ( ( l = (*pMap)(nRow, nCol)) < 0 )  // also boundary
                        continue;
                    if ( l != l2 && l2 >= 0 )
                        break;

                    l2 = l;
                }

                if ( i >= 8 && l2 >= 0 )  // inside the same class
                    (*pMap)(r, c) = l2;
            }
        }
    py::array_t<int> aMap = GrayImage2PyArray(pMap);
    delete pMap;
    return aMap;
}

py::array_t<int> PyRAG::KmeansInitialization() {

    // corresponds to m_nInitMode == KM_REGION in GaussianRegionClustering.cpp
    TRegionBasedClustering* pInitRegionClustering = new TRegionBasedClustering(m_pGraph, m_nClassCnt);

    TGrayImage<int>* pRet = 0;

    double dbInitRatio = 0.000001;
    int nInit1 = 25;
    int nInit2 = 20;
    int nInitIterFinal = 500;
//    int nInitInitMode = 1; // TRegionClustering::RAND2
    int nInitInitMode = 1;
    bool bInitSucFlag;

    pRet = pInitRegionClustering->RegionKMClustering(nInitInitMode, bInitSucFlag, dbInitRatio, nInit1, nInit2, nInitIterFinal);

    if ( pRet )
    {
        int i;
        TVector** ppMu = pInitRegionClustering->GetMean();
        TMatrix** ppCovInv = pInitRegionClustering->GetCovInv();
        double* dbPrior = pInitRegionClustering->GetPrior();
        double* dbCovDet = pInitRegionClustering->GetCovDet();
        for ( i = 0; i < m_nClassCnt; i++ )
        {
            gmm->m_ppMu[i]->CopyFrom(ppMu[i]);
            gmm->m_ppCovInv[i]->CopyFrom(ppCovInv[i]);
            gmm->m_dbPrior[i] = dbPrior[i];
            gmm->m_dbCovDet[i] = dbCovDet[i];
            gmm->m_dbLnPrior[i] = log(dbPrior[i]);
            gmm->m_dbLnSigma[i] = log(dbCovDet[i]);
        }

        delete pRet; // this HAS to be deleted.
        delete pInitRegionClustering;
    }

    // FinalClassify function from GaussianRegionClustering.cpp
    double U, dbMin;
    int i, nMinIdx;
    int* nLabelFlag = new int[m_nClassCnt];
    memset(nLabelFlag, 0, m_nClassCnt * sizeof(int));
    TRegion* p = (TRegion*) m_pGraph->GetVertices();
    while (p)
    {
        dbMin = 999999999999;  // big enough
        for ( i = 0; i < m_nClassCnt; i++ )
        {
            if ( ( U = GetUnaryEnergy(p, i) ) < dbMin )
            {
                dbMin = U;
                nMinIdx = i;
            }
        }
        p->m_nClassLabel = nMinIdx;

        if ( nLabelFlag[nMinIdx] == 0 )
            nLabelFlag[nMinIdx] = 1;

        p = (TRegion*) p->m_pNext;
    }
    return GetRegionClassLabels();
}

py::array_t<int> PyRAG::MergeRegions() {

    m_pGraph->SortEdge();
    std::vector<std::pair<int,int>> MergedRegions = m_pGraph->MergeAndReturnIndex();
    int n_merged;

    if (MergedRegions.empty()) n_merged=0;
    else n_merged = MergedRegions.size();

    py::array_t<int> aMerged ({n_merged, 2});

    if (n_merged > 0){
        auto ra = aMerged.mutable_unchecked<2>();
        //py::print("NUMBER OF MERGED VERTICES:");
        //py::print(n_merged);
        //py::print("nice!!!");
        for(int i=0; i<n_merged; i++){
            ra(i,0) = MergedRegions[i].first;
            ra(i,1) = MergedRegions[i].second;
        }

        // update the consecutive mapping
        TEstimateUtil pEstUtil;
        int nRegions;
        int* pMapping = pEstUtil.GenerateConsecutiveMapping(nRegions, m_pGraph->GetMap(), m_pMask);
        for(int i=0; i<nRegions; i++){
            m_vConsecutiveMapping[i] = pMapping[i];
        }
        m_pBoundaryLabelMap->SetValue(-2);
        // update the boundary map
        int edgeId=0;
        TStrengthPenaltyBoundary* pCurBoundary = (TStrengthPenaltyBoundary*) m_pGraph->GetEdges();
        while (pCurBoundary)
        {

            TChainCode *pChain = pCurBoundary->m_pChain;
            POINTex Pt;
            int nDir;

            while (pChain)
            {
                pChain->GetFirstSite(Pt, nDir);
                do {
                    (*m_pBoundaryLabelMap)(Pt.Row, Pt.Col) = edgeId;

                } while (pChain->GetNextSite(Pt, nDir));

                pChain = pChain->m_pNext;
            }
            edgeId++;
            pCurBoundary = (TStrengthPenaltyBoundary*) pCurBoundary->m_pNext;
        }


    } else {
        //py::print("NO REGIONS MERGED");
    }

    return aMerged;
}


double PyRAG::GetUnaryEnergy(TRegion *pRegion, int nClassIndex) const
{
    TGaussianRegion* p = (TGaussianRegion*) pRegion;
    int nFeatureCnt = gmm->m_ppMu[nClassIndex]->GetElementCnt();
    double *dbMu = gmm->m_ppMu[nClassIndex]->m_ppData[0];
    double **dbCoMu = gmm->m_ppCovInv[nClassIndex]->m_ppData;
    double U = 0;
    for (int i = 0; i < nFeatureCnt; i ++)
    {
        U += ((*p->m_pCoMean)(i, i) - 2 * (*p->m_pMean)(i) * dbMu[i] +
              dbMu[i] * dbMu[i]) * dbCoMu[i][i];
        for (int j = i + 1; j < nFeatureCnt; j ++)
            U += ((*p->m_pCoMean)(i, j) - (*p->m_pMean)(i) * dbMu[j] -
                  (*p->m_pMean)(j) * dbMu[i] + dbMu[i] * dbMu[j]) *
                 dbCoMu[i][j] * 2;
    }
    U = U / 2 + gmm->m_dbLnSigma[nClassIndex];  //Note that the determinant of COV is squared in the program
    return ( U * p->m_nPixelCnt / nFeatureCnt );
}

double PyRAG::GetBinaryEnergy(int nLabel1, int nLabel2, TBoundary *pBoundary) const {
    double penalty=0;
    penalty = (nLabel1 == nLabel2)? 0.0 : ((TStrengthPenaltyBoundary*) pBoundary)->m_dbPenaltyEnergy;
    return penalty;
}