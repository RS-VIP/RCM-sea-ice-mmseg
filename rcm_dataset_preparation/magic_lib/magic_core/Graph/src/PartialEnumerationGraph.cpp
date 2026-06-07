//
// Created by max on 2021-12-14.
//

#include "PartialEnumerationGraph.h"
#include <algorithm>
#include "search_util.h"

///
/// TClusterNode implementation
///
TClusterNode::TClusterNode(TGraph* pGraph, int Label): TVertex(pGraph, Label)
{
    m_nCenterPoints = 0;
    m_nMemberRegions = 0;
    m_nMemberBoundaries = 0;
}

TClusterNode::~TClusterNode(){}

///
/// TClusterEdge Implementation
///
TClusterEdge::TClusterEdge(TGraph* pGraph, int nLabel1, int nLabel2, int nBoundaryIndex):
                TEdge(pGraph),  TAVLTreeNode(new TBoundaryKey(nLabel1, nLabel2)) {
    m_nLabel1 = nLabel1;
    m_nLabel2 = nLabel2;
    m_nBoundaryIndex = nBoundaryIndex;
}

TClusterEdge::~TClusterEdge() {}

void TClusterEdge::AddSharedRegion(int SharedRegionIdx) {
    if (m_vSharedRegions.empty()){
        m_vSharedRegions.push_back(SharedRegionIdx);
    }
    else if (*std::find(m_vSharedRegions.begin(), m_vSharedRegions.end(), SharedRegionIdx) != SharedRegionIdx){
        m_vSharedRegions.push_back(SharedRegionIdx);
    }
}

///
/// TPartialEnumerationGraph implementation
///
TPartialEnumerationGraph::TPartialEnumerationGraph(TRegionAdjacencyGraph *pRAG,
                                                   TGrayImage<int>* pWatershedMap,
                                                   TGrayImage<int>* pBoundaryLabelMap,
                                                   TGrayImage<int>* pClusterCenterMap,
                                                   TMonoImage* pMask): TGraph(false) {

    // TODO: get rid of this attribute in TRegionAdjacencyGraph, I'm not using it anymore
    if (pRAG->GetBoundaryMode() != TRegionAdjacencyGraph::NO_INTERSECTIONS)
        throw std::runtime_error("PartialEnumerationGraph: RAG mode must be constructed in NO_INTERSECTIONS mode.");

    m_pRAG = pRAG;
    m_pWatershedMap = pWatershedMap;
    m_pBoundaryLabelMap = pBoundaryLabelMap;
    m_pClusterCenterMap = pClusterCenterMap;
    m_pMask = pMask;
    m_pClusterLUT = nullptr;

    // TODO: fix the depth first search to find cut vertices in the RAG. I need this to make sure I don't miss any clique
    // constraint edges, but for some reason the current implementation segfaults on some (but not all) large images...?
    // ideas on why: non-consecutive region indices maybe?
    // anyway the there are few enough cut vertices that pretty much everything works fine without this but there are
    // a couple of inconsistencies in the clique labels.
    //    m_vRAGCutVertices = find_cut_vertices(m_pRAG);
}

TPartialEnumerationGraph::~TPartialEnumerationGraph() {
    delete[] m_pClusterLUT;
    m_EdgeList.m_pRoot = 0;
}


void TPartialEnumerationGraph::ConstructGraph()
{
    TEstimateUtil EstUtil;
    // make the labels consecutive. They should be already, but it's probably a good idea to make sure.
    int nVertexCnt = EstUtil.RelabelMap(m_pClusterCenterMap, m_pMask);
    if(nVertexCnt<1)
        return;

    std::vector<int> directions;

    TChainCode* pChain;
    POINTex Pt;
    int sz, offset=0, nDir;

    for (TBoundary* pB=(TBoundary*)m_pRAG->GetEdges(); pB; pB=(TBoundary*)pB->m_pNext)
    {
        pChain = pB->m_pChain;
        sz=0;

        m_vStartPoints.push_back(pChain->m_StartPt);
        m_vEndPoints.push_back(pChain->m_EndPt);

        pChain->GetFirstSite(Pt, nDir);
        while (pChain)
        {
            pChain->GetFirstSite(Pt, nDir);
            while (pChain->GetNextSite(Pt, nDir))
            {
                m_vChainCodes.push_back(nDir);
                if(nDir<0 | nDir > 3 ){
                    std::cerr << "Uh Oh... invalid chain code value " << nDir << std::endl;
                }
                sz++;
            }
            pChain = pChain->m_pNext;
        }
        m_vChainCodeLengths.push_back(sz);
        m_vChainCodeOffsets.push_back(offset);
        offset += sz;
    }
    // add the end point (n+1st bin edge)
    m_vChainCodeOffsets.push_back(offset);

    if (m_pClusterLUT) delete[] m_pClusterLUT;

    m_pClusterLUT = new TClusterNode*[nVertexCnt];

    int i;
    for (i = 0; i < nVertexCnt; i ++)
        m_pClusterLUT[i] = new TClusterNode(this, i);

    // load the cluster centers into the clusters
    LoadClusterCenters();

    int nWidth = m_pWatershedMap->Width();
    int nHeight = m_pWatershedMap->Height();

    TClusterNode* pCur;
    int r,c, wcls, bcls; POINTex bS, bE;
    for (i = 0; i < nVertexCnt; i ++)
    {
        pCur = m_pClusterLUT[i];
        for(int j=0; j < pCur->m_nCenterPoints; j++){
            Pt = pCur->m_vCenterPoints[j];

            for(auto k : neighbor8){ // loop over the 8-connected pixels around the current center point

                r = Pt.Row +k.Row; c=Pt.Col + k.Col;

                // make sure we're within the bounds of the image
                if ((r<0) || (c<0) || (r>nHeight-1) || (c>nWidth-1)) continue;

                // make sure we're not on a masked pixel
                if (m_pMask && !(*m_pMask)(r, c))
                    continue;

                // add the adjacent region to the list of regions in this cluster if it's not there already
                wcls = (*m_pWatershedMap)(r,c);
                if (wcls>=0 && !std::count(pCur->m_vMemberRegions.begin(), pCur->m_vMemberRegions.end(), wcls)){
                    pCur->m_vMemberRegions.push_back(wcls);
                    pCur->m_nMemberRegions++;
                }
            }

            for(int k=0; k<4; k++){
                r = Pt.Row + neighbor4[k].Row; c=Pt.Col + neighbor4[k].Col;

                // make sure we're within the bounds of the image
                if ((r<0) || (c<0) || (r>nHeight-1) || (c>nWidth-1)) continue;

                // make sure we're not on a masked pixel
                if (m_pMask && !(*m_pMask)(r, c))
                    continue;

                // add the adjacent boundary to the list of boundaries in this cluster if it's not there already
                bcls = (*m_pBoundaryLabelMap)(r,c);
                if (bcls>=0 && !std::count(pCur->m_vMemberBoundaries.begin(), pCur->m_vMemberBoundaries.end(), bcls)){

                    pCur->m_vMemberBoundaries.push_back(bcls);
                    pCur->m_vMemberBoundaryDirs.push_back(k);
                    pCur->m_nMemberBoundaries++;

                    bS = m_vStartPoints[bcls];
                    bE = m_vEndPoints[bcls];
                    if ((r==bS.Row && c==bS.Col) || (Pt.Row==bS.Row && Pt.Col==bS.Col)){
                        pCur->m_vMemberBoundaryHead.push_back(1);
                        pCur->m_vMemberBoundaryDisplacement.push_back(std::pair<int, int>(bE.Row-bS.Row, bE.Col-bS.Col));
                    } else {
                        pCur->m_vMemberBoundaryHead.push_back(0);
                        pCur->m_vMemberBoundaryDisplacement.push_back(std::pair<int, int>(bS.Row-bE.Row, bS.Col-bE.Col));
                    }
                }
            }

            // sort the member regions and member boundaries for each cluster
            std::sort(pCur->m_vMemberRegions.begin(), pCur->m_vMemberRegions.end());
//            std::sort(pCur->m_vMemberBoundaries.begin(), pCur->m_vMemberBoundaries.end());
        }
    }

    // Get the edges!!!
    std::vector<int> contains_cut_vertex;
    int rr,cc, l, b, nNeighbLabel, nNeighbCount;

    for (i = 0; i < nVertexCnt; i ++){

        pCur = m_pClusterLUT[i];
        int nCurLabel = pCur->m_nLabel;

        // wherever there is not a bridge in the RAG, clique edges have corresponding RAG boundaries
        // so we can just loop over all the RAG boundaries and create clique edges for all the ones
        // that have clique centers on either side. After this is done, we will have to treat the
        // places in the RAG where there is a bridge as special case.
        for (int j=0; j<pCur->m_nMemberBoundaries; j++){

            nNeighbCount = 0;
            nNeighbLabel = -1;

            TBoundary* pCurBoundary = (TBoundary*) m_pRAG->GetEdge(pCur->m_vMemberBoundaries[j]);
            TChainCode* pChain = pCurBoundary->m_pChain;

            // we find neighbour clusters by looking around the 4-connected neighbourhood at the
            // start and end points of the chain corresponding to the current boundary

            r = pChain->m_StartPt.Row; c = pChain->m_StartPt.Col;
            b = (*m_pBoundaryLabelMap)(r,c);

            for (int k=0; k<5; k++){
                rr = r + neighbor5[k].Row;
                cc = c + neighbor5[k].Col;
                if ((rr<0) || (cc<0) || (rr>nHeight-1) || (cc>nWidth-1)) continue;
                if (m_pMask && !(*m_pMask)(rr, cc)) continue;
                l = (*m_pClusterCenterMap)(rr,cc);
                if ((l >= 0) && (l!=nCurLabel)){
                    nNeighbCount++;
                    nNeighbLabel = l;
                }
            }

            r = pChain->m_EndPt.Row; c = pChain->m_EndPt.Col;
            for (int k=0; k<5; k++){
                rr = r + neighbor5[k].Row;
                cc = c + neighbor5[k].Col;
                if ((rr<0) || (cc<0) || (rr>nHeight-1) || (cc>nWidth-1)) continue;
                if (m_pMask && !(*m_pMask)(rr, cc)) continue;
                l = (*m_pClusterCenterMap)(rr,cc);
                if ((l >= 0) && (l!=nCurLabel)){
                    nNeighbCount++;
                    nNeighbLabel = l;
                }
            }

            if (nNeighbLabel>=0){ // found a neighbor region!! Add it to the list

                TClusterEdge* pEdge = IsEdgeInList(nCurLabel, nNeighbLabel);

                if (!pEdge){
                    pEdge = new TClusterEdge(this, nCurLabel, nNeighbLabel, b);
                    pEdge->AddSharedRegion(pCurBoundary->m_nLabel1);
                    pEdge->AddSharedRegion(pCurBoundary->m_nLabel2);
                    m_EdgeList.Insert(pEdge);
                }
                else{
                    pEdge->AddSharedRegion(pCurBoundary->m_nLabel1);
                    pEdge->AddSharedRegion(pCurBoundary->m_nLabel2);
                }
            }
        }

    }

    // inefficient, will have to fix later I think??
    for(int c=0; c<m_vRAGCutVertices.size(); c++){

        int cv = m_vRAGCutVertices.at(c);

        std::vector<int> cut_vertex_parents;

        for (int i = 0; i < nVertexCnt; i ++){
            pCur = m_pClusterLUT[i];
            for(int j=0; j<pCur->m_nMemberRegions; j++){
                if (pCur->m_vMemberRegions.at(j) == cv) cut_vertex_parents.push_back(pCur->m_nLabel);
            }
        }

        // add edges between all pairs of cliques that share a cut vertex
        for(int i=0; i< cut_vertex_parents.size(); i++){
            for(int j=i+1; j<cut_vertex_parents.size(); j++){

                int src = cut_vertex_parents.at(i);
                int dst = cut_vertex_parents.at(j);
                if (src != dst){
                    TClusterEdge* pEdge = IsEdgeInList(src, dst);

                    if (!pEdge){
                        pEdge = new TClusterEdge(this, src, dst, -1);
                        pEdge->AddSharedRegion(cv);
                        m_EdgeList.Insert(pEdge);
                    }
                    else{
                        pEdge->AddSharedRegion(cv);
                    }
                }

            }
        }
    }

    for (i = 0; i < nVertexCnt; i ++)
        AddVertex(m_pClusterLUT[i]);
}

//std::vector<int> TPartialEnumerationGraph::GetBoundaryPairChainCodes(){
//
//    TChainCode* pChain;
//    POINTex Pt;
//    int edgeId=0, nDir;
//    for (TBoundary* pB=(TBoundary*)m_pRAG->GetEdges(); pB; pB=(TBoundary*)pB->m_pNext)
//    {
//        pChain = pB->m_pChain;
//        pChain->GetFirstSite(Pt, nDir);
//        while (pChain)
//        {
//            pChain->GetFirstSite(Pt, nDir);
//            while (pChain->GetNextSite(Pt, nDir))
//            {
//            }
//            pChain = pChain->m_pNext;
//        }
//        edgeId++;
//    }
//    std::vector<int> PairChainCodes;
//    return PairChainCodes;
//}

void TPartialEnumerationGraph::LoadClusterCenters()
{
    int nWidth = m_pWatershedMap->Width();
    int nHeight = m_pWatershedMap->Height();

    // Load all the center pixels for each cluster
    int r, c, k;
    POINTex Pt;
    for (r = 0; r < nHeight; r ++) {
        for (c = 0; c < nWidth; c++) {
            if (m_pMask && !(*m_pMask)(r, c))
                continue;

            k = (*m_pClusterCenterMap)(r, c);
            if (k >= 0) {
                Pt.Row = r;
                Pt.Col = c;
                m_pClusterLUT[k]->AddCenterPoint(Pt);
            }
        }
    }
}


TVertex *TPartialEnumerationGraph::FirstNeighbor(TVertex *pVertex) {
    int nLabel = pVertex->m_nLabel;
    TBoundaryKey lKey(nLabel, 0);
    TBoundaryKey hKey(nLabel, nLabel-1);  // since segments with label greater than
    // current is not in graph yet, those segments are excluded.

    TClusterEdge* pCur = (TClusterEdge*) m_EdgeList.FindFirstInRange(lKey, hKey);
    if (pCur)
    {
        if (pCur->m_nLabel1 == nLabel)
            return m_pClusterLUT[pCur->m_nLabel2];
        else
            return m_pClusterLUT[pCur->m_nLabel1];
    }
    return 0;
}

//
// Get the next neighbor for constructing the graph
//
TVertex* TPartialEnumerationGraph::NextNeighbor(TVertex* pVertex)
{
    TClusterEdge* pCur = (TClusterEdge*) m_EdgeList.FindNextInRange();
    if (pCur)
    {
        if (pCur->m_nLabel1 == pVertex->m_nLabel)
            return m_pClusterLUT[pCur->m_nLabel2];
        else
            return m_pClusterLUT[pCur->m_nLabel1];
    }

    return 0;
}

TEdge* TPartialEnumerationGraph::ConstructEdge(TVertex* pNeighbor1, TVertex* pNeighbor2){
    // In the graph construction stage, use the existing elements
    // in the list directly
    TBoundaryKey Key(pNeighbor1->m_nLabel, pNeighbor2->m_nLabel);
    TAVLTreeNode* p = m_EdgeList.Find(Key);
    if (p)
    {
        ((TClusterEdge*) p)->SetVertexNeighbor(pNeighbor1, pNeighbor2);
        return (TClusterEdge*) p;
    }
    else
        return 0;
}

TClusterEdge* TPartialEnumerationGraph::IsEdgeInList(int nLabel1, int nLabel2)
{
    TBoundaryKey Key(nLabel1, nLabel2);
    TAVLTreeNode* p = m_EdgeList.Find(Key);
    if (p)
        return (TClusterEdge*) p;
    else
        return 0;
}

std::vector<int> TPartialEnumerationGraph::GetBoundaryPairChainCode(TClusterNode* pC, int b1, int b2){

    int nLabel = pC->m_nLabel;
    if (b1>=pC->m_nMemberBoundaries || b2>=pC->m_nMemberBoundaries){
        throw std::runtime_error("Boundary index is out of range for this clique!!");
    }
    std::vector<int> pair_chain_code;

    int pt1x, pt1y, pt2x, pt2y, dy, dx;
    int bLabel1 = pC->m_vMemberBoundaries[b1];
    int bLabel2 = pC->m_vMemberBoundaries[b2];

    if (pC->m_vMemberBoundaryHead[b1]>0){
        pt1x = m_vStartPoints[bLabel1].Col; pt1y = m_vStartPoints[bLabel1].Row;
        for (int i=m_vChainCodeOffsets[bLabel1+1] - 1; i>=m_vChainCodeOffsets[bLabel1]; i--){
            pair_chain_code.push_back((m_vChainCodes[i] + 2) & 3);
        }
    }else{
        pt1x = m_vEndPoints[bLabel1].Col; pt1y = m_vEndPoints[bLabel1].Row;
        for (int i=m_vChainCodeOffsets[bLabel1]; i<m_vChainCodeOffsets[bLabel1+1]; i++){
            pair_chain_code.push_back(m_vChainCodes[i]);
        }
    }

    if (pC->m_vMemberBoundaryHead[b2]>0){
        pt2x =m_vStartPoints[bLabel2].Col; pt2y =m_vStartPoints[bLabel2].Row;
    } else {
        pt2x =m_vEndPoints[bLabel2].Col; pt2y =m_vEndPoints[bLabel2].Row;
    }

    int rr,cc;

    int niter=0;
    while((pt1y!=pt2y) || (pt1x!=pt2x)){
        dy = pt2y - pt1y; dx = pt2x - pt1x;
        niter++;

        if (niter > 5){
            break;
        }

        if(dy > 0){
            rr = pt1y+1;cc=pt1x;
            //if ( !(m_pMask && !(*m_pMask)(rr, cc)) && (*m_pClusterCenterMap)(rr,cc) == nLabel)
            if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == nLabel)){
                pt1y=rr;
                pair_chain_code.push_back(3);
                continue;
            }
        }

        if(dx > 0){
            rr = pt1y;cc=pt1x+1;
            if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == nLabel)){
                pt1x=cc;
                pair_chain_code.push_back(2);
                continue;
            }
        }

        if(dy < 0){
            rr = pt1y-1;cc=pt1x;
            if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == nLabel)){
                pt1y=rr;
                pair_chain_code.push_back(1);
                continue;
            }
        }

        if(dx < 0){
            rr = pt1y;cc=pt1x-1;
            if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == nLabel)){
                pt1x=cc;
                pair_chain_code.push_back(0);
                continue;
            }
        }

        // ok so if we get to this point then it means that we couldn't reduce the inter-point distance in a single step, so
        // we just take any available step that puts us on top of the clique center
        for(int i=0; i<4; i++){
            rr = pt1y+neighbor4[i].Row; cc=pt1x+neighbor4[i].Col;
            if ((rr<0) || (cc<0) || (rr>m_pWatershedMap->Height()-1) || (cc>m_pWatershedMap->Width()-1)) continue;
            if((*m_pClusterCenterMap)(rr,cc) == nLabel){
                pt1x=cc;
                pt1y=rr;
                pair_chain_code.push_back(i);
                break;
            }
        }
    }

    if (pC->m_vMemberBoundaryHead[b2]==0){
        for (int i=m_vChainCodeOffsets[bLabel2+1] - 1; i>=m_vChainCodeOffsets[bLabel2]; i--){
            pair_chain_code.push_back((m_vChainCodes[i] + 2) & 3);
        }
    }else{
        for (int i=m_vChainCodeOffsets[bLabel2]; i<m_vChainCodeOffsets[bLabel2+1]; i++){
            pair_chain_code.push_back(m_vChainCodes[i]);
        }
    }

    return pair_chain_code;
}

void TPartialEnumerationGraph::AppendBoundaryPairChainCode(std::vector<int> &codes, std::vector<int> &pair_idx,
                                                           std::vector<int> &clique_idx, TClusterNode *pC, int b1,
                                                           int b2) {
    int pt1x, pt1y, pt2x, pt2y, dy, dx;
    int bLabel1 = pC->m_vMemberBoundaries[b1];
    int bLabel2 = pC->m_vMemberBoundaries[b2];
    int pair_idx_val = b1*pC->m_nMemberBoundaries + b2 + 1;
    int clique_idx_val = pC->m_nLabel;

    if (pC->m_vMemberBoundaryHead[b1]>0){
        pt1x = m_vStartPoints[bLabel1].Col; pt1y = m_vStartPoints[bLabel1].Row;
        for (int i=m_vChainCodeOffsets[bLabel1+1] - 1; i>=m_vChainCodeOffsets[bLabel1]; i--){
            codes.push_back((m_vChainCodes[i] + 2) & 3);
            pair_idx.push_back(pair_idx_val);
            clique_idx.push_back(clique_idx_val);
        }
    }else{
        pt1x = m_vEndPoints[bLabel1].Col; pt1y = m_vEndPoints[bLabel1].Row;
        for (int i=m_vChainCodeOffsets[bLabel1]; i<m_vChainCodeOffsets[bLabel1+1]; i++){
            codes.push_back(m_vChainCodes[i]);
            pair_idx.push_back(pair_idx_val);
            clique_idx.push_back(clique_idx_val);
        }
    }

    if (pC->m_vMemberBoundaryHead[b2]>0){
        pt2x =m_vStartPoints[bLabel2].Col; pt2y =m_vStartPoints[bLabel2].Row;
    } else {
        pt2x =m_vEndPoints[bLabel2].Col; pt2y =m_vEndPoints[bLabel2].Row;
    }

    int rr,cc;

    //std::cerr << pt1x << std::endl;
    //std::cerr << pt1y << std::endl;
    //std::cerr << pt2x << std::endl;
    //std::cerr << pt2y << std::endl;

    int niter=0, nWidth=m_pWatershedMap->Width(), nHeight=m_pWatershedMap->Height();
    while(pt1y!=pt2y || pt1x!=pt2x){

        dy = pt2y - pt1y;
        dx = pt2x - pt1x;
        niter++;

        if (niter > 5){
//            std::cerr << "THere are some problems going on over here." << std::endl;
//            std::cerr << "Clique Label: " << pC->m_nLabel << std::endl;
//            std::cerr << "Boundary Labels: " << bLabel1 << " , " << bLabel2 << std::endl;
//            std::cerr << dx << " , " << dy << std::endl;
            break;
        }

        if(dy > 0){
            rr=pt1y+1; cc=pt1x;

            //if ( !(m_pMask && !(*m_pMask)(rr, cc)) && (*m_pClusterCenterMap)(rr,cc) == nLabel)
            //if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val)){
            if ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val){
                pt1y=rr;
                codes.push_back(3);
                pair_idx.push_back(pair_idx_val);
                clique_idx.push_back(clique_idx_val);
                continue;
            }
        }

        if(dx > 0){
            rr = pt1y;cc=pt1x+1;
            //if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val)){
            if ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val){
                pt1x=cc;
                codes.push_back(2);
                pair_idx.push_back(pair_idx_val);
                clique_idx.push_back(clique_idx_val);
                continue;
            }
        }

        if(dy < 0){
            rr = pt1y-1;cc=pt1x;
            //if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val)){
            if ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val){
                pt1y=rr;
                codes.push_back(1);
                pair_idx.push_back(pair_idx_val);
                clique_idx.push_back(clique_idx_val);
                continue;
            }
        }

        if(dx < 0){
            rr = pt1y;cc=pt1x-1;
            //if ( !(m_pMask && !(*m_pMask)(rr, cc)) && ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val)){
            if ((*m_pBoundaryLabelMap)(rr,cc)==bLabel2 || (*m_pClusterCenterMap)(rr,cc) == clique_idx_val){
                pt1x=cc;
                codes.push_back(0);
                pair_idx.push_back(pair_idx_val);
                clique_idx.push_back(clique_idx_val);
                continue;
            }
        }

        // ok so if we get to this point then it means that we couldn't reduce the inter-point distance in a single step, so
        // we just take any available step that puts us on top of the clique center
        for(int i=0; i<4; i++){
            rr = pt1y+neighbor4[i].Row; cc=pt1x+neighbor4[i].Col;
            if ((rr<0) || (cc<0) || (rr>m_pWatershedMap->Height()-1) || (cc>m_pWatershedMap->Width()-1)) continue;
            if((*m_pClusterCenterMap)(rr,cc) == clique_idx_val){
                pt1x=cc;
                pt1y=rr;
                codes.push_back(i);
                pair_idx.push_back(pair_idx_val);
                clique_idx.push_back(clique_idx_val);
                // std::cerr << "fff " << pC->m_nLabel << " , " << b1 << " , " << b2 << std::endl;
                break;
            }
        }
    }

    if (pC->m_vMemberBoundaryHead[b2]==0){
        for (int i=m_vChainCodeOffsets[bLabel2+1] - 1; i>=m_vChainCodeOffsets[bLabel2]; i--){
            codes.push_back((m_vChainCodes[i] + 2) & 3);
            pair_idx.push_back(pair_idx_val);
            clique_idx.push_back(clique_idx_val);
        }
    }else{
        for (int i=m_vChainCodeOffsets[bLabel2]; i<m_vChainCodeOffsets[bLabel2+1]; i++){
            codes.push_back(m_vChainCodes[i]);
            pair_idx.push_back(pair_idx_val);
            clique_idx.push_back(clique_idx_val);
        }
    }
}

std::vector<int> TPartialEnumerationGraph::GetInnerBoundaryChainCode(TClusterNode *pC, int b) {

    if (b>pC->m_nMemberBoundaries){
        throw std::runtime_error("Boundary index is out of range for this clique!!");
    }
    std::vector<int> chain_code;

    int bLabel = pC->m_vMemberBoundaries[b];

    std::cerr << "Supposed length of chain code: " << m_vChainCodeLengths[bLabel] << std::endl;

    if (pC->m_vMemberBoundaryHead[b]>0){
        std::cerr << "backward boundary, flipping direction." << std::endl;
        std::cerr << "length: " << m_vChainCodeOffsets[bLabel+1] - m_vChainCodeOffsets[bLabel] << std::endl;
        for (int i=m_vChainCodeOffsets[bLabel+1] - 1; i>=m_vChainCodeOffsets[bLabel]; i--){
            chain_code.push_back((m_vChainCodes[i] + 2) & 3);
        }
    }else{
        std::cerr << "forward boundary, maintaining orientation." << std::endl;
        std::cerr << "length: " << m_vChainCodeOffsets[bLabel+1] - m_vChainCodeOffsets[bLabel] << std::endl;
        for (int i=m_vChainCodeOffsets[bLabel]; i<m_vChainCodeOffsets[bLabel+1]; i++){
            chain_code.push_back(m_vChainCodes[i]);
        }
    }

    return chain_code;
}

