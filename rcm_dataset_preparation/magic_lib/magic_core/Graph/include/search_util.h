//
// Created by max on 2022-01-31.
//

#include "RegionAdjacencyGraph.h"

#ifndef MAGIC_PY_FIND_BRIDGES_H
#define MAGIC_PY_FIND_BRIDGES_H

void DFS_cutvertex(TGraph* pGraph, std::vector<int>& found_itr, std::vector<int>& low, std::vector<int>& visited, std::vector<int>& cut_vertices, int u, int parent, int& itr){
    int children = 0;
    visited.at(u) = 1;
    found_itr.at(u) = low.at(u) = ++itr;

    TRegion* pU = ((TRegionAdjacencyGraph*) pGraph)->GetRegion(u);
    for (int i=0; i<pU->m_nNeighborCnt; i++){
        int v = pU->GetNeighborVertex(pU->m_pEdges[i])->m_nLabel;

        if (visited.at(v) == -1){
            children++;
            DFS_cutvertex(pGraph, found_itr, low, visited, cut_vertices, v, u, itr);
            low.at(u) = std::min(low.at(u), low.at(v));

            if(parent !=-1 && low.at(v) >= found_itr.at(u)) cut_vertices.push_back(u);

        } else if ( v!= parent) {
            low.at(u) = std::min(low.at(u), found_itr.at(v));
        }
    }

    if(parent !=-1 && children > 1) cut_vertices.push_back(u);
}

void DFS_cut(TGraph* pGraph, std::vector<int>& visit, std::vector<int>& cut, std::vector<int>& low, std::vector<int>& parent, int& count, TVertex* pV){

    int v = pV->m_nLabel;

    visit[v] = ++count;
    low[v] = count;

    for(int i=0; i<pV->m_nNeighborCnt; i++){
        TVertex* pW = pV->GetNeighborVertex(pV->m_pEdges[i]);
        int w = pW->m_nLabel;

        if (visit[w] == 0){
            parent[w] = v;
            DFS_cut(pGraph, visit, cut, low, parent, count, pW);

            low[v] = min(low[v], low[w]);

            if (low[w] >= visit[v]){
                cut[v] = 1;
            }
        } else if (parent[w] != v) {
            low[v] = min(low[v], visit[w]);
        }
    }
}

std::vector<int> find_cut_vertices(TGraph* pGraph){

    std::vector<int> visit;
    std::vector<int> cut;
    std::vector<int> low;
    std::vector<int> parent;

    int count=0;

    for(int i=0; i< pGraph->GetVertexCnt(); i++){
        visit.push_back(0);
        cut.push_back(0);
        low.push_back(0);
        parent.push_back(0);
    }

    for(TVertex* pV = pGraph->GetVertices(); pV; pV = pV->m_pNext){
        int i = pV->m_nLabel;
        if (visit[i] == 0){
            parent[i] = i;
            DFS_cut(pGraph, visit, cut, low, parent, count, pV);
        }
    }

    std::vector<int> cut_vertices;
    for(TVertex* pV = pGraph->GetVertices(); pV; pV = pV->m_pNext){
        int i = pV->m_nLabel;
        if (cut[i] == 1){
            cut_vertices.push_back(i);
        }
    }

    return cut_vertices;
}

void DFS_bridges(TGraph* pGraph, std::vector<int>& low, std::vector<int>& pre, int u, int v, int& cnt, std::vector<int>& bridges){
    pre.at(v) = cnt++;
    low.at(v) = pre.at(v);
    TRegion* pV = ((TRegionAdjacencyGraph*) pGraph)->GetRegion(v);
    for (int i=0; i<pV->m_nNeighborCnt; i++)
    {
        int w = pV->GetNeighborVertex(pV->m_pEdges[i])->m_nLabel;

        if (pre.at(w) == -1){

            DFS_bridges(pGraph, low, pre, v, w, cnt, bridges);
            low.at(v) = std::min(low.at(v), low.at(w));

            if(low.at(w) == pre.at(w)){
                //py::print("Found a bridge!!"); py::print(v); py::print(w);
                bridges.push_back(v); bridges.push_back(w);
            }
        } else if (w != u){
            low.at(v) = std::min(low.at(v), pre.at(w));
        }
    }

}

std::vector<int> find_bridges(TGraph* pGraph){
    int cnt=0;
    std::vector<int> low;
    std::vector<int> pre;
    std::vector<int> bridges;
    for(int i=0; i< pGraph->GetVertexCnt(); i++){
        low.push_back(-1);
        pre.push_back(-1);
    }
    for(int v=0; v< pGraph->GetVertexCnt(); v++){
        if(pre.at(v) == -1){
            DFS_bridges(pGraph, low, pre, v, v, cnt, bridges);
        }
    }

    return bridges;

}

// depth first search for finding connected clique centers
void DFS_cliques(TGrayImage<int>* ClusterCenterMap, int r, int c, int count)
{
    (*ClusterCenterMap)(r,c) = count;

    int rr,cc;

    // loop over 4-connected neighbours
    for (int i = 0; i < 4; i++) {

        rr = r + neighbor4[i].Row;
        cc = c + neighbor4[i].Col;

        if (rr < 0 || rr >= ClusterCenterMap->Height() ||
            cc < 0 || cc >= ClusterCenterMap->Width()) { // neighbour out of valid range
            continue;
        }

        if ((*ClusterCenterMap)(rr,cc) == -1){
            // recurse
            DFS_cliques(ClusterCenterMap, rr, cc, count);
        }
    }
}



#endif //MAGIC_PY_FIND_BRIDGES_H
