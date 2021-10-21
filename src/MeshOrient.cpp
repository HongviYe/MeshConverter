#include "MeshOrient.h"
#include "triMesh.h"

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <queue>

using namespace std;
using namespace MESHIO;

int MESHIO::resetOrientation(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &M){

	vector<vector<double>> point_list(V.rows(), vector<double>(V.cols()));
	vector<vector<int>> facet_list(F.rows(), vector<int>(F.cols()));
	vector<int> blockMark(M.rows());

	for (int i = 0; i < V.rows(); i++) {
		for (int j = 0; j < V.cols(); j++) {
			point_list[i][j] = V(i,j);
		}
	}
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			facet_list[i][j] = F(i, j);
		}
	}
	for (int i = 0; i < M.rows(); i++) {
		blockMark[i] = M(i,0);
	}

	sfMesh mesh(point_list, facet_list);
	if (!mesh.isManifold) {
		cout << "Reset Orientation failed. Input mesh is non-manifold." << endl;
		return 0;
	}
	mesh.resetOrientation();
	blockMark.resize(facet_list.size());
	for (int i = 0; i < facet_list.size(); i++) {
		for (int j = 0; j < 3; j++)
			facet_list[i][j] = mesh.facets[i].form[j];
		blockMark[i] = mesh.facets[i].blockId;
	}
	for (int i = 0; i < V.rows(); i++) {
		for (int j = 0; j < V.cols(); j++) {
			V(i, j)= point_list[i][j] ;
		}
	}
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			 F(i, j)= facet_list[i][j];
		}
	}
	for (int i = 0; i < M.rows(); i++) {
		M(i, 0) = blockMark[i];
	}

	return 1;
}


sfMesh::sfMesh(std::vector<std::vector<double>> plist, std::vector<std::vector<int>> flist) {
    this->Init(plist, flist);
    return;
}

void sfMesh::Init(std::vector<std::vector<double>> plist, std::vector<std::vector<int>> flist) {
    this->points.clear();
    this->facets.clear();
    for(int i = 0; i < plist.size(); i++) {
        point p(plist[i][0], plist[i][1], plist[i][2]);
        p.id = 2;
        this->points.push_back(p);
    }
    unordered_map<int, set<int>> conn_v_t;
    for(int i = 0; i < flist.size(); i++) {
        facet t;
        t.id = i;
        t.form[0] = flist[i][0];
        t.form[1] = flist[i][1];
        t.form[2] = flist[i][2];
        facets.push_back(t);
        conn_v_t[t.form[0]].insert(i);
        conn_v_t[t.form[1]].insert(i);
        conn_v_t[t.form[2]].insert(i);
    }

    this->isManifold = true;
    Block facet2block;
    facet2block.Init(facets.size());
    for(int i = 0; i < facets.size(); i++) {
        for(int j = 0; j < 3; j++) {
            int v1 = facets[i].form[(j + 1) % 3];
            int v2 = facets[i].form[(j + 2) % 3];
            set<int> &v1_set = conn_v_t[v1];
            set<int> &v2_set = conn_v_t[v2];
            set<int> inct;
            set_intersection(std::begin(v1_set), std::end(v1_set), std::begin(v2_set), std::end(v2_set), std::inserter(inct, std::begin(inct)));
            if(inct.size() != 2) {
                this->isManifold = false;
                return;
            }
            int nb_tri = i;
            for(auto f : inct) {
                if(f != i) {
                    nb_tri = f;
                    break;
                }
            }
            facets[i].neig[j] = nb_tri;
            facet2block.Union(i, nb_tri);
        }
    }

    unordered_map<int, int> tmpblocks;
    int bcnt = 0;
    for(int i = 0; i < facets.size(); ++i) {
        int parent = facet2block.Find(i);
        if(tmpblocks.find(parent) == tmpblocks.end())
            tmpblocks[parent] = bcnt++;
        facets[i].blockId = tmpblocks[parent];
    }
    this->nBlock = bcnt;
    return;
}

void sfMesh::resetOrientation() {
    vector<int> blockStart;
    for(int bcnt = 0; bcnt < this->nBlock; bcnt++) {
        for(int i = 0; i < facets.size(); i++) {
            if(facets[i].blockId == bcnt) {
                blockStart.push_back(i);
                break;
            }
        }
    }
    for(int start : blockStart)
        resetBlockOrientation(start);
    return;
}

void sfMesh::resetBlockOrientation(int start) {
    queue<int> Q;
    Q.push(start);
    unordered_set<int> elemVisited;
    elemVisited.insert(start);
    vector<int> blockFacets;
    while(!Q.empty()) {
        facet &curFacet = facets[Q.front()];
        blockFacets.push_back(curFacet.id);
        // curVerts : vertex global Id - vertex local Id
        unordered_map<int, int> curVerts;
        for(int i = 0; i < 3; i++)
            curVerts[curFacet.form[i]] = i;
        for(int i = 0; i < 3; i++) {
            // Current neighbor element.
            if(elemVisited.find(curFacet.neig[i]) != elemVisited.end())
                continue;
            facet &nghFacet = facets[curFacet.neig[i]];
            elemVisited.insert(curFacet.neig[i]);
            int tmpV = 0;
            while(curVerts.find(nghFacet.form[tmpV]) != curVerts.end())
                tmpV++;
            int v1 = nghFacet.form[(tmpV + 1) % 3]; // global
            int v2 = nghFacet.form[(tmpV + 2) % 3]; // global
            int localCurV1 = curVerts[v1];
            int localCurV2 = curVerts[v2];
            if((localCurV1 + 1) % 3 == localCurV2) // Neighbor element's normal direction need to be reset.
                swap(nghFacet.form[(tmpV + 1) % 3], nghFacet.form[(tmpV + 2) % 3]);
            Q.push(curFacet.neig[i]);
        }
        Q.pop();
    }
    // calculated volume
    double volume = 0;
    for(auto tri : blockFacets) {
        point v1 = points[facets[tri].form[0]];
        point v2 = points[facets[tri].form[1]];
        point v3 = points[facets[tri].form[2]];
        Eigen::Vector3d e21 = v2.coord - v1.coord;
        Eigen::Vector3d e31 = v3.coord - v1.coord;
        double tVol = e21.cross(e31).dot(v1.coord);
        volume += tVol;
    }

    if(volume < 0) {
        for(auto tri : blockFacets) {
            swap(facets[tri].form[0], facets[tri].form[1]);
        }
    }

    return;
}

