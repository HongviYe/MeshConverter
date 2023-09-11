#include "MeshOrient.h"
#include "triMesh.h"

#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <igl/winding_number.h>

using namespace std;

void resetOrientation(sfMesh& mesh);
void resetOrientationBFS(sfMesh& mesh, int start);
void resetOrientationVisualDriven(sfMesh& mesh);
void InOutFiltering(sfMesh& mesh);

int TIGER::resetOrientation(vector<vector<double>>& point_list, vector<vector<int>>& facet_list, vector<int>& block_mark) {
	sfMesh mesh(point_list, facet_list);
	if (!mesh.isManifold) {
		cout << "Non-manifold edge found! ignore corresponding edges." << endl;
	}
	resetOrientation(mesh);
	block_mark.resize(facet_list.size());
	for (int i = 0; i < facet_list.size(); i++) {
		for (int j = 0; j < 3; j++)
			facet_list[i][j] = mesh.facets[i].form[j];
		block_mark[i] = mesh.facets[i].blockId;
	}
	return 1;
}

int TIGER::resetOrinetation(Eigen::MatrixXd& point_list, Eigen::MatrixXi& facet_list, Eigen::VectorXi& block_mark) {
	sfMesh mesh(point_list, facet_list);
	if (!mesh.isManifold) {
		cout << "Non-manifold edge found! ignore corresponding edges." << endl;
	}
	resetOrientation(mesh);
	block_mark.resize(facet_list.size());
	for (int i = 0; i < facet_list.size(); i++) {
		for (int j = 0; j < 3; j++)
			facet_list(i, j) = mesh.facets[i].form[j];
		block_mark(i) = mesh.facets[i].blockId;
	}
	return 1;
}

void resetOrientation(sfMesh& mesh) {
	vector<int> blockStart;
	for (int bcnt = 0; bcnt < mesh.nBlock; bcnt++) {
		for (int i = 0; i < mesh.facets.size(); i++) {
			if (mesh.facets[i].blockId == bcnt) {
				blockStart.push_back(i);
				break;
			}
		}
	}
	for (int start : blockStart)
		resetOrientationBFS(mesh, start);

	InOutFiltering(mesh);
}

void resetOrientationBFS(sfMesh& mesh, int start) {
	queue<int> Q;
	Q.push(start);
	unordered_set<int> elemVisited;
	elemVisited.insert(start);
	vector<int> blockFacets;
	// BFS all the facet and reorient the facets by manifold criterion
	while (!Q.empty()) {
		facet& curFacet = mesh.facets[Q.front()];
		blockFacets.push_back(curFacet.id);
		// curVerts : vertex global Id - vertex local Id
		unordered_map<int, int> curVerts;
		for (int i = 0; i < 3; i++)
			curVerts[curFacet.form[i]] = i;
		for (int i = 0; i < 3; i++) {
			// Current neighbor element.
			if (curFacet.neig[i] == NON_MANIFOLD_EDGE_TAG) {
				continue;
			}
			if (elemVisited.find(curFacet.neig[i]) != elemVisited.end())
				continue;
			facet& nghFacet = mesh.facets[curFacet.neig[i]];
			elemVisited.insert(curFacet.neig[i]);
			int tmpV = 0;
			while (curVerts.find(nghFacet.form[tmpV]) != curVerts.end())
				tmpV++;
			int v1 = nghFacet.form[(tmpV + 1) % 3]; // global
			int v2 = nghFacet.form[(tmpV + 2) % 3]; // global
			int localCurV1 = curVerts[v1];
			int localCurV2 = curVerts[v2];
			if ((localCurV1 + 1) % 3 == localCurV2) // Neighbor element's normal direction need to be reset.
				swap(nghFacet.form[(tmpV + 1) % 3], nghFacet.form[(tmpV + 2) % 3]);
			Q.push(curFacet.neig[i]);
		}
		Q.pop();
	}
	// calculated volume
	double volume = 0;
	for (auto tri : blockFacets) {
		point v1 = mesh.points[mesh.facets[tri].form[0]];
		point v2 = mesh.points[mesh.facets[tri].form[1]];
		point v3 = mesh.points[mesh.facets[tri].form[2]];
		Eigen::Vector3d e21 = v2.coord - v1.coord;
		Eigen::Vector3d e31 = v3.coord - v1.coord;
		double tVol = e21.cross(e31).dot(v1.coord);
		volume += tVol;
	}

	if (volume < 0) {
		for (auto tri : blockFacets) {
			swap(mesh.facets[tri].form[0], mesh.facets[tri].form[1]);
		}
	}
}

void InOutFiltering(sfMesh& mesh) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd O;
	V.resize(mesh.points.size(), 3);
	F.resize(mesh.facets.size(), 3);
	O.resize(mesh.points.size(), 3);
	for (int i = 0; i < mesh.points.size(); i++) {
		V(i, 0) = mesh.points[i].coord(0);
		V(i, 1) = mesh.points[i].coord(1);
		V(i, 2) = mesh.points[i].coord(2);
		O(i, 0) = mesh.points[i].coord(0);
		O(i, 1) = mesh.points[i].coord(1);
		O(i, 2) = mesh.points[i].coord(2);
	}
	for (int i = 0; i < mesh.facets.size(); i++) {
		F(i, 0) = mesh.facets[i].form[0];
		F(i, 1) = mesh.facets[i].form[1];
		F(i, 2) = mesh.facets[i].form[2];
	}
	Eigen::VectorXd W;

	igl::winding_number(V, F, O, W);

	vector<double> block_winding_number_avg(mesh.nBlock, 0);
	vector<int> block_point_count(mesh.nBlock, 0);
	vector<int> visited(mesh.points.size(), 0);

	for (int i = 0; i < mesh.facets.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int pId = mesh.facets[i].form[j];
			int bId = mesh.facets[i].blockId;
			if (visited[pId] == 0) {
				block_winding_number_avg[bId] += W(pId);
				block_point_count[bId]++;
				visited[pId] = 1;
			}
		}
	}

	unordered_set<int> to_orient;
	for (int i = 0; i < mesh.nBlock; i++) {
		if (block_point_count[i] > 0) {
			block_winding_number_avg[i] /= block_point_count[i];
		}
		if (block_winding_number_avg[i] > 0.5) {
			to_orient.insert(i);
		}
	}

	for (int i = 0; i < mesh.facets.size(); i++) {
		if (to_orient.count(mesh.facets[i].blockId)) {
			swap(mesh.facets[i].form[0], mesh.facets[i].form[1]);
		}
	}
}