#include "triMesh.h"

using namespace std;

sfMesh::sfMesh(std::vector<std::vector<double>> &plist, std::vector<std::vector<int>> &flist) {
    this->dataCopy(plist, flist);
    this->Initialize();
}

sfMesh::sfMesh(Eigen::MatrixXd &plist, Eigen::MatrixXi &flist) {
    this->dataCopy(plist, flist);
    this->Initialize();
}

void sfMesh::dataCopy(std::vector<std::vector<double>> &plist, std::vector<std::vector<int>> &flist) {
    this->points.clear();
    this->facets.clear();
    for(auto & i : plist) {
        point p(i[0], i[1], i[2]);
        p.id = 2;
        this->points.push_back(p);
    }
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
}

void sfMesh::dataCopy(Eigen::MatrixXd &plist, Eigen::MatrixXi &flist) {
    this->points.clear();
    this->facets.clear();
    for (int i = 0; i < plist.rows(); i++) {
        point p(plist(i, 0), plist(i, 1), plist(i, 2));
        p.id = 2;
        this->points.push_back(p);
    }
    for (int i = 0; i < flist.size(); i++) {
        facet t;
        t.id = i;
        t.form[0] = flist(i, 0);
        t.form[1] = flist(i, 1);
        t.form[2] = flist(i, 2);
        facets.push_back(t);
        conn_v_t[t.form[0]].insert(i);
        conn_v_t[t.form[1]].insert(i);
        conn_v_t[t.form[2]].insert(i);
    }
}

void sfMesh::Initialize() {
    Block facet2block;
    facet2block.Init(facets.size());
    for(int i = 0; i < facets.size(); i++) {
        for(int j = 0; j < 3; j++) {
            int v1 = facets[i].form[(j + 1) % 3];
            int v2 = facets[i].form[(j + 2) % 3];
            set<int> &v1_set = conn_v_t[v1];
            set<int> &v2_set = conn_v_t[v2];
            set<int> inct;

            bool is_current_edge_manifold = true;
            set_intersection(std::begin(v1_set), std::end(v1_set), std::begin(v2_set), std::end(v2_set), std::inserter(inct, std::begin(inct)));
            if(inct.size() != 2) {
                this->isManifold = false;
                is_current_edge_manifold = false;
            }
            int nb_tri = i;
            for(auto f : inct) {
                if(f != i) {
                    nb_tri = f;
                    break;
                }
            }
            if (is_current_edge_manifold) {
                facets[i].neig[j] = nb_tri;

                /* yhf: should this line be outside "if"? Should we distinguish blocks connected by non-manifold edge?  */
                facet2block.Union(i, nb_tri);
            }
            else {
                facets[i].neig[j] = NON_MANIFOLD_EDGE_TAG;
            }
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
}
