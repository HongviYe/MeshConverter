#ifndef TIGER_TRI_MESH_H
#define TIGER_TRI_MESH_H

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>

#define NON_MANIFOLD_EDGE_TAG -1

class point {
public:
    int id = 0;
    Eigen::Vector3d coord;
    point(double x, double y, double z) : coord(Eigen::Vector3d(x, y, z)) {}
};

class facet {
public:
    int id = 0;
    int blockId = -1;
    std::array<int, 3> form = {0, 0, 0};
    std::array<int, 3> neig = {0, 0, 0};
};

class sfMesh {
public:
    int nBlock = 0;
    std::vector<point> points;
    std::vector<facet> facets;
    std::unordered_map<int, std::set<int>> conn_v_t;
    bool isManifold = true;
    sfMesh(std::vector<std::vector<double>> &plist, std::vector<std::vector<int>> &flist);
    sfMesh(Eigen::MatrixXd &plist, Eigen::MatrixXi &flist);
    void dataCopy(std::vector<std::vector<double>> &plist, std::vector<std::vector<int>> &flist);
    void dataCopy(Eigen::MatrixXd &plist, Eigen::MatrixXi &flist);
    void Initialize();
};

class Block {
public:
    int nSet;
    std::vector<int> parents;
    std::vector<int> rank;
    inline void Init(int n) {
        this->nSet = n;
        this->rank.resize(n);
        for(int i = 0; i < n; i++) {
            parents.push_back(i);
        }
    }
    inline int Find(int x) {
        return (x == parents[x]) ? x : Find(parents[x]);
    }
    inline void Union(int x1, int x2) {
        int f1 = Find(x1);
        int f2 = Find(x2);
        if(f1 == f2) return;
        if(rank[f1] < rank[f2])
            parents[f1] = f2;
        else if(rank[f1] > rank[f2])
            parents[f2] = f1;
        else {
            parents[f2] = f1;
            rank[f1]++;
        }
        nSet--;
    }
};

#endif  // TIGER_TRI_MESH_H