#ifndef TIGER_TRI_MESH_H
#define TIGER_TRI_MESH_H

#include <Eigen/Dense>
#include <vector>
#include <array>

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
    std::array<int, 3> form;
    std::array<int, 3> neig;
};

class sfMesh {
public:
    int nBlock;
    std::vector<point> points;
    std::vector<facet> facets;
    bool isManifold;
    sfMesh(std::vector<std::vector<double>> plist, std::vector<std::vector<int>> flist);
    void Init(std::vector<std::vector<double>> plist, std::vector<std::vector<int>> flist);
    void resetOrientation();
    void resetBlockOrientation(int start);
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
        return;
    }
};

#endif  // TIGER_TRI_MESH_H