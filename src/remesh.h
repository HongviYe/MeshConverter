#ifndef REMESHALGORITHM_H
#define REMESHALGORITHM_H

#include "Halfedge/AABB_Tree.h"
#include "Halfedge/PolyMesh.h"
#include "Halfedge/PolyMesh_Base.h"
#include <Eigen/Dense>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include "meshIO.h"

using namespace MESHIO::polymesh;

namespace MESHIO
{
  void split_long_edges(PolyMesh *mesh, double high);
  void collapse_short_edges(PolyMesh *mesh, double high, double low);
  void delete_lowdegree(PolyMesh *mesh);
  void equalize_valences(PolyMesh *mesh);
  void tangential_relaxation(PolyMesh *mesh);
  void project_to_surface(PolyMesh *mesh, AABB_Tree *abtree);
  void get_AABB_tree(PolyMesh *mesh, AABB_Tree *&abtree);
  double calculateTargetEdgeLength(PolyMesh *mesh);

  void removeHangingFace(Mesh& mesh);

  void remesh(Mesh& mesh);
}

#endif 