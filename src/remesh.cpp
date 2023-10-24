
// #ifdef USE_GEOGRAM
// #include <geogram/mesh/mesh.h>
// #include <geogram/mesh/mesh_fill_holes.h>
// #include <geogram/basic/command_line_args.h>
// #include <geogram/basic/command_line.h>
// #endif

#include "remesh.h"
#include "meshAlgorithm.h"
#include "igl/remove_unreferenced.h"

#include <vector>
#include <array>
#include <queue>
#include <set>

using namespace MESHIO::polymesh;



void MESHIO::split_long_edges(PolyMesh *mesh, double high)
{
  int NE = mesh->numEdges();
  for (int i = 0; i < NE; i++)
  {
    MEdge *e = mesh->edge(i);
    double len = e->length();
    if (len > high)
      mesh->splitEdgeTriangle(e);
  }
}

void MESHIO::collapse_short_edges(PolyMesh *mesh, double high, double low)
{
  int NE = mesh->numEdges();
  for (int i = NE - 1; i >= 0; i--)
  {
    if (i > mesh->numEdges() - 1)
      continue;
    MEdge *e = mesh->edge(i);
    MHalfedge *he = e->halfEdge();
    MVert *p0 = he->fromVertex();
    MVert *p1 = he->toVertex();
    if (!mesh->is_collapse_ok(he))
      continue;
    if (mesh->isBoundary(p0) || mesh->isBoundary(p1))
      continue;
    double len = e->length();
    if (len < low)
    {
      bool is_collapse = true;
      for (VertexVertexIter vv_it = mesh->vv_iter(p0); vv_it.isValid(); ++vv_it)
      {
        MVert *vv = *vv_it;
        double len = (p1->position() - vv->position()).norm();
        if (len > high)
        {
          is_collapse = false;
          break;
        }
      }
      if (is_collapse)
        mesh->collapseTriangle(he);
    }
  }
}

void MESHIO::delete_lowdegree(PolyMesh *mesh)
{
  int NE = mesh->numEdges();
  for (int i = NE - 1; i >= 0; i--)
  {
    if (i > mesh->numEdges() - 1)
      continue;

    MEdge *e = mesh->edge(i);
    MHalfedge *he = e->halfEdge();
    MVert *p0 = he->fromVertex();
    MVert *p1 = he->toVertex();

    if (
        (!mesh->isBoundary(e)) && (!mesh->isBoundary(p0)) && (!mesh->isBoundary(p1)) &&
        ((mesh->valence(p0) == 3) || (mesh->valence(p1) == 3)))
    {
      // std::cout << "delete non-boundary 3 degree point" << std::endl;
      // std::cout << "collapse edge p0-p1£º "
      //<< p0->index() << " " << p1->index() << std::endl;

      mesh->collapseTriangle(he);
    }
  }
}

void MESHIO::equalize_valences(PolyMesh *mesh)
{
  std::vector<int> target_valence;
  int deviation_pre, deviation_post;
  for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
  {
    if (mesh->isBoundary(*v_it))
      target_valence.push_back(4);
    else
      target_valence.push_back(6);
  }

  for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
  {
    if (mesh->isBoundary(*e_it) || !mesh->is_flip_ok_Triangle(*e_it))
      continue;

    MHalfedge *he1 = (*e_it)->halfEdge(); // edge 1->2
    MVert *v1 = he1->fromVertex();
    MVert *v2 = he1->toVertex();

    MHalfedge *he2 = (*e_it)->halfEdge()->next(); // edge 2->3
    MVert *v3 = he2->toVertex();
    MHalfedge *he3 = (*e_it)->halfEdge()->pair()->next(); // edge 1->4
    MVert *v4 = he3->toVertex();

    /*angle judge*/
    // *** notation(original configuration)
    //
    //                        1        4
    //                       *---------*
    //						/ \		  /
    //					   /   \  E2 /
    //					  /     \   /
    //				     /	E1	 \ /
    //				    *---------*
    //                  3         2
    //
    //      Diagonal swapping. '
    /**
     *        a
     *      /
     *    /
     * b ----- c
     * calculate the angle of (ba and bc)
     */
    auto vectorAngle = [](const Eigen::Vector3d &a,
                          const Eigen::Vector3d &b,
                          const Eigen::Vector3d &c)
    {
      Eigen::Vector3d ba = (a - b).normalized();
      Eigen::Vector3d bc = (c - b).normalized();
      double cos_theta = ba.dot(bc);
      double ang = std::acos(cos_theta) * 180.0 / M_PI;
      return ang;
    };

    double a123 = vectorAngle(v1->position(), v2->position(), v3->position());
    double a231 = vectorAngle(v2->position(), v3->position(), v1->position());
    double a142 = vectorAngle(v1->position(), v4->position(), v2->position());
    double a421 = vectorAngle(v4->position(), v2->position(), v1->position());
    double a312 = 180 - a123 - a231;
    double a214 = 180 - a142 - a421;
    std::vector<double> BeforeSwap{a123, a231, a142, a421, a312, a214};

    auto MinAngle_before = std::min_element(BeforeSwap.begin(), BeforeSwap.end());

    double a234 = vectorAngle(v2->position(), v3->position(), v4->position());
    double a423 = vectorAngle(v4->position(), v2->position(), v3->position());
    double a143 = vectorAngle(v1->position(), v4->position(), v3->position());
    double a314 = vectorAngle(v3->position(), v1->position(), v4->position());
    double a342 = 180 - a234 - a423;
    double a431 = 180 - a143 - a314;
    std::vector<double> AfterSwap{a234, a423, a143, a314, a342, a431};

    auto MinAngle_after = std::min_element(AfterSwap.begin(), AfterSwap.end());

    // if the min angle get smaller after flip, this operation should be banned!
    if ((*MinAngle_after) < (*MinAngle_before))
      continue;

    /*degree judge*/
    deviation_pre = abs(int(mesh->valence(v1) - target_valence[v1->index()])) + abs(int(mesh->valence(v2) - target_valence[v2->index()])) + abs(int(mesh->valence(v3) - target_valence[v3->index()])) + abs(int(mesh->valence(v4) - target_valence[v4->index()]));
    deviation_post = abs(int(mesh->valence(v1) - 1 - target_valence[v1->index()])) + abs(int(mesh->valence(v2) - 1 - target_valence[v2->index()])) + abs(int(mesh->valence(v3) + 1 - target_valence[v3->index()])) + abs(int(mesh->valence(v4) + 1 - target_valence[v4->index()]));

    // non-boundary point shouldn't have degree less than 3
    if ((!mesh->isBoundary(v1)) && (mesh->valence(v1) - 1 <= 3))
      continue;
    if ((!mesh->isBoundary(v2)) && (mesh->valence(v2) - 1 <= 3))
      continue;
    if ((!mesh->isBoundary(v3)) && (mesh->valence(v3) + 1 <= 3))
      continue;
    if ((!mesh->isBoundary(v4)) && (mesh->valence(v4) + 1 <= 3))
      continue;

    if (deviation_pre > deviation_post || (*MinAngle_before) < 1)
      mesh->flipEdgeTriangle(*e_it);
  }

  /*boundary judge*/
  // if two points are boundary but the edge is none boundary,then must flip!
  // for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
  //{
  //	if (mesh->isBoundary(*e_it) || !mesh->is_flip_ok_Triangle(*e_it)) continue;

  //	MHalfedge* he1 = (*e_it)->halfEdge(); //edge 1->2
  //	MVert* v1 = he1->fromVertex();
  //	MVert* v2 = he1->toVertex();

  //	MHalfedge* he2 = (*e_it)->halfEdge()->next();//edge 2->3
  //	MVert* v3 = he2->toVertex();
  //	MHalfedge* he3 = (*e_it)->halfEdge()->pair()->next();//edge 1->4
  //	MVert* v4 = he3->toVertex();

  //	if ((mesh->isBoundary(v1) && mesh->isBoundary(v2))
  //		&& (mesh->isBoundary(v3) && mesh->isBoundary(v4)))
  //		continue;

  //	if ((!mesh->isBoundary(*e_it)) && (mesh->isBoundary(v1) && mesh->isBoundary(v2)))

  //		mesh->flipEdgeTriangle(*e_it);
  //}
}

void MESHIO::tangential_relaxation(PolyMesh *mesh)
{
  mesh->updateMeshNormal();

  auto get_TriFace_Area = [](const Eigen::Vector3d &v0, const Eigen::Vector3d &v1, const Eigen::Vector3d &v2)
  {
    return sqrt(0.5 * ((v1 - v0).cross((v2 - v0))).squaredNorm());
  };

  for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
  {
    if (mesh->isBoundary(*v_it))
      continue;
    double count = 0.0;
    Eigen::Vector3d p = (*v_it)->position();
    Eigen::Vector3d q(0.0, 0.0, 0.0);

    // CPT
    double area = 0, sum_area = 0;
    for (VertexFaceIter vf_it = mesh->vf_iter(*v_it); vf_it.isValid(); ++vf_it)
    {
      std::vector<Eigen::Vector3d> vertexlist;
      for (FaceVertexIter fv_it = mesh->fv_iter(*vf_it); fv_it.isValid(); ++fv_it)
      {
        vertexlist.push_back((*fv_it)->position());
      }

      area = get_TriFace_Area(vertexlist[0], vertexlist[1], vertexlist[2]);
      sum_area += area;

      q += ((*vf_it)->getFaceCenter()) * area;
    }
    q /= sum_area;

    // Laplace
    // for (VertexVertexIter vv_it = mesh->vv_iter(*v_it); vv_it.isValid(); ++vv_it)
    //{
    //	q += (*vv_it)->position();
    //	++count;
    // }
    // q /= count;
    Eigen::Vector3d n = (*v_it)->normal();
    n.normalize();

    Eigen::Vector3d newPos = q + (n.dot(p - q)) * n;

    //// For continue surface.
    // Eigen::Vector2d res_uv = {(*v_it)->u(), (*v_it)->v()};
    // double dis;
    // surface->OrthogonalProjection(newPos, 1e-5, res_uv, dis);
    // surface->param_to_coord(res_uv, newPos);

    (*v_it)->setPosition(newPos);
  }
}

void MESHIO::get_AABB_tree(PolyMesh *mesh, AABB_Tree*& abtree)
{
  std::vector<Vector3f> point_set;
  point_set.clear();
  for (FaceIter f_it = mesh->polyfaces_begin(); f_it != mesh->polyfaces_end(); ++f_it)
  {
    for (FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.isValid(); ++fv_it)
    {
      MVert *fv = *fv_it;
      Vector3f p;
      p[0] = float(fv->x());
      p[1] = float(fv->y());
      p[2] = float(fv->z());
      point_set.push_back(p);
    }
  }
  abtree = new AABB_Tree(point_set);
}

void MESHIO::project_to_surface(PolyMesh *mesh, AABB_Tree *abtree)
{
  for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
  {
    Vector3f p;
    p[0] = float((*v_it)->x());
    p[1] = float((*v_it)->y());
    p[2] = float((*v_it)->z());
    Vector3f ab_nearst_point;
    abtree->findNearstPoint(p, ab_nearst_point);
    Eigen::Vector3d new_point;
    new_point[0] = double(ab_nearst_point[0]);
    new_point[1] = double(ab_nearst_point[1]);
    new_point[2] = double(ab_nearst_point[2]);
    (*v_it)->setPosition(new_point);
  }
}

double MESHIO::calculateTargetEdgeLength(PolyMesh *mesh)
{

  double target_edge_length = 0.0;
  std::cout << mesh->numEdges() << std::endl;
  for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
  {
    target_edge_length += (*e_it)->length();
  }
  target_edge_length /= mesh->numEdges();
  return target_edge_length;
}

void MESHIO::removeHangingFace(Mesh &mesh)
{
  std::cout << "trying to remove the hanging face" << std::endl;
  auto &T = mesh.Topo;
  auto &V = mesh.Vertex;

  std::vector<std::set<int>> point_to_triangle(V.rows());
  for (int i = 0; i < T.rows(); i++)
  {
    for (int j = 0; j < 3; j++)
    {
      point_to_triangle[T(i, j)].insert(i);
    }
  }

  auto judge_hanging = [&](int tri_id)
  {
    int edge_1 = 0;
    int edge_3 = 0;
    for (int j = 0; j < 3; ++j)
    {
      int v1 = T(tri_id, j);
      int v2 = T(tri_id, (j + 1) % 3);
      if (v1 > v2)
        std::swap(v1, v2);
      std::vector<int> edge_facet;
      std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(edge_facet));

      if (edge_facet.size() == 1)
        edge_1++;
      if (edge_facet.size() > 2)
        edge_3++;
    }
    if (edge_1 > 0 && edge_3 > 0)
    {
      return true;
    }
    return false;
  };

  std::queue<int> que;
  std::vector<int> removed_tri(T.rows(), 0);
  for (int i = 0; i < T.rows(); ++i)
  {
    if (judge_hanging(i))
    {
      que.push(i);
      removed_tri[i] = 1;
    }
  }

  while (!que.empty())
  {
    int cur_t = que.front();
    que.pop();
    for (int i = 0; i < 3; ++i)
    {
      point_to_triangle[T(cur_t, i)].erase(cur_t);
    }
    for (int i = 0; i < 3; ++i)
    {
      int v1 = T(cur_t, i);
      int v2 = T(cur_t, (i + 1) % 3);
      std::vector<int> nxt_tri;
      std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
      for (int j = 0; j < nxt_tri.size(); ++j)
      {
        if (removed_tri[nxt_tri[j]] == 0 && judge_hanging(nxt_tri[j]))
        {
          que.push(nxt_tri[j]);
          removed_tri[nxt_tri[j]] = 1;
        }
      }
    }
  }


  /////////////////////////////

  std::vector<std::array<int, 2>> non_manifold_edge;
  for (int i = 0; i < T.rows(); i++)
  {
    if (removed_tri[i])
      continue;
    for (int j = 0; j < 3; j++)
    {
      int v1 = T(i, j);
      int v2 = T(i, (j + 1) % 3);
      std::vector<int> nxt_tri;
      std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
      if(nxt_tri.size() != 2){
        non_manifold_edge.emplace_back(std::array<int, 2>{v1, v2});
      }
    }
  }

  std::vector<std::vector<int>> loops;
  get_boundary_loop(non_manifold_edge, loops);
  std::vector<int> warn_facet;
  std::vector<int> tmp_edge;
  for(int i = 0; i < loops.size(); ++i){
    int edge_1 = 0;
    int edge_3 = 0;
    for(int j = 0; j < loops[i].size(); ++j){
      int v1 = loops[i][j];
      int v2 = loops[i][(j + 1) % loops[i].size()];
      std::vector<int> nxt_tri;
      std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
      tmp_edge.insert(tmp_edge.end(), nxt_tri.begin(), nxt_tri.end());

      if(nxt_tri.size() == 1) edge_1 ++;
      if(nxt_tri.size() > 2) edge_3 ++;
    }
    if(edge_1 > 0 && edge_3 > 0){
      warn_facet.insert(warn_facet.end(), tmp_edge.begin(), tmp_edge.end());
    }
  }
  for(int i = 0; i < warn_facet.size(); ++i){
    removed_tri[warn_facet[i]] = true;
    for (int j = 0; j < 3; ++j)
    {
      point_to_triangle[T(warn_facet[i], j)].erase(warn_facet[i]);
    }
  }

  {
    std::vector<std::array<int, 2>> non_manifold_edge;
  for (int i = 0; i < T.rows(); i++)
  {
    if (removed_tri[i])
      continue;
    for (int j = 0; j < 3; j++)
    {
      int v1 = T(i, j);
      int v2 = T(i, (j + 1) % 3);
      std::vector<int> nxt_tri;
      std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
      if(nxt_tri.size() != 2){
        non_manifold_edge.emplace_back(std::array<int, 2>{v1, v2});
      }
    }
  }

  std::vector<std::vector<int>> loops;
  get_boundary_loop(non_manifold_edge, loops);
  std::vector<int> warn_facet;
  std::vector<int> tmp_edge;
  for(int i = 0; i < loops.size(); ++i){
    int edge_1 = 0;
    int edge_3 = 0;
    for(int j = 0; j < loops[i].size(); ++j){
      int v1 = loops[i][j];
      int v2 = loops[i][(j + 1) % loops[i].size()];
      std::vector<int> nxt_tri;
      std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
      tmp_edge.insert(tmp_edge.end(), nxt_tri.begin(), nxt_tri.end());

      if(nxt_tri.size() == 1) edge_1 ++;
      if(nxt_tri.size() > 2) edge_3 ++;
    }
    if(edge_1 > 0 && edge_3 > 0){
      warn_facet.insert(warn_facet.end(), tmp_edge.begin(), tmp_edge.end());
    }
  }
  for(int i = 0; i < warn_facet.size(); ++i){
    removed_tri[warn_facet[i]] = true;
    for (int j = 0; j < 3; ++j)
    {
      point_to_triangle[T(warn_facet[i], j)].erase(warn_facet[i]);
    }
  }
  }

  // {
  // std::vector<std::array<int, 2>> non_manifold_edge;
  // for (int i = 0; i < T.rows(); i++)
  // {
  //   if (removed_tri[i])
  //     continue;
  //   for (int j = 0; j < 3; j++)
  //   {
  //     int v1 = T(i, j);
  //     int v2 = T(i, (j + 1) % 3);
  //     std::vector<int> nxt_tri;
  //     std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
  //     if(nxt_tri.size() == 1){
  //       non_manifold_edge.emplace_back(std::array<int, 2>{v1, v2});
  //     }
  //   }
  // }

  // std::map<int, int> mp_du;
  // for(int i = 0; i < non_manifold_edge.size(); ++i){
  //   mp_du[non_manifold_edge[i][0]]++;
  //   mp_du[non_manifold_edge[i][1]]++;
  // }
  // for(auto& t : mp_du){
  //   if(t.second > 2){
  //     for(auto& t : point_to_triangle[t.first]){
  //       removed_tri[] = true;
  //     }
  //   }
  // }

  // std::vector<std::vector<int>> loops;
  // get_boundary_loop(non_manifold_edge, loops);
  // std::vector<int> warn_facet;
  // std::vector<int> tmp_edge;
  // for(int i = 0; i < loops.size(); ++i){
  //   int edge_1 = 0;
  //   int edge_3 = 0;
  //   for(int j = 0; j < loops[i].size(); ++j){
  //     int v1 = loops[i][j];
  //     int v2 = loops[i][(j + 1) % loops[i].size()];
  //     std::vector<int> nxt_tri;
  //     std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
  //     tmp_edge.insert(tmp_edge.end(), nxt_tri.begin(), nxt_tri.end());

  //     if(nxt_tri.size() == 1) edge_1 ++;
  //     if(nxt_tri.size() > 2) edge_3 ++;
  //   }
  //   if(edge_1 > 0 && edge_3 > 0){
  //     warn_facet.insert(warn_facet.end(), tmp_edge.begin(), tmp_edge.end());
  //   }
  // }
  // for(int i = 0; i < warn_facet.size(); ++i){
  //   removed_tri[warn_facet[i]] = true;
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     point_to_triangle[T(warn_facet[i], j)].erase(warn_facet[i]);
  //   }
  // }
  // }

  ///// for one body
  std::vector<int> tri_color(T.rows(), -1);
  int init_color = 0;
  for (int i = 0; i < T.rows(); ++i)
  {
    if (removed_tri[i] || tri_color[i] > -1)
      continue;
    std::queue<int> que;
    que.push(i);
    tri_color[i] = init_color;
    while (!que.empty())
    {
      int cur_t = que.front();
      que.pop();
      for (int i = 0; i < 3; ++i)
      {
        int v1 = T(cur_t, i);
        int v2 = T(cur_t, (i + 1) % 3);
        std::vector<int> nxt_tri;
        std::set_intersection(point_to_triangle[v1].begin(), point_to_triangle[v1].end(), point_to_triangle[v2].begin(), point_to_triangle[v2].end(), std::back_inserter(nxt_tri));
        for (int j = 0; j < nxt_tri.size(); ++j)
        {
          if (tri_color[nxt_tri[j]] == init_color)
            continue;
          que.push(nxt_tri[j]);
          tri_color[nxt_tri[j]] = init_color;
        }
      }
    }
    init_color++;
  }

  std::vector<int> color_num(init_color, 0);
  for (int i = 0; i < T.rows(); ++i)
  {
    if (removed_tri[i])
      continue;
    color_num[tri_color[i]]++;
  }
  int max_color = 0;
  int idx = 0;
  for (int i = 0; i < color_num.size(); ++i)
  {
    if (color_num[i] > max_color)
    {
      max_color = color_num[i];
      idx = i;
    }
  }
  for (int i = 0; i < T.rows(); ++i)
  {
    if (removed_tri[i])
      continue;
    if (tri_color[i] != idx)
      removed_tri[i] = true;
  }

  ////////// for body end.

  Eigen::MatrixXi newT;
  int newT_num = std::count(removed_tri.begin(), removed_tri.end(), 0);
  newT.resize(newT_num, 3);
  int count = 0;
  for (int i = 0; i < T.rows(); ++i)
  {
    if (removed_tri[i])
      continue;
    newT.row(count++) = T.row(i);
  }
  auto oldV = V;
  Eigen::VectorXi I;
  igl::remove_unreferenced(oldV, newT, V, T, I);


{
   std::map<std::pair<int, int>, std::set<int>>
        edge_to_facet;  // edge to facet
    for (int i = 0; i < T.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            int v1 = T(i, j);
            int v2 = T(i, (j + 1) % 3);
            if (v1 > v2) std::swap(v1, v2);
            edge_to_facet[std::make_pair(v1, v2)].insert(i);
        }
    }

    int open_edge_num = 0;
    std::vector<std::array<int, 2>> open_edge;
    std::map<std::pair<int, int>, bool> mp_open_edge;
    for (auto &t : edge_to_facet) {
        if (t.second.size() == 1) {
            open_edge.push_back({{t.first.first, t.first.second}});
            mp_open_edge[std::make_pair(t.first.first, t.first.second)] = 1;
            open_edge_num += 1;
        }
    }
    std::vector<std::vector<int>> loops;
    get_boundary_loop(open_edge, loops);
    for(int i = 0; i < loops.size(); ++i){
      for(int j = 0; j < loops[i].size(); ++j){
        int a = loops[i][j];
        int b = loops[i][(j + 1) % loops[i].size()];
        if(a > b) std::swap(a, b);
        mp_open_edge[std::make_pair(a, b)] = 0;
      }
    }

    if (loops.size() > 0) {

        std::vector<std::array<int, 2>> lines;
        for(int i = 0; i < loops.size(); ++i){
            for(int j = 0; j < loops[i].size(); ++j){
                int e1 = loops[i][j];
                int e2 = loops[i][(j + 1) % loops[i].size()];
                lines.emplace_back(std::array<int, 2>{e1, e2});
            }
        }
        writeLineVTK("./open_boundary.vtk", V, lines);
        writeLineVTK("./open_edge.vtk", V, open_edge);
    }
}
  // GEO::initialize();
  // GEO::Mesh geomesh;
  // geomesh.facets.create_triangles(T.rows());
  // geomesh.vertices.create_vertices(V.rows());
  // for (int i = 0; i < V.rows(); ++i)
  // {
  //   geomesh.vertices.point(i) = GEO::vec3(V(i, 0), V(i, 1), V(i, 2));
  // }

  // for (int i = 0; i < T.rows(); ++i)
  // {
  //   geomesh.facets.set_vertex(i, 0, T(i, 0));
  //   geomesh.facets.set_vertex(i, 1, T(i, 1));
  //   geomesh.facets.set_vertex(i, 2, T(i, 2));
  // }
  // geomesh.facets.connect();
  // // geomesh.facets.compute_borders();
  // GEO::CmdLine::import_arg_group("algo");
  // // GEO::CmdLine::set_arg("algo::hole_filling", "ear_cut");
  // // GEO::CmdLine::declare_arg(
  // //           "algo:hole_filling", "loop_split",
  // //           "Hole filling mode (loop_split, Nloop_split, ear_cut)");
  // GEO::fill_holes(geomesh, std::numeric_limits<double>::max());

  // size_t num_vertices = geomesh.vertices.nb();
  // size_t num_faces = geomesh.facets.nb();

  // V.resize(num_vertices, 3);
  // T.resize(num_faces, 3); // 假设是三角面

  // // 将顶点数据从geogram Mesh复制到Eigen矩阵
  // for (int i = 0; i < num_vertices; ++i)
  // {
  //   const GEO::vec3& vertex = geomesh.vertices.point(i);
  //   V(i, 0) = vertex.x;
  //   V(i, 1) = vertex.y;
  //   V(i, 2) = vertex.z;
  // }

  // // 将面数据从geogram Mesh复制到Eigen矩阵
  // for (int i = 0; i < num_faces; ++i)
  // {
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     T(i, j) = geomesh.facets.vertex(i, j);
  //   }
  // }
}

void MESHIO::remesh(Mesh& mesh){
  polymesh::PolyMesh half_mesh;
  for(int i = 0; i < mesh.Vertex.rows(); ++i){
    half_mesh.addVertex(mesh.Vertex.row(i), Eigen::Vector2d(0, 0));
  }
  for(int i = 0; i < mesh.Topo.rows(); ++i){
    auto t = std::vector<size_t>{size_t(mesh.Topo(i, 0)), size_t(mesh.Topo(i, 1)), size_t(mesh.Topo(i, 2))};
    half_mesh.addPolyFace(t);
  }
  AABB_Tree* abtree;
  get_AABB_tree(&half_mesh, abtree);
  double target_edge_length, high, low;
	target_edge_length = calculateTargetEdgeLength(&half_mesh)/2.0;
	high = 4.0 / 3.0 * target_edge_length;
	low = 4.0 / 5.0 * target_edge_length;
  for (int i = 0; i < 10; i++)
	{
		split_long_edges(&half_mesh, high);
		collapse_short_edges(&half_mesh, high, low);
		equalize_valences(&half_mesh);
		tangential_relaxation(&half_mesh);
		project_to_surface(&half_mesh, abtree);
	}

  mesh.Vertex.resize(half_mesh.numVertices(), 3);
  int v_id = 0, f_id = 0;
  for(VertexIter iter = half_mesh.vertices_begin(); iter != half_mesh.vertices_end(); iter++){
    mesh.Vertex.row(v_id++) = (*iter)->position();
  }
  mesh.Topo.resize(half_mesh.numPolygons(), 3);
  mesh.Masks.resize(half_mesh.numPolygons(), 1);
  for(FaceIter iter = half_mesh.polyfaces_begin(); iter != half_mesh.polyfaces_end(); iter++){
    std::vector<int> index;
		for (FaceVertexIter fv_it = half_mesh.fv_iter(*iter); fv_it.isValid(); ++fv_it)
		{
			MVert* fv = *fv_it;
			index.push_back(fv->index());
		}
    mesh.Topo.row(f_id) = Eigen::Vector3i(index[0], index[1], index[2]);
    mesh.Masks(f_id++, 0) = 0;
  }

}