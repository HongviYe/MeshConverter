
#include "remesh.h"

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

// void main(int argc, char** argv)
// {
// 	if (argc != 3)
// 	{
// 		std::cout << "========== Hw11 Usage  ==========\n";
// 		std::cout << std::endl;
// 		std::cout << "Input:	ACAM_mesh_HW11.exe	input_mesh.obj	output_mesh.obj\n";
// 		std::cout << std::endl;
// 		std::cout << "=================================================\n";
// 		return;
// 	}

// 	mesh = new PolyMesh();
// 	//¶ÁÈëÍø¸ñ
// 	std::string mesh_path = argv[1];
// 	loadMesh(mesh_path, mesh);

// 	std::string out_path = argv[2];

// 	//mesh load an write , now only support obj/off
// 	//loadMesh("small_bunny.obj", mesh);
// 	double target_edge_length, high, low;
// 	target_edge_length = calculateTargetEdgeLength()/2.0;
// 	high = 4.0 / 3.0 * target_edge_length;
// 	low = 4.0 / 5.0 * target_edge_length;
// 	get_AABB_tree();
// 	for (int i = 0; i < 10; i++)
// 	{
// 		split_long_edges(high);
// 		collapse_short_edges(high, low);
// 		equalize_valences();
// 		tangential_relaxation();
// 		project_to_surface();
// 	}
// 	writeMesh(out_path, mesh);

// }
