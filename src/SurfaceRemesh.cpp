
#include "remesh.h"
#include "SurfaceRemesh.h"
#include <vector>
#include <array>
#include <queue>
#include <set>




void SurfaceRemesh::iso_split_long_edges(PolyMesh* mesh, double high)
{
	int NE = mesh->numEdges();
	for (int i = 0; i < NE; i++)
	{
		MEdge* e = mesh->edge(i);
		MHalfedge* he = e->halfEdge();
		MVert* p0 = he->fromVertex();
		MVert* p1 = he->toVertex();
		if (mesh->isBoundary(p0) && mesh->isBoundary(p1)) continue;
		double len = e->length();
		if (len > high)
			mesh->splitEdgeTriangle(e);
	}
}

void SurfaceRemesh::iso_collapse_short_edges(PolyMesh* mesh, double high, double low)
{
	int NE = mesh->numEdges();
	for (int i = NE - 1; i >= 0; i--)
	{
		if (i > mesh->numEdges() - 1)
			continue;
		MEdge* e = mesh->edge(i);
		MHalfedge* he = e->halfEdge();
		MVert* p0 = he->fromVertex();
		MVert* p1 = he->toVertex();
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
				MVert* vv = *vv_it;
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


void SurfaceRemesh::iso_equalize_valences(PolyMesh* mesh)
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

		MHalfedge* he1 = (*e_it)->halfEdge(); // edge 1->2
		MVert* v1 = he1->fromVertex();
		MVert* v2 = he1->toVertex();

		MHalfedge* he2 = (*e_it)->halfEdge()->next(); // edge 2->3
		MVert* v3 = he2->toVertex();
		MHalfedge* he3 = (*e_it)->halfEdge()->pair()->next(); // edge 1->4
		MVert* v4 = he3->toVertex();

		/*degree judge*/
		deviation_pre = abs(int(mesh->valence(v1) - target_valence[v1->index()])) + abs(int(mesh->valence(v2) - target_valence[v2->index()])) + abs(int(mesh->valence(v3) - target_valence[v3->index()])) + abs(int(mesh->valence(v4) - target_valence[v4->index()]));
		deviation_post = abs(int(mesh->valence(v1) - 1 - target_valence[v1->index()])) + abs(int(mesh->valence(v2) - 1 - target_valence[v2->index()])) + abs(int(mesh->valence(v3) + 1 - target_valence[v3->index()])) + abs(int(mesh->valence(v4) + 1 - target_valence[v4->index()]));

		if (deviation_pre > deviation_post)
			mesh->flipEdgeTriangle(*e_it);
	}
}

void SurfaceRemesh::iso_tangential_relaxation(PolyMesh* mesh)
{
	mesh->updateMeshNormal();
	for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
	{
		if (mesh->isBoundary(*v_it)) continue;
		double count = 0.0;
		Eigen::Vector3d p = (*v_it)->position();
		Eigen::Vector3d q(0.0, 0.0, 0.0);
		for (VertexVertexIter vv_it = mesh->vv_iter(*v_it); vv_it.isValid(); ++vv_it)
		{
			q += (*vv_it)->position();
			++count;
		}
		q /= count;
		Eigen::Vector3d n = (*v_it)->normal();
		n.normalize();
		(*v_it)->setPosition(q + (n.dot(p - q)) * n);
	}
}

double SurfaceRemesh::iso_calculate_target_lenth_shortest(PolyMesh* mesh)
{
	double edge_min = std::numeric_limits<double>::max();
	double edge_max = std::numeric_limits<double>::min();

	for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		edge_min = std::min(edge_min, (*e_it)->length());
		edge_max = std::max(edge_max, (*e_it)->length());
	}
	double target_edge_length = edge_min * 1.5;
	target_edge_length = std::min(target_edge_length, edge_max);
	return target_edge_length;
}


double SurfaceRemesh::iso_calculateTargetEdgeLength(PolyMesh* mesh)
{

	double target_edge_length = 0.0;
  //std::cout << mesh->numEdges() << std::endl;
	for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		target_edge_length += (*e_it)->length();
	}
	target_edge_length /= mesh->numEdges();
	return target_edge_length;
}





void SurfaceRemesh::iso_remesh(Mesh& mesh) {
	//std::cout << "use iso_remesh\n";
	PolyMesh half_mesh;
	for (int i = 0; i < mesh.Vertex.rows(); ++i) {
		half_mesh.addVertex(mesh.Vertex.row(i), Eigen::Vector2d(0, 0));
	}
	for (int i = 0; i < mesh.Topo.rows(); ++i) {
		auto t = std::vector<size_t>{ size_t(mesh.Topo(i, 0)), size_t(mesh.Topo(i, 1)), size_t(mesh.Topo(i, 2)) };
		half_mesh.addPolyFace(t);
	}
	double target_edge_length, high, low;
	target_edge_length = iso_calculate_target_lenth_shortest(&half_mesh);
	high = 4.0 / 3.0 * target_edge_length;
	low = 4.0 / 5.0 * target_edge_length;
	for (int i = 0; i < 10; i++)
	{
		SurfaceRemesh::iso_split_long_edges(&half_mesh, high);
		SurfaceRemesh::iso_collapse_short_edges(&half_mesh, high, low);
		SurfaceRemesh::iso_equalize_valences(&half_mesh);
		SurfaceRemesh::iso_tangential_relaxation(&half_mesh);
	}

	mesh.Vertex.resize(half_mesh.numVertices(), 3);
	int v_id = 0, f_id = 0;
	for (VertexIter iter = half_mesh.vertices_begin(); iter != half_mesh.vertices_end(); iter++) {
		mesh.Vertex.row(v_id++) = (*iter)->position();
	}
	mesh.Topo.resize(half_mesh.numPolygons(), 3);
	mesh.Masks.resize(half_mesh.numPolygons(), 1);
	for (FaceIter iter = half_mesh.polyfaces_begin(); iter != half_mesh.polyfaces_end(); iter++) {
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