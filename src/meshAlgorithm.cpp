#include <assert.h>
#include <igl/bfs_orient.h>
#include <igl/remove_unreferenced.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
#include <igl/PI.h>
#include <igl/flipped_triangles.h>
#include <igl/topological_hole_fill.h>
#include <igl/upsample.h>
#include <igl/list_to_matrix.h>
#include <igl/writeSTL.h>

#include <set>
#include <unordered_map>
#include <cmath>
#include <algorithm> // std::move_backward
#include <random> // std::default_random_engine
#include <chrono> // std::chrono::system_clock
#include <iostream>
#include <set>
#include <array>
#include "meshAlgorithm.h"
#undef min
#undef max
using namespace std;

/**
 * the start point is (start_x, start_y, start_z), the orient is (end_x, end_y, end_z), angle is PI * angle. angle \in (0, 2).
 * @param rotateVec is the param of rotate.
 * @param V
 * @param T
 * @return 1/0
 */
bool MESHIO::rotatePoint(vector<double> rotateVec, Mesh& mesh)
{
	if (rotateVec.size() == 0) return 0;
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	double start_x = 0.0, start_y = 0.0, start_z = 0.0;
	double end_x = 0.0, end_y = 0.0, end_z = 0.0;
	double angle = 0.0;
	if (rotateVec.size() > 0) {
		if (rotateVec.size() == 4) {
			end_x = rotateVec[0];
			end_y = rotateVec[1];
			end_z = rotateVec[2];
			angle = rotateVec[3];
		}
		else if (rotateVec.size() == 7) {
			start_x = rotateVec[0]; start_y = rotateVec[1]; start_z = rotateVec[2];
			end_x = rotateVec[3]; end_y = rotateVec[4]; end_z = rotateVec[5];
			angle = rotateVec[6];
		}
		else {
			std::cout << "The format is Error.Format is (start_x, start_y, start_z, end_x, end_y, end_z, angle) or (end_x, end_y, end_z, angle). angle value scale is (0, 2).";
			return -1;
		}
	}

	if (rotateVec.size() == 7) {
		std::cout << "Rotating\n";
		for (int i = 0; i < V.rows(); i++) {
			V(i, 0) -= start_x;
			V(i, 1) -= start_y;
			V(i, 2) -= start_z;
		}
	}
	if (rotateVec.size() > 0) {
		Eigen::AngleAxisd rotationVector(M_PI * angle, Eigen::Vector3d(end_x, end_y, end_z));
		Eigen::Matrix3d rotationMatrix = Eigen::Matrix3d::Identity();
		rotationMatrix = rotationVector.toRotationMatrix();
		Eigen::MatrixXd tmp = rotationMatrix * V.transpose();
		V = tmp.transpose();
	}
	if (rotateVec.size() == 7) {
		for (int i = 0; i < V.rows(); i++) {
			V(i, 0) += start_x;
			V(i, 1) += start_y;
			V(i, 2) += start_z;
		}
		std::cout << "Rotated\n";
	}
	return 1;
}

/**
 * @brief Create a Box object
 *
 * @param create_box
 * @param mesh
 * @return true
 * @return false
 */
bool MESHIO::createBox(std::vector<double> create_box, Mesh& mesh)
{
	/***************************
	 *          | Z
	 * 			|
	 * 			|
	 * 			|__ __ __ __ Y
	 * 		   /o
	 * 		  /
	 *       / X
	/***************************
	 *    4 _________2
	 *     /|        /|
	 *   3/_|______1/ |
	 *    | |8_____|__|6
	 *    | /      | /
	 *  7 |/_______|/5
	 ****************************/

	vector<int> boxTri = {
			1, 4, 2, // up
			1, 3, 4, // up
			5, 6, 8, // down
			5, 8, 7, // down
			5, 3, 1, // front
			5, 7, 3, // front
			6, 2, 4, // behind
			6, 4, 8, // behind
			7, 8, 4, // left
			7, 4, 3, // left
			5, 1, 2, // right
			5, 2, 6  // right
	};

	int box_num = create_box.size() / 6;

	mesh.Vertex.resize(box_num * 8, 3);
	mesh.Topo.resize(boxTri.size() / 3 * box_num, 3);
	mesh.Masks.resize(mesh.Topo.rows(), 1);

	std::cout << "Create the number of box is " << box_num << std::endl;
	std::cout << "Create the number of point is " << mesh.Vertex.rows() << std::endl;
	std::cout << "Create the number of tritopo is " << mesh.Topo.rows() << std::endl;

	int base_v = 0;
	int base_t = 0;
	int base_c = 0;

	for (int i = 0; i < create_box.size(); i += 6)
	{
		Eigen::Vector3d minn, maxx;
		minn = Eigen::Vector3d(create_box[i + 0], create_box[i + 1], create_box[i + 2]).cwiseMin(Eigen::Vector3d(create_box[i + 3], create_box[i + 4], create_box[i + 5]));
		maxx = Eigen::Vector3d(create_box[i + 0], create_box[i + 1], create_box[i + 2]).cwiseMax(Eigen::Vector3d(create_box[i + 3], create_box[i + 4], create_box[i + 5]));

		Eigen::Vector3d V_lst[8];

		mesh.Vertex.row(base_v + 0) = maxx;
		mesh.Vertex.row(base_v + 1) = Eigen::Vector3d(minn.x(), maxx.y(), maxx.z());
		mesh.Vertex.row(base_v + 2) = Eigen::Vector3d(maxx.x(), minn.y(), maxx.z());
		mesh.Vertex.row(base_v + 3) = Eigen::Vector3d(minn.x(), minn.y(), maxx.z());
		mesh.Vertex.row(base_v + 4) = Eigen::Vector3d(maxx.x(), maxx.y(), minn.z());
		mesh.Vertex.row(base_v + 5) = Eigen::Vector3d(minn.x(), maxx.y(), minn.z());
		mesh.Vertex.row(base_v + 6) = Eigen::Vector3d(maxx.x(), minn.y(), minn.z());
		mesh.Vertex.row(base_v + 7) = minn;

		for (int j = 0; j < boxTri.size() / 3; j++)
		{
			mesh.Masks(base_t + j, 0) = base_c;
			mesh.Topo.row(base_t + j) = Eigen::Vector3i(base_v + boxTri[j * 3 + 0] - 1,
				base_v + boxTri[j * 3 + 1] - 1,
				base_v + boxTri[j * 3 + 2] - 1);
		}
		base_v += 8;
		base_t += 12;
		base_c += 1;
	}

	return true;
}

Eigen::Vector3d MESHIO::ave_normal(Mesh &mesh)
{
	Eigen::Vector3d res(0, 0, 0);
	double all_area = 0.0;

	for(int i = 0; i < mesh.Topo.rows(); ++i){
		int a = mesh.Topo(i, 0);
		int b = mesh.Topo(i, 1);
		int c = mesh.Topo(i, 2);
		Eigen::Vector3d ba = mesh.Vertex.row(b) - mesh.Vertex.row(a);
		Eigen::Vector3d ca = mesh.Vertex.row(c) - mesh.Vertex.row(a);
		Eigen::Vector3d normal = ba.cross(ca);
		double area = normal.norm() * 0.5;
		normal = normal.normalized();

		res += normal * area;
		all_area += area;
	}

	return res / all_area;
}

bool MESHIO::Normalize(Mesh& mesh) {
	Eigen::Vector3d maxC(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
	Eigen::Vector3d minC = -maxC;
	for (int j = 0; j < mesh.Vertex.rows(); j++)
		for (int i = 0; i < 3; i++) {
			maxC(i) = std::max(mesh.Vertex(j, i), maxC(i));
			minC(i) = std::min(mesh.Vertex(j, i), minC(i));
		}
	double scale = 0;
	for (int i = 0; i < 3; i++) {
		scale = std::max(scale, maxC(i) - minC(i));
	}
	for (int j = 0; j < mesh.Vertex.rows(); j++)
		for (int i = 0; i < 3; i++) {
			mesh.Vertex(j, i) = (mesh.Vertex(j, i) - minC(i)) / scale;
		}
	return true;
}
bool MESHIO::removeBox(Mesh& mesh)
{
	Eigen::VectorXi component;
	Eigen::MatrixXi tmp = mesh.Topo;
	igl::bfs_orient(tmp, mesh.Topo, component);

	map<int, double> vols;
	for (int j = 0; j < component.rows(); j++) {
		if (vols.find(component(j)) == vols.end())
			vols[component(j)] = 0;
		vols[component(j)] += std::fabs(igl::volume_single(
			static_cast<Eigen::RowVector3d>(mesh.Vertex.row(mesh.Topo(j, 0))),
			static_cast<Eigen::RowVector3d>(mesh.Vertex.row(mesh.Topo(j, 1))),
			static_cast<Eigen::RowVector3d>(mesh.Vertex.row(mesh.Topo(j, 2))),
			Eigen::RowVector3d(0, 0, 0)));
	}

	auto maxi = std::max_element(vols.begin(), vols.end(),
		[](std::pair<int, double> t1, std::pair<int, double> t2) {return t1.second < t2.second; })->first;
	std::map<int, int> index_map; // new to old
	std::map<int, int> index_rmap; // old to new
	int topo_count = 0;
	for (int j = 0; j < component.rows(); j++) {
		if (component(j) != maxi) {
			topo_count++;
			for (int k = 0; k < 3; k++) {
				int index = mesh.Topo(j, k);
				if (index_rmap.find(index) == index_rmap.end()) {
					int newid = index_rmap.size();
					index_rmap[index] = newid;
					index_map[newid] = index;
				}
			}
		}
	}
	Eigen::MatrixXd new_vertex;
	Eigen::MatrixXi new_topo;
	new_topo.resize(topo_count, 3);
	new_vertex.resize(index_map.size(), 3);
	for (int i = 0; i < index_map.size(); i++) {
		new_vertex.row(i) = static_cast<Eigen::RowVector3d>(mesh.Vertex.row(index_map[i]));
	}
	topo_count = 0;
	for (int j = 0; j < component.rows(); j++) {
		if (component(j) != maxi) {
			for (int k = 0; k < 3; k++) {
				new_topo(topo_count, k) = index_rmap[mesh.Topo(j, k)];
			}
			topo_count++;
		}
	}

	if (mesh.Masks.rows() == mesh.Vertex.rows()) {
		for (int i = 0; i < index_map.size(); i++) {
			mesh.Masks.row(i) = mesh.Masks.row(index_map[i]);
		}
		mesh.Masks.conservativeResize(mesh.Vertex.rows(), mesh.Masks.cols());
	}
	if (mesh.Masks.rows() == mesh.Topo.rows()) {
		topo_count = 0;
		for (int j = 0; j < component.rows(); j++) {
			if (component(j) != maxi) {
				mesh.Masks.row(topo_count) = mesh.Masks.row(j);
				topo_count++;
			}
		}
		mesh.Masks.conservativeResize(mesh.Topo.rows(), mesh.Masks.cols());
	}

	mesh.Vertex = new_vertex;
	mesh.Topo = new_topo;

	return true;
}
bool MESHIO::scale(Mesh& mesh, double factor)
{
	mesh.Vertex = mesh.Vertex * factor;
	return true;
}

/**
 * Add bounding box.
 * @param boxVec
 * @param V
 * @param T
 * @return
 */
bool MESHIO::addBox(vector<double> boxVec, Mesh& mesh)
{
	if (boxVec.size() == 0) return 0;
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	if (boxVec.size() != 3) {
		std::cout << "the number of boxVec's param is not equal to 3." << std::endl;
		return 0;
	}

	std::cout << "Boxing\n";

	double len = boxVec[0], width = boxVec[1], hight = boxVec[2];

	vector<double> tmpx, tmpy, tmpz;

	for (int i = 0; i < V.rows(); i++)
	{
		tmpx.push_back(V(i, 0));
		tmpy.push_back(V(i, 1));
		tmpz.push_back(V(i, 2));
	}

	std::sort(tmpx.begin(), tmpx.end());
	std::sort(tmpy.begin(), tmpy.end());
	std::sort(tmpz.begin(), tmpz.end());

	int tmpid = V.rows() / 2;

	double midx = tmpx[tmpid];
	double midy = tmpy[tmpid];
	double midz = tmpz[tmpid];

	double maxlength = std::max((tmpx.back() - tmpx.front()) / 2, (tmpy.back() - tmpy.front()) / 2);
	maxlength = std::max(maxlength, (tmpz.back() - tmpz.front()) / 2);
	double lengthx = maxlength * boxVec[0];
	double lengthy = maxlength * boxVec[1];
	double lengthz = maxlength * boxVec[2];

	std::cout << midx << " " << midy << " " << midz << '\n';
	std::cout << tmpx[tmpx.size() - 1] - tmpx[0] << " " << tmpy[tmpy.size() - 1] - tmpy[0] << " " << \
		tmpz[tmpz.size() - 1] - tmpz[0] << '\n';

	tmpx.clear(); tmpy.clear(); tmpz.clear();
	Mesh box_mesh;

	vector<int> boxTri = {
			1, 4, 2, // up
			1, 3, 4, // up
			5, 6, 8, // down
			5, 8, 7, // down
			5, 3, 1, // front
			5, 7, 3, // front
			6, 2, 4, // behind
			6, 4, 8, // behind
			7, 8, 4, // left
			7, 4, 3, // left
			5, 1, 2, // right
			5, 2, 6  // right
	};

	box_mesh.Topo.resize(12, 3);
	for (int i = 0; i < 12; i++)
		for (int k = 0; k < 3; k++)
			box_mesh.Topo(i, k) = boxTri[3 * (i)+k] - 1;
	box_mesh.Vertex.resize(8, 3);
	int V_index = 0;
	for (double i = lengthx; i > -2 * lengthx; i -= 2 * lengthx) // up and down
	{
		for (double j = lengthy; j > -2 * lengthy; j -= 2 * lengthy) // left and right
		{
			for (double k = lengthz; k > -2 * lengthz; k -= 2 * lengthz) // forward and back
			{
				box_mesh.Vertex(V_index, 0) = midx + i;
				box_mesh.Vertex(V_index, 1) = midy + j;
				box_mesh.Vertex(V_index, 2) = midz + k;
				V_index++;
			}
		}
	}
	Mesh refine_box_mesh;

	int refine_time = 4;
	if (boxVec.size() >= 3) {
		refine_time = int(boxVec[3]);
		if (refine_time <= 0 || refine_time >= 10)
			refine_time = 4;
	}

	igl::upsample(box_mesh.Vertex, box_mesh.Topo, refine_box_mesh.Vertex, refine_box_mesh.Topo, refine_time);
	//igl::writeSTL("test.stl", refine_box_mesh.Vertex, refine_box_mesh.Topo);

	int num_vertex = mesh.Vertex.rows();
	Eigen::MatrixXd init_V = mesh.Vertex;
	Eigen::MatrixXi init_F = mesh.Topo;
	igl::cat(1, init_V, refine_box_mesh.Vertex, mesh.Vertex);
	Eigen::MatrixXi offset_topo = refine_box_mesh.Topo.array() + num_vertex;
	igl::cat(1, init_F, offset_topo, mesh.Topo);

	std::cout << "Boxed\n";

	return 1;
}

bool MESHIO::reverseOrient(Eigen::MatrixXi& T) {
	std::cout << "Reversing\n";
	for (int i = 0; i < T.rows(); i++)
	{
		int t = T(i, 0);
		T(i, 0) = T(i, 2);
		T(i, 2) = t;
	}
	std::cout << "reversed\n";
	return 1;
}

bool MESHIO::repair(Mesh& mesh)
{
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	std::cout << "Vertex number is  " << V.rows() << " X " << V.cols() << "  before clean. \n";
	std::cout << "Cell number is  " << T.rows() << " X " << T.cols() << "  before clean. \n";
	std::cout << "Attribute number is  " << M.rows() << " X " << M.cols() << "  before clean. \n";
	// 检测面积为0的单元个数

	Eigen::MatrixXd Tri;
	Tri.resize(3, 3);
	vector<int> emptyTri;
	for (int i = 0; i < T.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Tri(j, 0) = V(T(i, j), 0);
			Tri(j, 1) = V(T(i, j), 1);
			Tri(j, 2) = V(T(i, j), 2);
		}
		Eigen::Vector3d lin1(Tri(1, 0) - Tri(0, 0), Tri(1, 1) - Tri(0, 1), Tri(1, 2) - Tri(0, 2));
		Eigen::Vector3d lin2(Tri(2, 0) - Tri(0, 0), Tri(2, 1) - Tri(0, 1), Tri(2, 2) - Tri(0, 2));

		Eigen::Vector3d crossResult = lin1.cross(lin2);
		double area = crossResult.norm();
		if (area < 1e-8)
		{
			emptyTri.push_back(i);
		}
	}

	std::cout << "There are " << emptyTri.size() << " cells whose area is equal to zero.\n";

	// 去除面积为0的单元

	struct node
	{
		int oldid;
		Eigen::Vector3d point;
	};
	vector<node> vec(V.rows(), node());
	map<int, int> mpid;

	for (int i = 0; i < V.rows(); i++)
	{
		vec[i].oldid = i;
		vec[i].point = V.row(i);
	}

	std::sort(vec.begin(), vec.end(), [](const node& a, const node& b)
		{
			if (a.point.x() != b.point.x()) return a.point.x() < b.point.x();
			if (a.point.y() != b.point.y()) return a.point.y() < b.point.y();
			if (a.point.z() != b.point.z()) return a.point.z() < b.point.z();
			return a.oldid < b.oldid;
		});

	int curid = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		while (i + 1 < vec.size() && (vec[i].point - vec[i + 1].point).norm() < 1e-8)
		{
			mpid[vec[i].oldid] = curid;
			i++;
		}
		mpid[vec[i].oldid] = curid;
		curid++;
	}

	V.resize(curid, V.cols());

	for (int i = 0; i < vec.size(); i++)
	{
		V.row(mpid[vec[i].oldid]) = vec[i].point;
		while (i + 1 < vec.size() && (vec[i].point - vec[i + 1].point).norm() < 1e-8) i++;
	}

	Eigen::MatrixXi T2;
	Eigen::MatrixXi M2;
	T2.resize(T.rows() - emptyTri.size(), T.cols());
	M2.resize(T.rows() - emptyTri.size(), M.cols());

	int locT2 = 0, locEmptyTri = 0;
	for (int i = 0; i < T.rows(); i++)
	{
		if (emptyTri[locEmptyTri] == i)
		{
			locEmptyTri++;
			continue;
		}
		for (int j = 0; j < T.cols(); j++)
		{
			T2(locT2, j) = mpid[T(i, j)];
		}
		locT2++;
	}

	T = T2;
	// M = M2;

	std::cout << "Vertex number is  " << V.rows() << " X " << V.cols() << "  after clean. \n";
	std::cout << "Cell number is  " << T.rows() << " X " << T.cols() << "  after clean. \n";
	std::cout << "Attribute number is  " << M.rows() << " X " << M.cols() << "  after clean. \n";
	std::cout << "Clean all cell whose area is equal to zero. " << '\n';
	return 0;
}

void MESHIO::checkOrientation(Mesh& mesh) {
	vector<vector<double>> point_list(mesh.Vertex.rows(), vector<double>(mesh.Vertex.cols(), 0));
	vector<vector<int>> facet_list(mesh.Topo.rows(), vector<int>(mesh.Topo.cols(), 0));
	for (int i = 0; i < mesh.Vertex.rows(); i++) {
		for (int j = 0; j < mesh.Vertex.cols(); j++) {
			point_list[i][j] = mesh.Vertex(i, j);
		}
	}
	for (int i = 0; i < mesh.Topo.rows(); i++) {
		for (int j = 0; j < mesh.Topo.cols(); j++) {
			facet_list[i][j] = mesh.Topo(i, j);
		}
	}

	//[](vector<int>& a,vector<int>& b) {}
	vector<vector<int>>  facet_list_copy = facet_list;
	vector<int> block_mark;
//	TIGER::resetOrientation(point_list, facet_list, block_mark);

	std::set<vector<int>> qs;
	for (auto i : facet_list) {
		int id = std::min_element(i.begin(), i.end()) - i.begin();
		qs.insert(vector<int>{i[id], i[(id + 1) % 3], i[(id + 2) % 3]});
	}

	int message_count = 0;
	for (int j = 0; j < facet_list_copy.size(); j++) {
		swap(facet_list_copy[j][0], facet_list_copy[j][1]);
		if (message_count > 100) {
			cout << "Too much oritation error! program finished.";
			break;
		}
		int id = std::min_element(facet_list_copy[j].begin(), facet_list_copy[j].end()) - facet_list_copy[j].begin();
		vector<int> vep = vector<int>{ facet_list_copy[j][id], facet_list_copy[j][(id + 1) % 3], facet_list_copy[j][(id + 2) % 3] };

		if (qs.find(vep) == qs.end()) {
			message_count++;
			cout << "ortient error found in face id=" << j << endl;
			for (auto k : vep) {
				cout << k << " ";
			}
			cout << endl;
			//	break;
		}
	}

	return;
}
bool MESHIO::resetOrientation(Mesh& mesh, bool reset_mask,std::vector<int> saveid) {
	vector<vector<double>> point_list(mesh.Vertex.rows(), vector<double>(mesh.Vertex.cols(), 0));
	vector<vector<int>> facet_list(mesh.Topo.rows(), vector<int>(mesh.Topo.cols(), 0));
	for (int i = 0; i < mesh.Vertex.rows(); i++) {
		for (int j = 0; j < mesh.Vertex.cols(); j++) {
			point_list[i][j] = mesh.Vertex(i, j);
		}
	}
	for (int i = 0; i < mesh.Topo.rows(); i++) {
		for (int j = 0; j < mesh.Topo.cols(); j++) {
			facet_list[i][j] = mesh.Topo(i, j);
		}
	}
	vector<int> block_mark;
//	TIGER::resetOrientation(point_list, facet_list, block_mark);
	for (int i = 0; i < mesh.Vertex.rows(); i++) {
		for (int j = 0; j < mesh.Vertex.cols(); j++) {
			mesh.Vertex(i, j) = point_list[i][j];
		}
	}
	for (int i = 0; i < mesh.Topo.rows(); i++) {
		for (int j = 0; j < mesh.Topo.cols(); j++) {
			mesh.Topo(i, j) = facet_list[i][j];
		}
	}
	if (reset_mask) {
		if (mesh.Masks.rows() != mesh.Topo.rows()) {
			mesh.Masks.resize(mesh.Topo.rows(), mesh.Masks.cols());
			for (int i = 0; i < mesh.Topo.rows(); i++) {
				mesh.Masks(i, 0) = block_mark[i];
			}
		}
		else {
			std::map<int,int> id_map;
			int count = 0;
			for (auto i : saveid)
				if(id_map.find(i)==id_map.end())
					id_map[i]= count++;
			for (int i = 0; i < mesh.Topo.rows(); i++) {
				if (id_map.find(mesh.Masks(i, 0)) != id_map.end()) {
					mesh.Masks(i, 0) = id_map[mesh.Masks(i, 0)];// block_mark[i];
				}
				else {
					mesh.Masks(i, 0)= block_mark[i]+ id_map.size();
				}
			
			}
		}
	}
	cout << "Orientation reset" << endl;
	return 1;
}
/**
 * sort eigen vector
 * @param in_vec
 * @param out_vec
 * @param ind   index
 */
void eigen_sort_3d(const Eigen::MatrixXd& in_vec, Eigen::MatrixXd& out_vec, Eigen::VectorXi& ind, double eps)
{
	ind = Eigen::VectorXi::LinSpaced(in_vec.rows(), 0, in_vec.rows() - 1);
	auto rule = [&in_vec, &eps](int i, int j)->bool {
		if (fabs(in_vec(i, 0) - in_vec(j, 0)) > eps) return in_vec(i, 0) < in_vec(j, 0);
		if (fabs(in_vec(i, 1) - in_vec(j, 1)) > eps) return in_vec(i, 1) < in_vec(j, 1);
		if (fabs(in_vec(i, 2) - in_vec(j, 2)) > eps) return in_vec(i, 2) < in_vec(j, 2);
		return i < j;
		};
	std::sort(ind.data(), ind.data() + ind.size(), rule);
	out_vec.resize(in_vec.rows(), in_vec.cols());
	for (int i = 0; i < in_vec.rows(); i++)
	{
		out_vec.row(i) << in_vec.row(ind(i));
	}
}

/**
 * Remove dulplicate point.
 * @param V
 * @param T
 */
bool MESHIO::removeDulplicatePoint(Eigen::MatrixXd& V, Eigen::MatrixXi& F, double eps) {
	Eigen::MatrixXd V_in = V;
	Eigen::MatrixXi F_in = F;
	Eigen::MatrixXd V_sort;
	Eigen::VectorXi ind;
	eigen_sort_3d(V_in, V_sort, ind, eps);
	int cnt = 0;
	std::unordered_map<int, int> mpid;
	std::vector<int> new_point_lst;
	for (int i = 0; i < V_sort.rows(); i++)
	{
		while (i + 1 < V_sort.rows())
		{
			if ((V_sort.row(i + 1) - V_sort.row(i)).norm() < eps)
			{
				mpid[ind[i]] = cnt;
				i++;
			}
			else break;
		}

		mpid[ind[i]] = cnt;
		new_point_lst.push_back(ind[i]);
		cnt++;
	}

	V.resize(new_point_lst.size(), 3);
	F.resize(F_in.rows(), 3);

	for (int i = 0; i < new_point_lst.size(); i++)
	{
		V.row(i) = V_in.row(new_point_lst[i]);
	}

	for (int i = 0; i < F_in.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			F(i, j) = mpid[F_in(i, j)];
		}
	}
	Eigen::VectorXi I;
	auto V_old = V;
	auto F_old = F;
	igl::remove_unreferenced(V_old, F_old, V, F, I);

	std::cout << "Remove " << V_in.rows() - V.rows() << " duplicate points\n";
	std::cout << "Remove " << F_in.rows() - F.rows() << " facets\n";
	return true;
}
void MESHIO::get_boundary_loop(std::vector<std::array<int, 2>> &edge,
                       std::vector<std::vector<int>> &loop) {
    std::vector<int> bv;  // boundary vertex.
    for (int i = 0; i < edge.size(); ++i) {
        bv.push_back(edge[i][0]);
        bv.push_back(edge[i][1]);
    }
    std::sort(bv.begin(), bv.end());
    bv.erase(std::unique(bv.begin(), bv.end()), bv.end());
    std::map<int, int> bv_to_nbv;  // to new boundary vertex.
    for (int i = 0; i < bv.size(); ++i) {
        bv_to_nbv[bv[i]] = i;
    }

    std::vector<std::vector<int>> G(bv.size(), std::vector<int>());
    std::vector<bool> vis(bv.size(), 0);
    std::vector<int> du(bv.size(), 0);
    for (int i = 0; i < edge.size(); ++i) {
        G[bv_to_nbv[edge[i][0]]].push_back(bv_to_nbv[edge[i][1]]);
        G[bv_to_nbv[edge[i][1]]].push_back(bv_to_nbv[edge[i][0]]);
        du[bv_to_nbv[edge[i][0]]]++;
        du[bv_to_nbv[edge[i][1]]]++;
    }

    for (int i = 0; i < bv.size(); ++i) {
        if (!vis[i]) {
            std::vector<int> path;
            std::vector<std::vector<int>> loop_lst;
            path.push_back(i);
            boundary_loop_dfs(i, -1, vis, du, G, path, loop_lst);
            for (auto &t : loop_lst) loop.push_back(t);
        }
    }

    // update vertex id
    for (int i = 0; i < loop.size(); ++i) {
        for (int j = 0; j < loop[i].size(); ++j) {
            loop[i][j] = bv[loop[i][j]];
        }
    }
}
void MESHIO::boundary_loop_dfs(int cur,
                       int pre,
                       std::vector<bool> &vis,
                       std::vector<int> &du,
                       std::vector<std::vector<int>> &G,
                       std::vector<int> &path,
                       std::vector<std::vector<int>> &loop_lst) {
    if (vis[cur]) {
        std::vector<int> tmp;
        for (int i = path.size() - 2; i >= 0; i--) {
            if (path[i] != cur) {
                tmp.push_back(path[i]);
            } else {
                tmp.push_back(path[i]);
                break;
            }
        }
        if (tmp.size() > 1) {
            loop_lst.push_back(tmp);
            for (auto &t : tmp) du[t] -= 2;
        }
        return;
    }

    vis[cur] = true;
    for (int i = 0; i < G[cur].size(); i++) {
        if (G[cur][i] == pre || du[G[cur][i]] == 0) continue;
        path.push_back(G[cur][i]);
        boundary_loop_dfs(G[cur][i], cur, vis, du, G, path, loop_lst);
        path.pop_back();
    }
};

// void MESHIO::dfs_get_loop2(int cur, int pre, std::vector<bool>& vis, std::vector<std::vector<int>>& G, std::vector<int>& path, std::vector<std::vector<int>>& loop_lst)
// {
// 	if (vis[cur])
// 	{
// 		std::vector<int> tmp;
// 		for (int i = path.size() - 2; i >= 0; i--)
// 		{
// 			if (path[i] != cur)
// 			{
// 				tmp.push_back(path[i]);
// 			}
// 			else
// 			{
// 				tmp.push_back(path[i]);
// 				break;
// 			}
// 		}
// 		loop_lst.push_back(tmp);
// 		return;
// 	}

// 	vis[cur] = 1;
// 	for (int i = 0; i < G[cur].size(); i++)
// 	{
// 		if (G[cur][i] == pre) continue;
// 		path.push_back(G[cur][i]);
// 		dfs_get_loop2(G[cur][i], cur, vis, G, path, loop_lst);
// 		path.pop_back();
// 	}
// 	vis[cur] = 0;
// }

void MESHIO::boundary_loop_by_dfs2(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& bnd)
{
	double eps = 1e-6;
	Eigen::MatrixXd V_tmp;
	Eigen::MatrixXi F_tmp;
	Eigen::VectorXi ind;

	V_tmp = V;
	F_tmp = F;

	bool b_dulplicated_point = false;
	if (b_dulplicated_point)
	{
		for (int i = 0; i < V.rows(); i++)
		{
			for (int j = i + 1; j < V.rows(); j++)
				if ((V.row(i) - V.row(j)).norm() < eps)
				{
					std::cout << "Warning : boundary loop have dulplicate point.\n";
					break;
				}
		}
	}

	// To find all boundary for building graph.
	std::vector<std::vector<int>> G(V.rows(), std::vector<int>());
	std::map<std::array<int, 2>, int> mp;
	std::vector<bool> vis(V.rows(), 0);
	std::vector<int> du(V.rows(), 0);
	std::vector<std::vector<int>> loop_lst;
	std::vector<int> path;
	std::set<int> num_p_set;

	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int a = F(i, j);
			int b = F(i, (j + 1) % 3);
			std::array<int, 2> tmp = { std::min(a, b), std::max(a, b) };
			mp[tmp]++;
		}
	}

	for (auto iter = mp.begin(); iter != mp.end(); iter++)
	{
		if (iter->second == 1)
		{
			G[iter->first[0]].push_back(iter->first[1]);
			G[iter->first[1]].push_back(iter->first[0]);
			// std::cout << iter->first[0] << " " << iter->first[1] << '\n';
			du[iter->first[0]] ++;
			du[iter->first[1]] ++;
			num_p_set.insert(iter->first[0]);
			num_p_set.insert(iter->first[1]);
		}
	}

	int start = *num_p_set.begin();
	path.push_back(start);
	// MESHIO::dfs_get_loop2(start, -1, vis, G, path, loop_lst);
	MESHIO::boundary_loop_dfs(start, -1, vis, du, G, path, loop_lst);

	for (int i = 0; i < loop_lst.size(); i++)
	{
		std::cout << "#BK The " << i << "th loop number is : " << loop_lst[i].size() << '\n';
	}

	sort(loop_lst.begin(), loop_lst.end(), [&V](const std::vector<int>& a, const std::vector<int>& b) {
		double len_a = 0.0;
		double len_b = 0.0;

		for (int i = 0; i < a.size(); i++)
		{
			len_a += (V.row(i) - V.row((i + 1) % a.size())).norm();
		}

		for (int i = 0; i < b.size(); i++)
		{
			len_b += (V.row(i) - V.row((i + 1) % b.size())).norm();
		}

		return len_a > len_b;
		});

	int s = loop_lst[0][0];
	int e = loop_lst[0][1];
	int ok = 0;
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int a = F(i, j);
			int b = F(i, (j + 1) % 3);
			if (a == s && b == e)
			{
				ok = 1;
			}
		}
		if (ok) break;
	}

	if (!ok)
	{
		std::reverse(loop_lst[0].begin(), loop_lst[0].end());
	}

	bnd.resize(loop_lst[0].size(), 1);
	for (int i = 0; i < loop_lst[0].size(); i++)
	{
		bnd(i) = loop_lst[0][i];
	}
}

bool MESHIO::buildTuttleParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv)
{
	Eigen::MatrixXi C;
	Eigen::VectorXi bnd;

	igl::bfs_orient(T_3d, T_3d, C);

	// remove duplicate point
	double minn_size = std::numeric_limits<double>::max();

	for (int i = 0; i < T_3d.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int a = T_3d(i, j);
			int b = T_3d(i, (j + 1) % 3);
			minn_size = std::min(minn_size, (V_3d.row(a) - V_3d.row(b)).norm());
		}
	}
	minn_size = sqrt(minn_size);
	Eigen::MatrixXd tmpV = V_3d;
	Eigen::MatrixXi tmpT = T_3d;
	Eigen::VectorXi SVI, SVJ;
	igl::remove_duplicate_vertices(tmpV, tmpT, minn_size / 1000.0, V_3d, SVI, SVJ, T_3d);

	// removeDulplicatePoint(V_3d, T_3d, minn_size / 1000.0);
	int V_N = V_3d.rows();

	boundary_loop_by_dfs2(V_3d, T_3d, bnd); // mine

	// Map the boundary to a circle, preserving edge proportions
	Eigen::MatrixXd bnd_uv;
	igl::map_vertices_to_circle(V_3d, bnd, bnd_uv);

	// build adj relation
	std::vector<std::vector<int>> adj(V_N);
	for (int i = 0; i < T_3d.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int a = T_3d(i, j);
			int b = T_3d(i, (j + 1) % 3);
			adj[a].push_back(b);
			adj[b].push_back(a);
		}
	}

	typedef Eigen::Triplet<double> T;
	typedef Eigen::SparseMatrix<double> SMatrix;

	std::vector<T> tripletlist;
	Eigen::VectorXd bu = Eigen::VectorXd::Zero(V_N);
	Eigen::VectorXd bv = Eigen::VectorXd::Zero(V_N);

	std::unordered_map<int, bool> is_boundary;
	std::unordered_map<int, Eigen::Vector2d> mp_bnd;
	for (int i = 0; i < bnd.rows(); i++)
	{
		is_boundary[bnd(i, 0)] = 1;
		mp_bnd[bnd(i, 0)] = bnd_uv.row(i);
	}

	for (int i = 0; i < V_N; ++i)
	{
		if (is_boundary[i])
		{
			tripletlist.emplace_back(i, i, 1);
			bu(i) = mp_bnd[i].x();
			bv(i) = mp_bnd[i].y();
		}
		else {
			for (int j = 0; j < adj[i].size(); j++)
			{
				tripletlist.emplace_back(i, adj[i][j], -1);
			}
			tripletlist.emplace_back(i, i, adj[i].size());
		}
	}

	SMatrix coff(V_N, V_N);
	coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
	Eigen::SparseLU<SMatrix> solver;
	solver.compute(coff);

	Eigen::VectorXd xu = solver.solve(bu);
	Eigen::VectorXd xv = solver.solve(bv);

	V_uv.resize(V_N, 2);

	V_uv.col(0) = xu;
	V_uv.col(1) = xv;

	std::cout << "Build tuttle parameter is success.\n";
	return 1;
}

bool MESHIO::buildHarmonicParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv)
{
	Eigen::VectorXi bnd;

	// remove duplicate point
	double minn_size = std::numeric_limits<double>::max();

	for (int i = 0; i < T_3d.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int a = T_3d(i, j);
			int b = T_3d(i, (j + 1) % 3);
			minn_size = std::min(minn_size, (V_3d.row(a) - V_3d.row(b)).norm());
		}
	}
	minn_size = sqrt(minn_size);
	Eigen::MatrixXd tmpV = V_3d;
	Eigen::MatrixXi tmpT = T_3d;
	Eigen::VectorXi SVI, SVJ;
	igl::remove_duplicate_vertices(tmpV, tmpT, minn_size / 1000.0, V_3d, SVI, SVJ, T_3d);

	boundary_loop_by_dfs2(V_3d, T_3d, bnd); // mine

	// Map the boundary to a circle, preserving edge proportions
	Eigen::MatrixXd bnd_uv;
	igl::map_vertices_to_circle(V_3d, bnd, bnd_uv);

	igl::harmonic(V_3d, T_3d, bnd, bnd_uv, 1, V_uv);

	// Scale UV to make the texture more clear
	// V_uv *= 10000;

	std::cout << "Build harmonic parameter is success.\n";
	return 1;
}

bool MESHIO::shuffleSurfaceid(int num, Mesh& mesh)
{
	num = std::min(100, num);
	num = std::max(1, num);

	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	std::vector<int> random_id;

	if (M.cols() != 1) return 0;

	std::unordered_map<int, int> mp_id;
	int cnt = 0;
	for (int i = 0; i < M.rows(); i++)
	{
		if (mp_id.count(M(i, 0))) continue;
		mp_id[M(i, 0)] = cnt++;
	}

	for (int i = 0; i < cnt; ++i)
	{
		random_id.push_back(i);
	}

	for (int i = 0; i < std::min(100, num); ++i)
	{
		std::shuffle(random_id.begin(), random_id.end(), std::default_random_engine(i));
	}

	set<int> tt;
	for (int i = 0; i < M.rows(); ++i)
	{
		if (M(i, 0) == random_id[mp_id[M(i, 0)]])
		{
			tt.insert(M(i, 0));
		}
		M(i, 0) = random_id[mp_id[M(i, 0)]];
	}
	std::cout << "After shuffle the number of dulplication is : " << tt.size() << "\n";
	return 1;
}

template <
	typename DerivedF,
	typename Derivedb,
	typename VectorIndex,
	typename DerivedF_filled>
void MESHIO::topological_hole_fill(
	const Eigen::MatrixBase<DerivedF>& F,
	const Eigen::MatrixBase<Derivedb>& b,
	const std::vector<VectorIndex>& holes,
	Eigen::PlainObjectBase<DerivedF_filled>& F_filled,
	Eigen::MatrixXd& V)
{
	int n_filled_faces = 0;
	int num_holes = holes.size();
	int real_F_num = F.rows();
	const int V_rows = F.maxCoeff() + 1;

	for (int i = 0; i < num_holes; i++)
		n_filled_faces += holes[i].size();
	Eigen::MatrixXd V_ = V;
	F_filled.resize(n_filled_faces + real_F_num, 3);
	F_filled.topRows(real_F_num) = F;
	V.resize(V_rows + num_holes, 3);
	V.topRows(V_rows) = V_;

	int new_vert_id = V_rows;
	int new_face_id = real_F_num;

	for (int i = 0; i < num_holes; i++, new_vert_id++)
	{
		Eigen::Vector3d hole_v(0, 0, 0);
		for (auto& bnd_id : holes[i])
		{
			hole_v += V.row(bnd_id);
		}
		hole_v /= holes[i].size();
		V.row(new_vert_id) = hole_v;

		int cur_bnd_size = holes[i].size();
		int it = 0;
		int back = holes[i].size() - 1;
		F_filled.row(new_face_id++) << holes[i][it], holes[i][back], new_vert_id;
		while (it != back)
		{
			F_filled.row(new_face_id++)
				<< holes[i][(it + 1)],
				holes[i][(it)], new_vert_id;
			it++;
		}
	}
	assert(new_face_id == F_filled.rows());
	assert(new_vert_id == V_rows + num_holes);
}

bool MESHIO::topoFillHole(Mesh& mesh)
{
	auto& V = mesh.Vertex;
	auto& F = mesh.Topo;
	std::vector<std::vector<int>> all_bnds;
	igl::boundary_loop(F, all_bnds);
	// boundary_loop_by_dfs2(V, F, all_bnds); // mine

	// Heuristic primary boundary choice: longest
	auto primary_bnd = std::max_element(all_bnds.begin(), all_bnds.end(), [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() > b.size(); });

	Eigen::VectorXi bnd = Eigen::Map<Eigen::VectorXi>(primary_bnd->data(), primary_bnd->size());

	// if there is a hole, fill it and erase additional vertices.
	// all_bnds.erase(primary_bnd);
	Eigen::MatrixXi F_filled;
	Eigen::MatrixXd& V_filled = V;
	topological_hole_fill(F, bnd, all_bnds, F_filled, V_filled);
	F = F_filled;
	mesh.Masks.resize(0, 0);
	return true;
}


bool MESHIO::dynamicFillHole(const Eigen::MatrixXd &V,const  Eigen::MatrixXi &T, const std::vector<std::set<int>>& point_to_tris, const std::vector<int> hole, std::vector<Eigen::Vector3i>& patches){

	struct Weight
	{
		Weight(double angle = std::numeric_limits<double>::max(),
					 double area = std::numeric_limits<double>::max())
				: angle_(angle), area_(area)
		{
		}

		Weight &operator=(const Weight &rhs)
		{
			if (this == &rhs) // 自赋值检查
			{
				return *this;
			}
			angle_ = rhs.angle_;
			area_ = rhs.area_;
			return *this;
		}

		Weight operator+(const Weight &rhs) const
		{
			return Weight(std::max(angle_, rhs.angle_), area_ + rhs.area_);
		}

		bool operator<(const Weight &rhs) const
		{
			return (angle_ < rhs.angle_ ||
							(angle_ == rhs.angle_ && area_ < rhs.area_));
		}

		double angle_;
		double area_;
	};	
	
	int n = hole.size();
	std::vector<std::vector<Weight>> weight(n, std::vector<Weight>(n, Weight()));
	std::vector<std::vector<int>> res(n, std::vector<int>(n, 0));
	Weight w;

	auto TriangleArea = [](const Eigen::Vector3d &a,
												 const Eigen::Vector3d &b,
												 const Eigen::Vector3d &c)
	{
		Eigen::Vector3d u = b - a;
		Eigen::Vector3d v = c - a;
		double area = 0.5 * u.cross(v).norm();
		return area;
	};

	auto computer_angle = [](const Eigen::Vector3d& n1, const Eigen::Vector3d& n2){
		return 1.0 - n1.dot(n2);
	};

	auto computer_normal = [&](int a, int b, int c){
		Eigen::Vector3d n1 = V.row(b) - V.row(a);
		Eigen::Vector3d n2 = V.row(c) - V.row(a);
		Eigen::Vector3d n3 = n1.cross(n2);
		return n3.normalized();
	};

	auto computer_weight = [&](int i, int j, int k)
	{
		Eigen::Vector3d a = V.row(i);
		Eigen::Vector3d b = V.row(j);
		Eigen::Vector3d c = V.row(k);
		double area = TriangleArea(a, b, c);
		double angle = 0;
		Eigen::Vector3d normal = computer_normal(hole[i], hole[j], hole[k]);
		Eigen::Vector3d adj_nomal;

		if(i + 1 == j){
			std::vector<int> facets;
			std::set_intersection(point_to_tris[hole[i]].begin(), point_to_tris[hole[i]].end(), point_to_tris[hole[j]].begin(), point_to_tris[hole[j]].end(), std::back_inserter(facets));
			Eigen::Vector3i f = T.row(facets[0]);
			adj_nomal = computer_normal(f[0], f[1], f[2]);
		} else {
			adj_nomal = computer_normal(hole[i], hole[res[i][j]], hole[j]);
		}
		angle = std::max(angle, computer_angle(normal, adj_nomal));

		if(j + 1 == k){
			std::vector<int> facets;
			std::set_intersection(point_to_tris[hole[j]].begin(), point_to_tris[hole[j]].end(), point_to_tris[hole[k]].begin(), point_to_tris[hole[k]].end(), std::back_inserter(facets));
			Eigen::Vector3i f = T.row(facets[0]);
			adj_nomal = computer_normal(f[0], f[1], f[2]);
		} else {
			adj_nomal = computer_normal(hole[j], hole[res[j][k]], hole[k]);
		}
		angle = std::max(angle, computer_angle(normal, adj_nomal));

		if(i == 0 && k == n - 1){
			std::vector<int> facets;
			std::set_intersection(point_to_tris[hole[i]].begin(), point_to_tris[hole[i]].end(), point_to_tris[hole[k]].begin(), point_to_tris[hole[k]].end(), std::back_inserter(facets));
			Eigen::Vector3i f = T.row(facets[0]);
			adj_nomal = computer_normal(f[0], f[1], f[2]);
			angle = std::max(angle, computer_angle(normal, adj_nomal));
		}
		return Weight(angle, area);
	};



	for (int i = 0; i < n - 1; ++i)
	{
		weight[i][i + 1] = Weight(0, 0);
		res[i][i + 1] = -1;
	}
	for (int i = 2; i < n; ++i)
	{
		for (int j = 0; j < n - i; ++j)/// [j, j + i]
		{ 
			Weight wmin;
			int index_min = j + 1;
			for (int k = j + 1; k < j + i; ++k)
			{
				/// Forbid area is zero.
				Eigen::Vector3d a = V.row(hole[j]);
				Eigen::Vector3d b = V.row(hole[k]);
				Eigen::Vector3d c = V.row(hole[j + i]);
				double area = TriangleArea(a, b, c);
				if(area < 1e-5) continue;

				w = weight[j][k] + weight[k][j + i] + computer_weight(j, k, j + i);
				if(w < wmin){
					wmin = w;
					index_min = k;
				}
			}

			weight[j][j + i] = wmin;
			res[j][j + i] = index_min;
		}
	}

	// add triangle to hole
	std::queue<Eigen::Vector2i> que;
	que.push(Eigen::Vector2i(0, n - 1));
	while(!que.empty()){
		Eigen::Vector2i front = que.front();
		que.pop();
		int start = front[0];
		int end = front[1];
		if(end - start < 2) continue;
		int target = res[start][end];
		patches.emplace_back(hole[start], hole[target], hole[end]);
		que.push({start, target});
		que.push({target, end});
	}
	return 0;
}

bool MESHIO::dynamicFillHole(Mesh &mesh){

	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	
	Eigen::MatrixXi C;
	igl::bfs_orient(T, T, C);

	std::vector<std::array<int, 2>> edge;
	std::vector<std::vector<int>> loop;
	std::vector<std::set<int>> point_to_tris(V.rows());


	for (int i = 0; i < mesh.Topo.rows(); ++i){
		for(int j = 0; j < 3; ++j){
			int a = mesh.Topo(i, j);
			point_to_tris[a].insert(i);
		}
	}
	for (int i = 0; i < mesh.Topo.rows(); ++i){
		for(int j = 0; j < 3; ++j){
			int a = mesh.Topo(i, j);
			int b = mesh.Topo(i, (j + 1) % 3);
			std::vector<int> facets;
			std::set_intersection(point_to_tris[a].begin(), point_to_tris[a].end(), point_to_tris[b].begin(), point_to_tris[b].end(), std::back_inserter(facets));
			if(facets.size() == 1){
				edge.push_back(std::array<int, 2>{a, b});
			}
		}
	}
	get_boundary_loop(edge, loop);

	/// fix boundary orient
	for(int i = 0; i < loop.size(); ++i){ 
		int a = loop[i][0];
		int b = loop[i][1];
		std::vector<int> facets;
		std::set_intersection(point_to_tris[a].begin(), point_to_tris[a].end(), point_to_tris[b].begin(), point_to_tris[b].end(), std::back_inserter(facets));
		int f = 0;
		for(int j = 0; j < 3; ++j){
			int v1 = mesh.Topo(facets[0], j);
			int v2 = mesh.Topo(facets[0], (j + 1) % 3);
			if(v1 == b && v2 == a){
				f = 1;
				break;
			}
		}
		if(!f){
			std::reverse(loop[i].begin(), loop[i].end());
		}
	}

	for(int i = 0; i < loop[0].size(); ++i){
		std::cout << loop[0][i] << " ";
	}
	std::cout << "\n";

	std::cout << "Hole number is " << loop.size() << "\n";

	std::vector<Eigen::Vector3i> all_patches;
	for(int i = 0; i < loop.size(); ++i){
		auto& hole = loop[i];
		std::vector<Eigen::Vector3i> patches;
		dynamicFillHole(V, T, point_to_tris, hole, patches);
		all_patches.insert(all_patches.end(), patches.begin(), patches.end());
	}

	int origin_tri_num = T.rows();
	T.conservativeResize(origin_tri_num + all_patches.size(), 3);
	int index = 0;
	for(int i = origin_tri_num; i < T.rows(); ++i){
		T.row(i) = all_patches[index++];
	}
	return 0;
}