#include "meshAlgorithm.h"
#include "MeshOrient.h"
#include <iostream>

using namespace std;

/**
 * the start point is (start_x, start_y, start_z), the orient is (end_x, end_y, end_z), angle is PI * angle. angle \in (0, 2).
 * @param rotateVec is the param of rotate.
 * @param V
 * @param T
 * @return 1/0
 */
bool MESHIO::rotatePoint(vector<double> rotateVec, Mesh &mesh)
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
bool MESHIO::createBox(std::vector<double> create_box, Mesh &mesh)
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

	for(int i = 0; i < create_box.size(); i += 6)
	{
		Eigen::Vector3d minn, maxx;
		minn = Eigen::Vector3d(create_box[i + 0], create_box[i + 1], create_box[i + 2]).cwiseMin(Eigen::Vector3d(create_box[i + 3], create_box[i + 4], create_box[i + 5]));
		maxx = Eigen::Vector3d(create_box[i + 0], create_box[i + 1], create_box[i + 2]).cwiseMax(Eigen::Vector3d(create_box[i + 3], create_box[i + 4], create_box[i + 5]));

		Eigen::Vector3d V_lst[8];
		
		mesh.Vertex.row(base_v + 0) = maxx;
		mesh.Vertex.row(base_v + 1) = Eigen::Vector3d( minn.x(), maxx.y(), maxx.z() );
		mesh.Vertex.row(base_v + 2) = Eigen::Vector3d( maxx.x(), minn.y(), maxx.z() );
		mesh.Vertex.row(base_v + 3) = Eigen::Vector3d( minn.x(), minn.y(), maxx.z() );
		mesh.Vertex.row(base_v + 4) = Eigen::Vector3d( maxx.x(), maxx.y(), minn.z() );
		mesh.Vertex.row(base_v + 5) = Eigen::Vector3d( minn.x(), maxx.y(), minn.z() );
		mesh.Vertex.row(base_v + 6) = Eigen::Vector3d( maxx.x(), minn.y(), minn.z() );
		mesh.Vertex.row(base_v + 7) = minn;

		for(int j = 0; j < boxTri.size() / 3; j++)
		{
			mesh.Masks(base_t + j, 0) = base_c;
			mesh.Topo.row(base_t + j) = Eigen::Vector3i( base_v + boxTri[j * 3 + 0] - 1,
														 base_v + boxTri[j * 3 + 1] - 1,
														 base_v + boxTri[j * 3 + 2] - 1);
		}
		base_v += 8;
		base_t += 12;
		base_c += 1;
	}
	
}

/**
 * Add bounding box.
 * @param boxVec
 * @param V
 * @param T
 * @return
 */
bool MESHIO::addBox(vector<double> boxVec, Mesh &mesh)
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

	std::cout << midx << " " << midy << " " << midz << '\n';
	std::cout << tmpx[tmpx.size() - 1] - tmpx[0] << " " << tmpy[tmpy.size() - 1] - tmpy[0] << " " << \
		tmpz[tmpz.size() - 1] - tmpz[0] << '\n';

	tmpx.clear(); tmpy.clear(); tmpz.clear();


	Eigen::MatrixXd tmpV;
	tmpV.resize(V.rows() + 8, V.cols());
	for (int i = 0; i < V.rows(); i++)
		for (int j = 0; j < V.cols(); j++)
			tmpV(i, j) = V(i, j);

	Eigen::MatrixXi tmpT;
	tmpT.resize(T.rows() + 12, T.cols());
	for (int i = 0; i < T.rows(); i++)
		for (int j = 0; j < T.cols(); j++)
			tmpT(i, j) = T(i, j);

	Eigen::MatrixXi tmpM;
	tmpM.resize(M.rows() + 12, M.cols());
	for (int i = 0; i < tmpM.cols(); i++)
		tmpM(i, 0) = M(i, 0);


	int V_index = V.rows();
	for (double k = hight / 2; k > -hight; k -= hight) // up and down
	{
		for (double j = width / 2; j > -width; j -= width) // left and right
		{
			for (double i = len / 2; i > -len; i -= len) // forward and back
			{
				tmpV(V_index, 0) = midx + i;
				tmpV(V_index, 1) = midy + j;
				tmpV(V_index, 2) = midz + k;
				V_index++;
			}
		}
	}


	/***************************
	 *    4 _________2
	 *     /|        /|
	 *   3/_|______1/ |
	 *    | |8_____|__|6
	 *    | /      | /
	 *  7 |/_______|/5
	 ****************************/

	V_index = V.rows() - 1;
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

	for (int i = T.rows(), j = 0; i < T.rows() + 12; i++, j += 3) {
		tmpT(i, 0) = V_index + boxTri[j + 0];
		tmpT(i, 1) = V_index + boxTri[j + 1];
		tmpT(i, 2) = V_index + boxTri[j + 2];
		tmpM(i, 0) = M(T.rows() - 1, 0) + 1;
	}


	V.resize(tmpV.rows(), tmpV.cols());
	T.resize(tmpT.rows(), tmpT.cols());
	V = tmpV;
	T = tmpT;

	std::cout << "Boxed\n";

	return 1;
}

bool MESHIO::reverseOrient(Eigen::MatrixXi &T) {
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

bool MESHIO::repair(Mesh &mesh)
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

bool MESHIO::resetOrientation(Mesh &mesh, bool reset_mask) {
	vector<vector<double>> point_list(mesh.Vertex.rows(), vector<double>(mesh.Vertex.cols(), 0));
	vector<vector<int>> facet_list(mesh.Topo.rows(), vector<int>(mesh.Topo.cols(), 0));
	for(int i = 0; i < mesh.Vertex.rows(); i++) {
		for(int j = 0; j < mesh.Vertex.cols(); j++) {
			point_list[i][j] = mesh.Vertex(i, j);
		}
	}
	for(int i = 0; i < mesh.Topo.rows(); i++) {
		for(int j = 0; j < mesh.Topo.cols(); j++) {
			facet_list[i][j] = mesh.Topo(i, j);
		}
	}
	vector<int> block_mark;
	TIGER::resetOrientation(point_list, facet_list, block_mark);
	for(int i = 0; i < mesh.Vertex.rows(); i++) {
		for(int j = 0; j < mesh.Vertex.cols(); j++) {
			mesh.Vertex(i, j) = point_list[i][j];
		}
	}
	for(int i = 0; i < mesh.Topo.rows(); i++) {
		for(int j = 0; j < mesh.Topo.cols(); j++) {
			mesh.Topo(i, j) = facet_list[i][j];
		}
	}
	if (reset_mask) {
			for (int i = 0; i < mesh.Topo.rows(); i++) {
				mesh.Masks(i,0) = block_mark[i];
			}
	}
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
	auto rule = [&in_vec, &eps](int i, int j)->bool{
		if( fabs(in_vec(i, 0) - in_vec(j, 0)) > eps) return in_vec(i, 0) < in_vec(j, 0);
		if( fabs(in_vec(i, 1) - in_vec(j, 1)) > eps) return in_vec(i, 1) < in_vec(j, 1);
		if( fabs(in_vec(i, 2) - in_vec(j, 2)) > eps) return in_vec(i, 2) < in_vec(j, 2);
		return i < j;
	};
	std::sort(ind.data(), ind.data() + ind.size(), rule);
	out_vec.resize(in_vec.rows(), in_vec.cols());
	for(int i = 0; i < in_vec.rows(); i++)
	{
		out_vec.row(i) << in_vec.row(ind(i));
	}
}

/**
 * Remove dulplicate point.
 * @param V
 * @param T
 */
bool MESHIO::removeDulplicatePoint(Eigen::MatrixXd& V, Eigen::MatrixXi& F, double eps){
	Eigen::MatrixXd V_in = V;
	Eigen::MatrixXi F_in = F;
	Eigen::MatrixXd V_sort;
	Eigen::VectorXi ind;
	eigen_sort_3d(V_in, V_sort, ind, eps);
	int cnt = 0;
	std::unordered_map<int, int> mpid;
	std::vector<int> new_point_lst;
	for(int i = 0; i < V_sort.rows(); i++)
	{
		while (i + 1 < V_sort.rows())
		{
			if( (V_sort.row(i + 1) - V_sort.row(i)).norm() < eps )
			{
				mpid[ ind[i] ] = cnt;
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

	for(int i = 0; i < new_point_lst.size(); i++)
	{
		V.row(i) = V_in.row(new_point_lst[i]);
	}

	for(int i = 0; i < F_in.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			F(i, j) = mpid[ F_in(i, j) ];
		}
	}

	std::cout << "Remove " << V_in.rows() - V.rows() << " dulplicate points\n";

}