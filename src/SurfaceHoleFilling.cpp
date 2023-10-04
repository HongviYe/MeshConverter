#include "SurfaceHoleFilling.h"
#include "meshIO.h"
#include "SurfaceRemesh.h"

#include "igl/boundary_loop.h"
#include "igl/cat.h"
#include "igl/slice.h"
#include "igl/upsample.h"
#include "igl/exact_geodesic.h"
#include "igl/all_pairs_distances.h"
#include "igl/per_vertex_normals.h"
#include "igl/cotmatrix.h"
#include "igl/grad.h"
#include "igl/harmonic.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Geometry"
#include <map>
#include <array>

using namespace Eigen;

SurfaceHoleFilling::SurfaceHoleFilling(Mesh& xyz_mesh) : xyz_mesh_(xyz_mesh){}


template <typename Derived, class MatC>
void SurfaceHoleFilling::my_cat(
	const Eigen::MatrixBase<Derived> & A,
	const Eigen::MatrixBase<Derived> & B,
	MatC & C)
{
	if (A.size() == 0)
	{
		C = B;
		return;
	}
	if (B.size() == 0)
	{
		C = A;
		return;
	}
	assert(A.cols() == B.cols());
	C.resize(A.rows() + B.rows(), A.cols());
	C << A, B;
}

template <
	typename DerivedX,
	typename DerivedR,
	typename DerivedC,
	typename DerivedY>
	void SurfaceHoleFilling::my_slice(
		const Eigen::DenseBase<DerivedX> & X,
		const Eigen::DenseBase<DerivedR> & R,//行索引矩阵
		const Eigen::DenseBase<DerivedC> & C,//列索引矩阵
		Eigen::PlainObjectBase<DerivedY> & Y)
{
	int x = R.size();
	int y = C.size();
	if (x == 0 || y == 0)
	{
		Y.resize(x, y);
		return;
	}

	int m = X.rows();
	int n = X.cols();
	assert(R.minCoeff() >= 0 &&R.maxCoeff ()< m);
	assert(C.minCoeff() >= 0 &&C.maxCoeff ()< n);
	Y.resize(x, y);
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			Y(i, j) = X(R(i, 0), C(j, 0));
		}
	}
}

void SurfaceHoleFilling::my_triangle_triangle_adjacency(
		const Eigen::MatrixXi& F,
		Eigen::MatrixXi& TT,
		Eigen::MatrixXi& TTi)
{
	std::map<std::pair<int, int>, std::vector<int>> edge_triangle;
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int v = F(i, j), vn = F(i, (j + 1) % 3);
			if (v > vn)std::swap(v, vn);
			std::pair<int, int>p(v, vn);
			edge_triangle[p].push_back(i);
		}
	}
	TT=MatrixXi::Constant(F.rows(), 3, -1);
	for (int i = 0; i < F.rows(); i++)//获得三角形面片边的相邻三角形
	{
		for (int j = 0; j < 3; j++)
		{
			int v = F(i, j), vn = F(i, (j + 1) % 3);
			if (v > vn)std::swap(v, vn);
			std::pair<int, int>p(v, vn);
			for (int k = 0; k < edge_triangle[p].size(); k++)
			{
				int fn = edge_triangle[p][k];
				if (fn == i)continue;
				TT(i, j) = fn;
				break;
			}
		}
	}
	TTi = MatrixXi::Constant(TT.rows(), TT.cols(), -1);//三角形面片i和相邻三角形fn的相邻边所在fn中的边序号
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int v = F(i, j), vn = F(i, (j + 1) % 3);
			int fn = TT(i, j);
			if (fn >= 0)
			{
				for (int k = 0; k < 3; k++)
				{
					int vk = F(fn, k), vkn = F(fn, (k + 1) % 3);
					if (v == vkn && vn == vk)
					{
						TTi(i, j) = k;
						break;
					}
				}
			}
		}
	}

}


//template <typename DerivedF, typename DerivedTT, typename DerivedTTi>
//void SurfaceHoleFilling::my_triangle_triangle_adjacency(
//	const Eigen::MatrixBase<DerivedF>& F,
//	Eigen::PlainObjectBase<DerivedTT>& TT,
//	Eigen::PlainObjectBase<DerivedTTi>& TTi)
//{
//	/*const int n = F.maxCoeff() + 1;
//	typedef Eigen::Matrix<typename DerivedTT::Scalar, Eigen::Dynamic, 1> VectorXI;
//	VectorXI VF, NI;
//	my_vertex_triangle_adjacency(F, n, VF, NI);//每个点获得该点的三角形*/
//	std::map<std::pair<int, int>, std::vector<int>> edge_triangle;
//	for (int i = 0; i < F.rows(); i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			int v = F(i, j), vn = F(i, (j + 1) % 3);
//			if (v > vn)std::swap(v, vn);
//			std::pair<int, int>p(v, vn);
//			edge_triangle[p].push_back(i);
//		}
//	}
//	TT = DerivedTT::Constant(F.rows(), 3, -1);
//	for (int i = 0; i < F.rows(); i++)//获得三角形面片边的相邻三角形
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			int v = F(i, j), vn = F(i, (j + 1) % 3);
//			if (v > vn)std::swap(v, vn);
//			std::pair<int, int>p(v, vn);
//			for (int k = 0; k < edge_triangle[p].size(); k++)
//			{
//				int fn = edge_triangle[p][k];
//				if (fn == i)continue;
//				TT(i, j) = fn;
//				break;
//			}
//			/*for (int k = NI[v]; k < NI[v + 1]; k++)
//			{
//				int fn = VF[k];
//				if (fn == i)continue;//不是当前面片
//				if (F(fn, 0) == vn || F(fn, 1) == vn || F(fn, 2) == vn)
//				{
//					TT(i, j) = fn;
//					break;
//				}
//			}*/
//		}
//	}
//	TTi = DerivedTTi::Constant(TT.rows(), TT.cols(), -1);//三角形面片i和相邻三角形fn的相邻边所在fn中的边序号
//	for (int i = 0; i < F.rows(); i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			int v = F(i, j),vn = F(i, (j + 1) % 3);
//			int fn = TT(i, j);
//			if (fn >= 0)
//			{
//				for (int k = 0; k < 3; k++)
//				{
//					int vk = F(fn, k), vkn = F(fn, (k + 1) % 3);
//					if (v == vkn && vn == vk)
//					{
//						TTi(i, j) = k;
//						break;
//					}
//				}
//			}
//		}
//	}
//}


void SurfaceHoleFilling::my_boundary_loop(
	const Eigen::MatrixXi& F,
	Eigen::VectorXi& L)
{
	if (F.rows() == 0)
		return;
	Eigen::MatrixXi TT, TTi;

	my_triangle_triangle_adjacency(F, TT, TTi);
	std::vector<std::vector<int> > VF, VFi;//VF是点所在面片的集合，VFi是点所在边的集合
	VF.clear(); VFi.clear();
	VF.resize(F.maxCoeff() + 1); VFi.resize(F.maxCoeff() + 1);
	//typedef typename DerivedF::Index Index;
	for (int i = 0; i < F.rows(); ++i)
	{
		for (int j = 0; j < F.cols(); ++j)
		{
			VF[F(i, j)].push_back(i);
			VFi[F(i, j)].push_back(j);
		}
	}
	std::vector<bool> unvisited_border_vertex(F.maxCoeff() + 1);
	std::set<int>rest_vertex;
	for (int i = 0; i < unvisited_border_vertex.size(); i++)
	{
		unvisited_border_vertex[i] = false;
	}
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			if (TT(i, j) == -1)
			{
				unvisited_border_vertex[F(i, j)] = true;
				unvisited_border_vertex[F(i, (j + 1) % F.cols())] = true;
			}
		}
	}
	for (int i = 0; i < unvisited_border_vertex.size(); i++)
	{
		if (unvisited_border_vertex[i])
			rest_vertex.insert(rest_vertex.end(), i);
	}
	std::vector<std::vector<int>>all_loop;
	while (rest_vertex.size())
	{
		std::vector<int>loop;
		auto start = *rest_vertex.begin();
		rest_vertex.erase(rest_vertex.begin());
		unvisited_border_vertex[start] = false;
		loop.push_back(start);
		while (1)
		{
			bool newedge = false;
			int v = loop[loop.size() - 1];
			int next;
			for (int i = 0; i < VF[v].size() && !newedge; i++)
			{
				int face = VF[v][i];
				if (TT.row(face).minCoeff() < 0)
				{
					int vLoc = -1;
					if (F(face, 0) == v) vLoc = 0;
					if (F(face, 1) == v) vLoc = 1;
					if (F(face, 2) == v) vLoc = 2;
					int vNext = F(face, (vLoc + 1) % F.cols());
					if (unvisited_border_vertex[vNext] && TT(face, vLoc) < 0)
					{
						next = vNext;
						newedge = true;
					}
				}
			}
			if (newedge)
			{
				loop.push_back(next);
				rest_vertex.erase(next);
				unvisited_border_vertex[next] = false;
			}
			else
				break;
		}
		all_loop.push_back(loop);
	}

	int max_loop_idx = -1;
	int max_len = 0;
	for (int i = 0; i < all_loop.size(); i++)
	{
		if (all_loop[i].size() > max_len)
		{
			max_loop_idx = i;
			max_len = all_loop[i].size();
		}
	}
	if (max_loop_idx == -1)
	{
		return;
	}
	L.resize(all_loop[max_loop_idx].size());
	for (int i = 0; i < all_loop[max_loop_idx].size(); i++)
	{
		L[i] = all_loop[max_loop_idx][i];
	}
}


//template <typename DerivedF, typename DerivedL>
//void SurfaceHoleFilling::my_boundary_loop(
//	const Eigen::MatrixBase<DerivedF> & F,
//	Eigen::PlainObjectBase<DerivedL>& L)
//{
//	using namespace std;
//	using namespace Eigen;
//
//	if (F.rows() == 0)
//		return;
//	Eigen::Matrix<typename DerivedF::Scalar, Eigen::Dynamic, Eigen::Dynamic> TT, TTi;
//
//	my_triangle_triangle_adjacency(F, TT, TTi);
//	vector<std::vector<int> > VF, VFi;//VF是点所在面片的集合，VFi是点所在边的集合
//	VF.clear(); VFi.clear(); 
//	VF.resize(F.maxCoeff() + 1); VFi.resize(F.maxCoeff() + 1);
//	//typedef typename DerivedF::Index Index;
//	for (int i = 0; i < F.rows(); ++i)
//	{
//		for (int j = 0; j < F.cols(); ++j)
//		{
//			VF[F(i, j)].push_back(i);
//			VFi[F(i, j)].push_back(j);
//		}
//	}
//	std::vector<bool> unvisited_border_vertex(F.maxCoeff() + 1);
//	set<int>rest_vertex;
//	for (int i = 0; i < unvisited_border_vertex.size(); i++)
//	{
//		unvisited_border_vertex[i] = false;
//	}
//	for (int i = 0; i < F.rows(); i++)
//	{
//		for (int j = 0; j < F.cols(); j++)
//		{
//			if (TT(i, j) == -1)
//			{
//				unvisited_border_vertex[F(i, j)] = true;
//				unvisited_border_vertex[F(i, (j + 1) % F.cols())] = true;
//			}
//		}
//	}
//	for (int i = 0; i < unvisited_border_vertex.size(); i++)
//	{
//		if(unvisited_border_vertex[i])
//			rest_vertex.insert(rest_vertex.end(), i);
//	}
//	vector<vector<int>>all_loop;
//	while (rest_vertex.size())
//	{
//		vector<int>loop;
//		auto start = *rest_vertex.begin();
//		rest_vertex.erase(rest_vertex.begin());
//		unvisited_border_vertex[start] = false;
//		loop.push_back(start);
//		while (1)
//		{
//			bool newedge = false;
//			int v = loop[loop.size() - 1];
//			int next;
//			for (int i = 0; i < VF[v].size()&&!newedge; i++)
//			{
//				int face = VF[v][i];
//				if (TT.row(face).minCoeff() < 0)
//				{
//					int vLoc = -1;
//					if (F(face, 0) == v) vLoc = 0;
//					if (F(face, 1) == v) vLoc = 1;
//					if (F(face, 2) == v) vLoc = 2;
//					int vNext = F(face, (vLoc + 1) % F.cols());
//					if (unvisited_border_vertex[vNext] && TT(face, vLoc) < 0)
//					{
//						next = vNext;
//						newedge = true;
//					}
//				}
//			}
//			if (newedge)
//			{
//				loop.push_back(next);
//				rest_vertex.erase(next);
//				unvisited_border_vertex[next] = false;
//			}
//			else
//				break;
//		}
//		all_loop.push_back(loop);
//	}
//
//	int max_loop_idx = -1;
//	int max_len = 0;
//	for (int i = 0; i < all_loop.size(); i++)
//	{
//		if (all_loop[i].size() > max_len)
//		{
//			max_loop_idx = i;
//			max_len = all_loop[i].size();
//		}
//	}
//	if (max_loop_idx == -1)
//	{
//		return;
//	}
//	L.resize(all_loop[max_loop_idx].size());
//	for (int i = 0; i < all_loop[max_loop_idx].size(); i++)
//	{
//		L[i] = all_loop[max_loop_idx][i];
//	}
//}


void SurfaceHoleFilling::my_cotmatrix(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	Eigen::SparseMatrix<double>& L)
{
	L.resize(V.rows(), V.rows());
	Eigen::Matrix<int, Eigen::Dynamic, 2> edges;
	assert(F.cols() == 3);
	L.reserve(10 * V.rows());
	edges.resize(3, 2);
	edges <<
		1, 2,
		2, 0,
		0, 1;
	Eigen::MatrixXd C(F.rows(), F.cols()) ;
	//cotmatrix_entries(V, F, C);
	Eigen::MatrixXd squared_len(F.rows(), 3) ;//边长的平方
	for (int i = 0; i < F.rows(); i++)
	{
		squared_len(i, 0) = (V.row(F(i, 1)) - V.row(F(i, 2))).squaredNorm();
		squared_len(i, 1) = (V.row(F(i, 2)) - V.row(F(i, 0))).squaredNorm();
		squared_len(i, 2) = (V.row(F(i, 0)) - V.row(F(i, 1))).squaredNorm();
	}

	Eigen::MatrixXd len(F.rows(), 3) ;
	len = squared_len.array().sqrt();

	Eigen::MatrixXd double_area(F.rows(), 1);
	//igl::doublearea(len,0.,double_area);

	for (int i = 0; i < F.rows(); i++)
	{
		double arg = (len(i, 0) + len(i, 1) + len(i, 2))*
			(len(i, 2) + len(i, 1) - len(i, 0))*
			(len(i, 2) + len(i, 0) - len(i, 1))*
			(len(i, 0) + len(i, 1) - len(i, 2));
		double_area(i) = 2*0.25*sqrt(arg);
		/*if(double_area(i)!=double_area(i))
		{
			double_area(i)=0.0;
		}*/
	}

	for (int i = 0; i < F.rows(); i++)
	{
		C(i, 0) = (squared_len(i, 1) + squared_len(i, 2) - squared_len(i, 0)) / double_area(i) / 4.0;
		C(i, 1) = (squared_len(i, 2) + squared_len(i, 0) - squared_len(i, 1)) / double_area(i) / 4.0;
		C(i, 2) = (squared_len(i, 0) + squared_len(i, 1) - squared_len(i, 2)) / double_area(i) / 4.0;
	}

	std::vector<Eigen::Triplet<double> > IJV;
	IJV.reserve(F.rows()*edges.rows() * 4);
	for (int i = 0; i < F.rows(); i++)
	{
		// loop over edges of element
		for (int e = 0; e < edges.rows(); e++)
		{
			int source = F(i, edges(e, 0));
			int dest = F(i, edges(e, 1));
			IJV.push_back(Eigen::Triplet<double>(source, dest, C(i, e)));
			IJV.push_back(Eigen::Triplet<double>(dest, source, C(i, e)));
			IJV.push_back(Eigen::Triplet<double>(source, source, -C(i, e)));
			IJV.push_back(Eigen::Triplet<double>(dest, dest, -C(i, e)));
		}
	}
	L.setFromTriplets(IJV.begin(), IJV.end());
}


bool SurfaceHoleFilling::my_harmonic(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		const Eigen::MatrixXi & b,
		const Eigen::MatrixXd & bc,
		const int k,
		Eigen::MatrixXd & W)
{
	//std::cout<<"use my_harmonic"<<std::endl;
	Eigen::SparseMatrix<double> L;
	Eigen::SparseMatrix<double,RowMajor> Lx;

	int n=V.rows();

	my_cotmatrix(V, F, L);
	//igl::cotmatrix(V,F,L);
	//std::cout<<L1-L;

	// Lx = L;
	// for (int i = 0; i < k; i++)
	// {
	// 	Lx = Lx * L;
	// }

	Eigen::SparseMatrix<double> M;

	if(k>1)
	{
		//igl::massmatrix(V,F,igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT,M);
		Eigen::MatrixXd squared_len(F.rows(), 3) ;//边长的平方
		for (int i = 0; i < F.rows(); i++)
		{
		squared_len(i, 0) = (V.row(F(i, 1)) - V.row(F(i, 2))).squaredNorm();
		squared_len(i, 1) = (V.row(F(i, 2)) - V.row(F(i, 0))).squaredNorm();
		squared_len(i, 2) = (V.row(F(i, 0)) - V.row(F(i, 1))).squaredNorm();
		}

		Eigen::MatrixXd len(F.rows(), 3) ;
		len = squared_len.array().sqrt();

		Eigen::MatrixXd double_area(F.rows(), 1);
		//igl::doublearea(len,0.,double_area);
		for (int i = 0; i < F.rows(); i++)
		{
			double arg = (len(i, 0) + len(i, 1) + len(i, 2))*
				(len(i, 2) + len(i, 1) - len(i, 0))*
				(len(i, 2) + len(i, 0) - len(i, 1))*
				(len(i, 0) + len(i, 1) - len(i, 2));
			double_area(i) = 2*0.25*sqrt(arg);
			if(double_area(i)!=double_area(i))
			{
				double_area(i)=0.0;
			}
		}
		Matrix<int,Dynamic,1> MI;
		Matrix<int,Dynamic,1> MJ;
		Matrix<double,Dynamic,1> MV;

		MI.resize(F.rows()*3,1); MJ.resize(F.rows()*3,1); MV.resize(F.rows()*3,1);
		MI.block(0*F.rows(),0,F.rows(),1) = F.col(0);
		MI.block(1*F.rows(),0,F.rows(),1) = F.col(1);
		MI.block(2*F.rows(),0,F.rows(),1) = F.col(2);
		MJ = MI;

		Matrix<double,Dynamic,3> cosines(F.rows(),3);//每个三角形的余弦值
		cosines.col(0) = 
		  (len.col(2).array().pow(2)+len.col(1).array().pow(2)-len.col(0).array().pow(2))/(len.col(1).array()*len.col(2).array()*2.0);
		cosines.col(1) = 
		  (len.col(0).array().pow(2)+len.col(2).array().pow(2)-len.col(1).array().pow(2))/(len.col(2).array()*len.col(0).array()*2.0);
		cosines.col(2) = 
		  (len.col(1).array().pow(2)+len.col(0).array().pow(2)-len.col(2).array().pow(2))/(len.col(0).array()*len.col(1).array()*2.0);
		
		Matrix<double,Dynamic,3> barycentric = cosines.array() * len.array();//cosa*A
		
		barycentric = (barycentric.array().colwise() / barycentric.rowwise().sum().array()).eval();//?

		Matrix<double,Dynamic,3> partial = barycentric;
		partial.col(0).array() *= double_area.array() * 0.5;
		partial.col(1).array() *= double_area.array() * 0.5;
		partial.col(2).array() *= double_area.array() * 0.5;
		Matrix<double,Dynamic,3> quads(partial.rows(),partial.cols());
		quads.col(0) = (partial.col(1)+partial.col(2))*0.5;
		quads.col(1) = (partial.col(2)+partial.col(0))*0.5;
		quads.col(2) = (partial.col(0)+partial.col(1))*0.5;

		quads.col(0) = (cosines.col(0).array()<0).select( 0.25*double_area,quads.col(0));
		quads.col(1) = (cosines.col(0).array()<0).select(0.125*double_area,quads.col(1));
		quads.col(2) = (cosines.col(0).array()<0).select(0.125*double_area,quads.col(2));

		quads.col(0) = (cosines.col(1).array()<0).select(0.125*double_area,quads.col(0));
		quads.col(1) = (cosines.col(1).array()<0).select(0.25*double_area,quads.col(1));
		quads.col(2) = (cosines.col(1).array()<0).select(0.125*double_area,quads.col(2));

		quads.col(0) = (cosines.col(2).array()<0).select(0.125*double_area,quads.col(0));
		quads.col(1) = (cosines.col(2).array()<0).select(0.125*double_area,quads.col(1));
		quads.col(2) = (cosines.col(2).array()<0).select( 0.25*double_area,quads.col(2));

		MV.block(0*F.rows(),0,F.rows(),1) = quads.col(0);
		MV.block(1*F.rows(),0,F.rows(),1) = quads.col(1);
		MV.block(2*F.rows(),0,F.rows(),1) = quads.col(2);

		std::vector<Eigen::Triplet<double> > triplist;
		triplist.reserve(MI.size());
		for(int x = 0;x<MI.size();x++)
		{
			triplist.push_back(Eigen::Triplet<double>(MI(x),MJ(x),MV(x)));
		}
		M.resize(n,n);
		M.setFromTriplets(triplist.begin(),triplist.end());
	}

	Eigen::SparseMatrix<double> Q;
	Q = L;
	Eigen::SparseMatrix<double> Mi;
	//igl::invert_diag(M,Mi);
	Mi=M;
	for(int k=0; k<Mi.outerSize(); ++k)
	{
	// Iterate over inside
		for(Eigen::SparseMatrix<double>::InnerIterator it(Mi,k); it; ++it)
		{
		 if(it.col() == it.row())
		 {
			double v = it.value();
			assert(v != 0);
			v = ((double)1.0)/v;
			Mi.coeffRef(it.row(),it.col()) = v;
		 }
		}
	}

	for(int p = 1;p<k;p++)
	{
		Q = (Q*Mi*L).eval();
	}
	

	Lx = Q;
	W.resize(V.rows(),3);
	Eigen::SparseMatrix<double> I;
	I.resize(V.rows(), V.rows());
	for(int i=0;i<I.rows();i++)
	{
		I.insert(i,i)=1.0;
	}
	/*std::vector<Eigen::Triplet<double> > IJV;
	IJV.reserve(V.rows()*3);
	for(int i=0;i<b.size();i++)
	{
		IJV.push_back(Eigen::Triplet<double>(b(i),b(i),1.));
	}
	I.setFromTriplets(IJV.begin(), IJV.end());*/

	for(int i=0;i<b.size();i++)
	{
		Lx.row(b(i))=I.row(b(i));
	}

	Eigen::MatrixXd res(Lx.rows(),3);
	res.setZero(Lx.rows(),3);

	SparseLU<SparseMatrix<double>>  solver;
	solver.compute(Lx); 

	//SparseQR<SparseMatrix<double>,Eigen::NaturalOrdering<int>>  solver2;
	//solver2.compute(Lx);
	//std::cout<<"us QR";

	//solver.analyzePattern(Lx);
	//solver.factorize(Lx);

	//Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver1;
	//solver1.compute(Lx);
	//choleskyA.factorize(Lx);
	
	
	//choleskyA.compute(Lx);

	//Eigen::SimplicialCholesky<SparseMatrix<double>> chol(Lx); 
	//Eigen::VectorXd x = chol.solve(b);
	
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<b.size();j++)
		{
			res(b(j),i)=bc(j,i);
		}
		
		//W.col(i)=solver.solve(res.col(i));
	}
	//std::cout<<res<<"\n";
	Eigen::VectorXd x = solver.solve(res.col(0));
	Eigen::VectorXd y = solver.solve(res.col(1));
	Eigen::VectorXd z = solver.solve(res.col(2));
	//W = solver.solve(res);
	W.col(0)=x;
	W.col(1)=y;
	W.col(2)=z;
	//std::cout<<W-V<<"\n";
	//std::cout<<W<<std::endl<<"LU"<<std::endl;

}

void SurfaceHoleFilling::init_hole_fill()
{
	//std::cout << "init_hole\n";
	int upsampleN = 4;
	MatrixXd& originalV = this->xyz_mesh_.Vertex;
	MatrixXi& originalF = this->xyz_mesh_.Topo;

	//igl::boundary_loop(originalF, bnd_loop_);
	my_boundary_loop(originalF, bnd_loop_);
	if (bnd_loop_.size() == 0) {
		printf("Mesh has no hole!");
	exit(0);
	}

	int v0 = bnd_loop_[0];
	int v1 = bnd_loop_[1];
	int f = 0;
	for(int i = 0; i < originalF.rows(); ++i)
	{
		for(int j = 0; j < 3; ++j)
		{
			if(originalF(i, j) == v0 && originalF(i, (j + 1) % 3) == v1)
			{
				f = 1;
				break;
			}
		}
		if(f) break;
	}

	if(f)
	{
		std::cout << "Reverse normal\n";
		std::vector<int> bnd_loop_vec;
		for (int i = 0; i < bnd_loop_.size(); ++i)
			bnd_loop_vec.push_back(bnd_loop_[i]);
		std::reverse(bnd_loop_vec.begin(), bnd_loop_vec.end());
		for (int i = 0; i < bnd_loop_vec.size(); ++i)
			bnd_loop_[i] = bnd_loop_vec[i];
	}


	

	// compute boundary center.
	VectorXd bcenter(3);
	VectorXi R = bnd_loop_;
	VectorXi C(3);
	C << 0, 1, 2;
	MatrixXd B;
	MatrixXd A = this->xyz_mesh_.Vertex;
	//igl::slice(A, R, C, B);
	my_slice(A, R, C, B);
	bcenter = (1.0f / bnd_loop_.size()) * B.colwise().sum();

	hole_mesh_.Vertex = B;
	hole_mesh_.Topo = MatrixXi(bnd_loop_.size(), 3);
	for(int i = 0; i < bnd_loop_.size(); ++i)
	{
		hole_mesh_.Topo(i, 0) = hole_mesh_.Vertex.rows();
		hole_mesh_.Topo(i, 1) = i;
		hole_mesh_.Topo(i, 2) = (i + 1) % bnd_loop_.size();
	}
	hole_mesh_.Vertex.conservativeResize(B.rows() + 1, 3);
	hole_mesh_.Vertex.row(hole_mesh_.Vertex.rows() - 1) = bcenter;

	SurfaceRemesh::iso_remesh(hole_mesh_);
	//MESHIO::remesh(hole_mesh_);


	// MatrixXd all_V;
	// MatrixXi all_F;

	// MatrixXd new_v = MatrixXd(1, 3);
	// new_v.row(0) = bcenter;
	// igl::cat(1, originalV, new_v, all_V);

	// MatrixXi new_f = MatrixXi(bnd_loop_.size(), 3);
	// for(int i = 0; i < bnd_loop_.size(); ++i)
	// {
	// 	new_f(i, 0) = (int)all_V.rows() - 1;
	// 	new_f(i, 1) = bnd_loop_(i, 0);
	// 	new_f(i, 2) = bnd_loop_((i + 1) % bnd_loop_.size(), 0);
	// }

	// igl::cat(1, originalF, new_f, all_F);
	
	// this->xyz_mesh_.vertex = all_V;
	// this->xyz_mesh_.topo = all_F;

}

const Mesh& SurfaceHoleFilling::get_all_mesh()
{
	return this->all_mesh_;
}

const Eigen::MatrixXd& SurfaceHoleFilling::get_default_guide_h()
{
	return guide_h_;
}

const Mesh& SurfaceHoleFilling::get_xyz_mesh()
{
	return xyz_mesh_;
}

const Mesh& SurfaceHoleFilling::get_hole_mesh()
{
	return hole_mesh_;
}







//void SurfaceHoleFilling::cal_default_guide_h()
//{
//	std::cout << "Start -- Get guild vector field by line dis.\n";
//	this->guide_h_.resize(hole_mesh_.Topo.rows(), 3);
//
//	Eigen::MatrixXd all_normal;
//	igl::per_vertex_normals(this->xyz_mesh_.Vertex, this->xyz_mesh_.Topo, all_normal);
//	this->xyz_mesh_.v_normal = all_normal;
//
//	Eigen::MatrixXd bnd_normal;
//	my_slice(all_normal, this->bnd_loop_, Eigen::Vector3d(0, 1, 2), bnd_normal);
//	//igl::slice(all_normal, this->bnd_loop_, vec3d(0, 1, 2), bnd_normal);
//
//	const auto get_distance = [&](const int fid)
//	{
//		// Eigen::VectorXi VS, FS, VT, FT;
//		// FS.resize(1);
//		// FS << fid;
//		// VT.resizeLike(this->bnd_loop_);
//		// VT.setLinSpaced(this->bnd_loop_.size(), 0, this->bnd_loop_.size() - 1);
//		// Eigen::VectorXd d;
//		// igl::exact_geodesic(hole_mesh_.vertex, hole_mesh_.topo, VS, FS, VT, FT, d);
//
//		Eigen::VectorXd d;
//		Eigen::Vector3d center(0, 0, 0);
//		for(int i = 0; i < 3; ++i)
//		{
//			center += hole_mesh_.Vertex.row(hole_mesh_.Topo(fid, i));
//		}
//		center /= 3.0;
//
//		d.resize(this->bnd_loop_.size());
//		for(int i = 0; i < this->bnd_loop_.size(); ++i)
//		{
//			d[i] = (center - Eigen::Vector3d(this->hole_mesh_.Vertex.row(i)) ).norm();
//		}
//		return d;
//	};
//
//	for(int i = 0; i < hole_mesh_.Topo.rows(); i++)
//	{
//		auto dis = get_distance(i);
//		double all_dis = dis.sum();
//		// double re_all_dis = 1.0 / all_dis;
//
//		Eigen::Vector3d ave_normal(0, 0, 0);
//  	for(int j = 0; j < dis.size(); ++j)
//		{
//			// ave_normal += dis[j] * re_all_dis * bnd_normal.row(j);
//			ave_normal += all_dis / dis[j] * bnd_normal.row(j);
//		}
//		this->guide_h_.row(i) = ave_normal.normalized();
//	}
//	// this->guide_h_.col(2) *= -1;
//	this->hole_mesh_.f_normal = this->guide_h_;
//	std::cout << "End -- Get guild vector field by line dis.\n";
//
//}
//
//void SurfaceHoleFilling::possion_deformation(const Eigen::MatrixXd& guide_h)
//{
//	if(guide_h.rows() == 0)
//	{
//		std::cout << "guide vector field is not assign.\n";
//		return ;
//	}
//	// A * x = B
//	// using T = typename Eigen::Triplet<double>;
//	// using SMatrix = typename Eigen::SparseMatrix<double>;
//	int n = this->hole_mesh_.Vertex.rows();
//	int bnd_num = this->bnd_loop_.size();
//	
//	/// set A matrix
//	Eigen::SparseMatrix<double> cot_matrix, A(n, n);
//	my_cotmatrix(this->hole_mesh_.Vertex, this->hole_mesh_.Topo, cot_matrix);
//	std::vector<Eigen::Triplet<double>> triplets;
//
//	for(int i = 0; i < bnd_num; ++i)
//	{
//		triplets.emplace_back(i, i, 1.0);
//	}
//
//	for(int i = bnd_num; i < n; i++)
//	{
//		for(Eigen::SparseMatrix<double>::InnerIterator iter(cot_matrix, i); iter; ++iter)
//		{
//			int v = iter.index();
//			double w = iter.value();
//			triplets.emplace_back(i, v, w);
//			// std::cout << i << " " << v << " " << w << "\n";
//		}
//	}
//	A.setFromTriplets(triplets.begin(), triplets.end());
//
//	// set b matrix
//	Eigen::MatrixXd B(n, 3);
//	B.setZero();
//	for(int i = 0; i < bnd_num; ++i)
//	{
//		B.row(i) = this->hole_mesh_.Vertex.row(i);
//	}
//
//	// //// for test :
//	// std::vector<vec3d> r_v_vec;
//	// std::vector<vec3i> r_t_vec;
//	// int cnt = 0;
//	// ////
//
//	for(int i = 0; i < this->hole_mesh_.Topo.rows(); ++i)
//	{
//		Eigen::Vector3i vidxs = this->hole_mesh_.Topo.row(i);
//		std::array<Eigen::Vector3d, 3> v_array =
//				{this->hole_mesh_.Vertex.row(vidxs[0]),
//				 this->hole_mesh_.Vertex.row(vidxs[1]),
//				 this->hole_mesh_.Vertex.row(vidxs[2])};
//		Eigen::Vector3d center_coord = (v_array[0] + v_array[1] + v_array[2]) / 3.0;
//		
//		const Eigen::Vector3d target_normal = guide_h.row(i);
//		const Eigen::Vector3d origin_normal = (v_array[1] - v_array[0]).cross(v_array[2] - v_array[0]).normalized();
//
//		Eigen::Matrix4d to_zero = Eigen::Matrix4d::Identity();
//		to_zero.col(3).topRows(3) = -center_coord;
//		Eigen::Matrix4d re_to_zero = Eigen::Matrix4d::Identity();
//		re_to_zero.col(3).topRows(3) = center_coord;
//
//		Eigen::Matrix3d rotate_matrix_33 = Eigen::Quaternion<double>::
//						FromTwoVectors(origin_normal, target_normal).toRotationMatrix();
//		Eigen::Matrix4d rotate_matrix_44 = Eigen::Matrix4d::Identity();
//		rotate_matrix_44.block<3, 3>(0, 0) = rotate_matrix_33;
//
//		auto final_rotate_matrix = re_to_zero * rotate_matrix_44 * to_zero;
//		
//		Eigen::Vector4d homogeneoue_coord(0., 0., 0., 1.);
//		std::array<Eigen::Vector3d, 3> v_rotated_array;
//		std::array<Eigen::Vector3d, 3> grad_array;
//
//		for(int j = 0; j < 3; ++j)
//		{
//			homogeneoue_coord.topRows(3) = v_array[j];
//			v_rotated_array[j] = (final_rotate_matrix * homogeneoue_coord).topRows(3);
//		}
//
//		double tri_area;
//		double len_a = (v_rotated_array[0] - v_rotated_array[1]).norm();
//		double len_b = (v_rotated_array[1] - v_rotated_array[2]).norm();
//		double len_c = (v_rotated_array[2] - v_rotated_array[0]).norm();
//		double p = (len_a + len_b + len_c) / 2.0;
//		tri_area = std::sqrt(p * (p - len_a) * (p - len_b) * (p - len_c));
//		/*CommonAlg::get_tri_area(v_rotated_array[0], 
//														v_rotated_array[1], 
//														v_rotated_array[2], area);*/
//
//
//		Eigen::Vector3d vecij = v_rotated_array[1] - v_rotated_array[0];
//		Eigen::Vector3d vecjk = v_rotated_array[2] - v_rotated_array[1];
//		Eigen::Vector3d vecki = v_rotated_array[0] - v_rotated_array[2];
//		Eigen::Vector3d facet_normal = vecij.cross(-vecki);
//		tri_area *= 2.0;
//
//		// (X_k - X_j)^90° / 2A;
//		grad_array[0] = facet_normal.cross(vecjk).normalized() * vecjk.norm() / tri_area;
//
//		// (X_i - X_k)^90° / 2A;
//		grad_array[1] = facet_normal.cross(vecki).normalized() * vecki.norm() / tri_area;
//
//		// (X_j - X_i)^90° / 2A;
//		grad_array[2] = facet_normal.cross(vecij).normalized() * vecij.norm() / tri_area;
//
//		//CommonAlg::get_tri_grad(v_rotated_array, grad_array);
//
//		Eigen::Vector3d wT = grad_array[0].cwiseProduct(v_rotated_array[0]) +
//							 grad_array[1].cwiseProduct(v_rotated_array[1]) +
//							 grad_array[2].cwiseProduct(v_rotated_array[2]);
//
//		/////////////////////////////////////////////////////////////
//		// output:
//		// 6.18122 12.9093 16.9896
//		// -14.9932 -15.913 -0.63831
//		// 8.81193 3.00369 -16.3513
//		// 6.18122 -14.9932 8.81193 
//		// 12.9093 -15.913 3.00369 
//		// 16.9896 -0.63831 -16.3513	
//		// for(int i = 0; i < grad_array.size(); ++i)
//		// {
//		// 	std::cout << grad_array[i].x() << " " << 
//		// 	grad_array[i].y() << " " << grad_array[i].z() << "\n";
//		// }
//		// igl::grad(A, T, G);
//		// Eigen::MatrixXd G_dence = Eigen::MatrixXd(G);
//
//		// for(int i = 0; i < G_dence.rows(); ++i)
//		// {
//		// 	for(int j = 0; j < G_dence.cols(); ++j)
//		// 	{
//		// 		std::cout << G_dence(i, j) << " ";
//		// 	}
//		// 	std::cout << "\n";
//		// }
//		/////////////////////////////////////////////////////////////
//
//
//		Eigen::Vector3d divergence;
//
//		for(int j = 0; j < 3; ++j)
//		{
//			int vid = vidxs[j];
//			if(vid < bnd_num) continue;
//			divergence = grad_array[j].cwiseProduct(wT);
//			divergence *= tri_area;
//			B.row(vid) -= divergence;
//		}
//
//	}
//
//	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//  solver.compute(A);
//	
//  Eigen::VectorXd res_x = solver.solve(B.col(0));
//  Eigen::VectorXd res_y = solver.solve(B.col(1));
//  Eigen::VectorXd res_z = solver.solve(B.col(2));
//
//	// std::cout << res_x.rows() << " " << res_x.cols() << "\n";
//	// std::cout << hole_mesh_.vertex.rows() << " " << hole_mesh_.vertex .cols() << "\n";
//
//	hole_mesh_.Vertex.col(0) = res_x;
//	hole_mesh_.Vertex.col(1) = res_y;
//	hole_mesh_.Vertex.col(2) = res_z;
//
//
//	// /// for test
//	// Mesh dmesh;
//	// dmesh.vertex.resize(r_v_vec.size(), 3);
//	// for(int i = 0; i < r_v_vec.size(); ++i)
//	// {
//	// 	dmesh.vertex.row(i) = r_v_vec[i];
//	// }
//	// dmesh.topo.resize(r_t_vec.size(), 3);
//	// for(int i = 0; i < r_t_vec.size(); ++i)
//	// {
//	// 	dmesh.topo.row(i) = r_t_vec[i];
//	// }
//	// writeVTK("./debug_rotate.vtk", dmesh);
//	// /// end test
//
//
//
//}



void SurfaceHoleFilling::fair(const int k)
{
	//std::cout << "fairing\n";
	// Eigen::MatrixXd bc, Z;
	// Eigen::VectorXi C(3);
	// C << 0, 1, 2;
	// Eigen::VectorXi R;
	// R.resizeLike(this->bnd_loop_);
	// R.setLinSpaced(this->bnd_loop_.rows(), 0, this->bnd_loop_.rows() - 1);
	// igl::slice(this->hole_mesh_.vertex, R, C, bc);
	// igl::harmonic(this->hole_mesh_.vertex, this->hole_mesh_.topo, R, bc, k, Z);
	// this->hole_mesh_.vertex = Z;
	
	// Eigen::MatrixXi hole_facet;
	// hole_facet.resizeLike(this->hole_mesh_.topo);
	// for(int i = 0; i < hole_facet.rows(); ++i)
	// {
	// 	for(int j = 0; j < 3; ++j)
	// 	{
	// 		int curF = this->hole_mesh_.topo(i, j);
	// 		if(curF < this->bnd_loop_.size())
	// 		{
	// 			hole_facet(i, j) = this->bnd_loop_(curF);
	// 		}
	// 		else
	// 		{
	// 			hole_facet(i, j) = curF - this->bnd_loop_.size() + this->xyz_mesh_.vertex.rows();
	// 		}
	// 	}
	// }
	// Eigen::MatrixXd hole_vertex, bc, Z;
	// Eigen::VectorXi R;
	// int R_n = this->hole_mesh_.vertex.rows() - this->bnd_loop_.size();
	// R.resize(R_n);
	// R.setLinSpaced(R_n, this->bnd_loop_.size(), this->hole_mesh_.vertex.rows() - 1);
	// igl::slice(this->hole_mesh_.vertex, R, vec3i(0, 1, 2), hole_vertex);

	// igl::cat(1, this->xyz_mesh_.vertex, hole_vertex, this->all_mesh_.vertex);
	// igl::cat(1, this->xyz_mesh_.topo, hole_facet, this->all_mesh_.topo);

	// igl::slice(this->all_mesh_.vertex, this->bnd_loop_, vec3i(0, 1, 2), bc);

	// igl::harmonic(this->all_mesh_.vertex, this->all_mesh_.topo, this->bnd_loop_, bc, k, Z);
	
	// this->all_mesh_.vertex = Z;

	Eigen::MatrixXi hole_facet;
	hole_facet.resizeLike(this->hole_mesh_.Topo);
	for(int i = 0; i < hole_facet.rows(); ++i)
	{
		for(int j = 0; j < 3; ++j)
		{
			int curF = this->hole_mesh_.Topo(i, j);
			if(curF < this->bnd_loop_.size())
			{
				hole_facet(i, j) = this->bnd_loop_(curF);
			}
			else
			{
				hole_facet(i, j) = curF - this->bnd_loop_.size() + this->xyz_mesh_.Vertex.rows();
			}
		}
	}
	Eigen::MatrixXd hole_vertex, bc, Z;
	Eigen::VectorXi R;
	int R_n = this->hole_mesh_.Vertex.rows() - this->bnd_loop_.size();
	R.resize(R_n);
	R.setLinSpaced(R_n, this->bnd_loop_.size(), this->hole_mesh_.Vertex.rows() - 1);
	//igl::slice(this->hole_mesh_.vertex, R, vec3i(0, 1, 2), hole_vertex);
	my_slice(this->hole_mesh_.Vertex, R, Eigen::Vector3i(0, 1, 2), hole_vertex);
	//igl::cat(1, this->xyz_mesh_.vertex, hole_vertex, this->all_mesh_.vertex);
	//igl::cat(1, this->xyz_mesh_.topo, hole_facet, this->all_mesh_.topo);

	my_cat(this->xyz_mesh_.Vertex, hole_vertex, this->all_mesh_.Vertex);
	my_cat(this->xyz_mesh_.Topo, hole_facet, this->all_mesh_.Topo);
	R.resize(this->xyz_mesh_.Vertex.rows());
	R.setLinSpaced(this->xyz_mesh_.Vertex.rows(), 0, this->xyz_mesh_.Vertex.rows() - 1);

	//igl::slice(this->all_mesh_.vertex, R, vec3i(0, 1, 2), bc);
	 my_slice(this->all_mesh_.Vertex, R, Eigen::Vector3i(0, 1, 2), bc);
	 my_harmonic(this->all_mesh_.Vertex, this->all_mesh_.Topo, R, bc, k, Z);
	//igl::harmonic(this->all_mesh_.vertex, this->all_mesh_.topo, R, bc, k, Z);
	
	this->all_mesh_.Vertex = Z;

	R.resizeLike(this->bnd_loop_);
	R.setLinSpaced(this->bnd_loop_.rows(), 0, this->bnd_loop_.rows() - 1);
	
	for(int i = 0; i < this->bnd_loop_.size(); i++)
	{
		this->hole_mesh_.Vertex.row(i) = this->all_mesh_.Vertex.row(this->bnd_loop_[i]);
	}

	for(int i = this->xyz_mesh_.Vertex.rows(), 
					j = this->bnd_loop_.size();
					
					i < this->all_mesh_.Vertex.rows(); 
					
					i++, j++)
	{
		this->hole_mesh_.Vertex.row(j) = this->all_mesh_.Vertex.row(i);
	}
	//MESHIO::remesh(this->all_mesh_);
}