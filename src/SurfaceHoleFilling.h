#pragma once
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "meshIO.h"
#include "Eigen/Geometry"


class SurfaceHoleFilling
{

public:
	SurfaceHoleFilling(Mesh& xyz_mesh);

	template <typename Derived, class MatC>
	void my_cat(
		const Eigen::MatrixBase<Derived>& A,
		const Eigen::MatrixBase<Derived>& B,
		MatC& C);

	template <
		typename DerivedX,
		typename DerivedR,
		typename DerivedC,
		typename DerivedY>
	void my_slice(
		const Eigen::DenseBase<DerivedX>& X,
		const Eigen::DenseBase<DerivedR>& R,//行索引矩阵
		const Eigen::DenseBase<DerivedC>& C,//列索引矩阵
		Eigen::PlainObjectBase<DerivedY>& Y);


	void SurfaceHoleFilling::my_triangle_triangle_adjacency(
		const Eigen::MatrixXi& F,
		Eigen::MatrixXi& TT,
		Eigen::MatrixXi& TTi);

	void SurfaceHoleFilling::my_boundary_loop(
		const Eigen::MatrixXi& F,
		Eigen::VectorXi& L);


	void my_cotmatrix(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		Eigen::SparseMatrix<double>& L);

	bool my_harmonic(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& b,
		const Eigen::MatrixXd& bc,
		const int k,
		Eigen::MatrixXd& W);

	// private:

	  /**
	   * @brief 洞的初始网格生成
	   *
	  */
	void init_hole_fill();


	//void cal_default_guide_h();

	/**
	 * @brief 通过possion方程来求解网格变形
	 *
	 * @param [in] guide_h 变形的指导梯度场
	*/
	//void possion_deformation(const Eigen::MatrixXd& guide_h);

	/**
	 * @brief 保证网格是k阶连续的
	 *
	 * @param k
	*/
	void fair(const int k);

	const Mesh& get_xyz_mesh();

	const Mesh& get_hole_mesh();
	const Mesh& get_all_mesh();
	const Eigen::MatrixXd& get_default_guide_h();

private:

	/**
	 * @brief 指导网格变形的向量场，默认通过距离插值得到，后期可以通过人为指定
	 *
	*/
	Eigen::MatrixXd guide_h_;

	/**
	 * @brief 网格数据结构存储
	 *
	*/
	Mesh& xyz_mesh_;

	Mesh hole_mesh_;

	Mesh all_mesh_;

	Eigen::VectorXi bnd_loop_;
};






