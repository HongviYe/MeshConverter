#ifndef MESHALGORITHM_H
#define MESHALGORITHM_H

#include <Eigen/Dense>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include "meshIO.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MESHIO{

	bool rotatePoint(std::vector<double> rotateVec, Mesh &mesh);
	bool addBox(std::vector<double> boxVec, Mesh &mesh);
	bool reverseOrient(Eigen::MatrixXi &T);
	bool repair(Mesh &mesh);
	bool resetOrientation(Mesh &mesh,bool reset_mask=false);
	bool removeDulplicatePoint(Eigen::MatrixXd& V, Eigen::MatrixXi& T, double eps);
	bool createBox(std::vector<double> create_box, Mesh &mesh);
	bool buildTuttleParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv);
	bool buildHarmonicParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv);
	bool shuffleSurfaceid(int num, Mesh& mesh);


  // called by the algorithm above.
	void dfs_get_loop2(
  int cur, int pre, 
  std::vector<bool>& vis, 
  std::vector<std::vector<int>>& G, 
  std::vector<int>& path, 
  std::vector<std::vector<int>>& loop_lst);
	void boundary_loop_by_dfs2(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& bnd);


}

#endif
 