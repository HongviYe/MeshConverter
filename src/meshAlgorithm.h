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

}

#endif
 