#include <vector>
#include <Eigen/Dense>

namespace TIGER {
	int resetOrientation(std::vector<std::vector<double>>& point_list, std::vector<std::vector<int>>& facet_list, std::vector<int>& blockMark);
	int resetOrinetation(Eigen::MatrixXd& point_list, Eigen::MatrixXi& facet_list, Eigen::VectorXi& block_mark);
};
