#include <Eigen/Dense>

// Futher: Attribute should contain these data structure: double,int, vector(double,3), tensor(double,9), and no quantity constrain should be in there
struct Mesh {
	Eigen::MatrixXd Vertex;
	Eigen::MatrixXi Topo;
	Eigen::MatrixXi Masks;
};
