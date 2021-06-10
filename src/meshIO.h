#ifndef MESHIO_H
#define MESHIO_H

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace MESHIO{
int readVTK(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M, std::string mark_pattern = "");
int writeVTK(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::MatrixXi M = Eigen::MatrixXi(), std::string mark_pattern = "");
int readMESH(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T);
int writeMESH(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
int readPLY(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T);
int writePLY(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);

}

#endif