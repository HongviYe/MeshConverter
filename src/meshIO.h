#ifndef MESHIO_H
#define MESHIO_H

#include <Eigen/Dense>
#include <map>
#include <string>
#include <iostream>
#include <vector>
#define DEBUG_INFO() printf("File:%s, Line:%d, Function:%s\n", __FILE__, __LINE__ , __FUNCTION__)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace std;
namespace MESHIO{
int readVTK(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M, std::string mark_pattern = "");
int readEPS(std::string filename, int& cou, std::map<int, double>& mpd, std::map<int, vector<int>>& mpi);
int writeVTK(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::MatrixXi M = Eigen::MatrixXi(), std::string mark_pattern = "");
int writeEpsVTK(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, int& cou, std::map<int, double>& mpd, std::map<int, vector<int>>& mpi, std::string mark_pattern = "");
int readMESH(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M);
int writeMESH(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
int readPLY(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T);
int writePLY(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
int writePLS(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::MatrixXi M);
int readPLS(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M);
int writeFacet(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::MatrixXi M );
bool rotatePoint(vector<double> rotateVec, Eigen::MatrixXd &V, Eigen::MatrixXi &T);
bool addBox(vector<double> boxVec, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M);
bool reverseOrient(Eigen::MatrixXi &T);
bool repair(Eigen::MatrixXd &V, Eigen::MatrixXi &T,Eigen::MatrixXi M = Eigen::MatrixXi());
}

#endif
 