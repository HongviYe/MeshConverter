#ifndef MESHIO_H
#define MESHIO_H

#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

#define DEBUG_INFO() printf("File:%s, Line:%d, Function:%s\n", __FILE__, __LINE__ , __FUNCTION__)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MESHIO{

int readVTK(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M, std::string mark_pattern = "");
int readEPS(std::string filename, int& cou, std::map<int, double>& mpd, std::map<int, std::vector<int>>& mpi);
int readMESH(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M);
int readPLS(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M);
int readPLY(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M); //TODO
int readOBJ(std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M);

int writeVTK(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::MatrixXi M = Eigen::MatrixXi(), std::string mark_pattern = "");
int writeEpsVTK(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, int& cou, std::map<int, double>& mpd, std::map<int, std::vector<int>>& mpi, std::string mark_pattern = "");
int writeMESH(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
int writePLY(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
int writePLS(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::MatrixXi &M);
int writeFacet(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::MatrixXi &M);
int writeOBJ(std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXi& T, const Eigen::MatrixXi &M);

bool rotatePoint(std::vector<double> rotateVec, Eigen::MatrixXd &V, Eigen::MatrixXi &T);
bool addBox(std::vector<double> boxVec, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXi &M);
bool reverseOrient(Eigen::MatrixXi &T);
bool repair(Eigen::MatrixXd &V, Eigen::MatrixXi &T,Eigen::MatrixXi M = Eigen::MatrixXi());

}

#endif
 