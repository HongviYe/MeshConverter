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

// Futher: Attribute should contain these data structure: double,int, vector(double,3), tensor(double,9), and no quantity constrain should be in there
struct Mesh {
	Eigen::MatrixXd Vertex;
	Eigen::MatrixXi Topo;
	Eigen::MatrixXi Masks;
};

namespace MESHIO{

int readVTK(std::string filename, Mesh &mesh, std::string mark_pattern = "");
int readEPS(std::string filename, int& cou, std::map<int, double>& mpd, std::map<int, std::vector<int>>& mpi);
int readMESH(std::string filename, Mesh &mesh);
int readPLS(std::string filename, Mesh &mesh);
int readPLY(std::string filename, Mesh &mesh); //TODO
int readOBJ(std::string filename, Mesh &mesh);
int readFacet(std::string filename, Mesh &mesh);
int readTetgen(std::string nodefilename, std::string elefilename, Mesh &mesh);

int writeVTK(std::string filename, const Mesh &mesh, std::string mark_pattern = "");
int writeEpsVTK(std::string filename, const Mesh &mesh, int& cou, std::map<int, double>& mpd, std::map<int, std::vector<int>>& mpi, std::string mark_pattern = "");
int writeMESH(std::string filename, const Mesh &mesh);
int writePLY(std::string filename, const Mesh &mesh);
int writePLS(std::string filename, const Mesh &mesh);
int writeFacet(std::string filename, const Mesh &mesh);
int writeOBJ(std::string filename, const Mesh &mesh);

}

#endif
 