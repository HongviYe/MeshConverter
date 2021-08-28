#include "meshIO.h"
#include "CLI11.hpp"
#include "MeshOrient.h"
#include "fstream"

#define _DEBUG_ 1

using namespace std;

int main(int argc, char** argv) {
    CLI::App app{"MeshConveter"};
    string input_filename;
	string input_filename_ex;
	bool exportMESH = false;
	bool exportVTK = false;
	bool exportPLY = false;
	bool exportPLS = false;
	bool exportfacet = false;
	bool resetoritation = false;
	bool exportEpsVTK = false;
	bool reverseFacetOrient = false;
	bool repairVtk = false;

	vector<double> rotateVec;
	vector<double> boxVec;
	app.add_option("-b", boxVec, "input bounding box. Format is (length, width, hight)");
	app.add_option("-r", rotateVec, "input rotate param. Format is (start_x, start_y, start_z, end_x, end_y, end_z, angle) or (end_x, end_y, end_z, angle). angle value scale is (0, 2).");
    app.add_option("-i", input_filename, "input filename. (string, required)")->required();
	app.add_option("-p", input_filename_ex, "input filename. (string, required)");
	app.add_flag("-k", exportVTK, "Write mesh in VTK format.");
	app.add_flag("-e", exportEpsVTK, "Set eps in VTK format.");
	app.add_flag("-o", reverseFacetOrient, "Reverse Facet Orient.");
	app.add_flag("-m", exportMESH, "Write mesh in MESH/MEDIT format.");
	app.add_flag("-y", exportPLY, "Write mesh in PLY format.");
	app.add_flag("-s", exportPLS, "Write mesh in PLS format.");
	app.add_flag("-f", exportfacet, "Write mesh in facet format.");
	app.add_flag("--orient", resetoritation, "Regularize oritation");
	app.add_flag("--repair", repairVtk, "Repair vtk file for the area is equal to zero.");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

	size_t input_dotpos = input_filename.find_last_of('.');
	string input_postfix = input_filename.substr(input_dotpos + 1, input_filename.length() - input_dotpos - 1);


	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXi M;
	int cou;
	std::map<int, double> mpd;
	std::map<int, vector<int>> mpi;

    if(input_postfix == "vtk")
        MESHIO::readVTK(input_filename, V, F, M);
    else if(input_postfix == "mesh")
        MESHIO::readMESH(input_filename, V, F, M);
	else if (input_postfix == "pls")
		MESHIO::readPLS(input_filename, V, F, M);
    else {
        cout << "Unsupported input format - " << input_postfix << endl;
        return -1;
    }

    if(exportEpsVTK)
    	MESHIO::readEPS(input_filename_ex, cou, mpd, mpi);


	//********* Rotate *********
	if(!rotateVec.empty())
	MESHIO::rotatePoint(rotateVec, V, F);
	//******** Rotated *********

	//********* Add Box *********
	if(!boxVec.empty())
		MESHIO::addBox(boxVec, V, F, M);
	//********* Add Box *********

	//********* Regularize mesh oritation *********
	if(resetoritation)
	MESHIO::resetOrientation(V, F, M);
	//********* Regularize mesh oritation *********

	//********* modify facet orient ******
	if(reverseFacetOrient){
		MESHIO::reverseOrient(F);
	}
	//********* modify facet orient ******

	//********** repair ********

	if(repairVtk){
		MESHIO::repair(V, F, M);
	}


	//********** repair ********





	
    if(exportVTK) {
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.vtk";
        MESHIO::writeVTK(output_filename, V, F, M);
    }
    if(exportMESH) {
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.mesh";
        MESHIO::writeMESH(output_filename, V, F);
    }
    if(exportPLY) { 
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.ply";
        MESHIO::writePLY(output_filename, V, F);
    }
    if(exportPLS){
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.pls";
        MESHIO::writePLS(output_filename, V, F, M);
    }
	if (exportfacet) {
		string output_filename = input_filename.substr(0, input_dotpos) + ".o.facet";
		MESHIO::writeFacet(output_filename, V, F, M);
	}
	if(exportEpsVTK){
		string output_filename = input_filename.substr(0, input_dotpos) + ".eps.vtk";
		MESHIO::writeEpsVTK(output_filename, V, F, cou, mpd, mpi);
	}
    return 0;
}