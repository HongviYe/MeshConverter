#include "meshIO.h"
#include "CLI11.hpp"
#include "meshAlgorithm.h"

#define _DEBUG_ 1

using namespace std;

int main(int argc, char** argv) {
    CLI::App app{"MeshConveter"};
    vector<string> input_filenames;
	string input_filename_ex;
	bool exportMESH = false;
	bool exportVTK = false;
	bool exportEpsVTK = false;
	bool exportPLY = false;
	bool exportPLS = false;
	bool exportFacet = false;
	bool exportOBJ = false;
	bool resetOritation = false;
	bool reverseFacetOrient = false;
	bool meshRepair = false;

	vector<double> rotateVec;
	vector<double> boxVec;
	app.add_option("-b", boxVec, "input bounding box. Format is (length, width, hight)");
	app.add_option("-r", rotateVec, "input rotate param. Format is (start_x, start_y, start_z, end_x, end_y, end_z, angle) or (end_x, end_y, end_z, angle). angle value scale is (0, 2).");
    app.add_option("-i", input_filenames, "input filename. (string, required, supported format: vtk, mesh, pls, obj)")->required()->expected(1, 3);
	app.add_option("-p", input_filename_ex, "input filename. (string, required)");
	app.add_flag("-k", exportVTK, "Write mesh in VTK format.");
	app.add_flag("-e", exportEpsVTK, "Set eps in VTK format.");
	app.add_flag("-m", exportMESH, "Write mesh in MESH/MEDIT format.");
	app.add_flag("-y", exportPLY, "Write mesh in PLY format.");
	app.add_flag("-s", exportPLS, "Write mesh in PLS format.");
	app.add_flag("-f", exportFacet, "Write mesh in facet format.");
	app.add_flag("-o", exportOBJ, "Write mesh in OBJ format.");
	app.add_flag("--reverse-orient", reverseFacetOrient, "Reverse Facet Orient.");
	app.add_flag("--reset-orient", resetOritation, "Regularize oritation");
	app.add_flag("--repair", meshRepair, "Repair vtk file for the area is equal to zero.");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

	string input_filename = input_filenames[0];
	size_t input_dotpos = input_filename.find_last_of('.');
	string input_postfix = input_filename.substr(input_dotpos + 1, input_filename.length() - input_dotpos - 1);

	Mesh mesh;
	int cou;
	std::map<int, double> mpd;
	std::map<int, vector<int>> mpi;

    if(input_postfix == "vtk")
        MESHIO::readVTK(input_filename, mesh);
    else if(input_postfix == "mesh")
        MESHIO::readMESH(input_filename, mesh);
	else if (input_postfix == "pls")
		MESHIO::readPLS(input_filename, mesh);
	else if (input_postfix == "ply")
		MESHIO::readPLY(input_filename, mesh);
	else if (input_postfix == "obj")
		MESHIO::readOBJ(input_filename, mesh);
	else if (input_postfix == "facet")
		MESHIO::readFacet(input_filename, mesh);
	else if (input_postfix == "node" || input_postfix == "ele" || input_postfix == "face") {
		string nodefilename = "";
		string facefilename = "";
		string elemfilename = "";
		for(auto t_filename : input_filenames) {
			size_t t_dotpos = t_filename.find_last_of('.');
			string t_postfix = t_filename.substr(t_dotpos + 1, t_filename.length() - t_dotpos - 1);
			if(t_postfix == "node")
				nodefilename = t_filename;
			else if(t_postfix == "ele")
				elemfilename = t_filename;
			else if(t_postfix == "face")
				facefilename = t_filename;
		}
		if(nodefilename.empty()) {
			cout << "Node file is required." << endl;
		}
		if(!elemfilename.empty()) {
			MESHIO::readTetgen(nodefilename, elemfilename, mesh);
		}
		else if(!facefilename.empty()) {
			MESHIO::readTetgen(nodefilename, facefilename, mesh);
		}
	}
    else {
        cout << "Unsupported input format - " << input_postfix << endl;
        return -1;
    }

	// cout for checking reading.
	cout << "Read " << mesh.Vertex.rows() << " points." << endl;
	cout << "Read " << mesh.Topo.rows() << " elements." << endl;

    if(exportEpsVTK)
    	MESHIO::readEPS(input_filename_ex, cou, mpd, mpi);

	//********* Rotate *********
	if(!rotateVec.empty())
	MESHIO::rotatePoint(rotateVec, mesh);

	//********* Add Box *********
	if(!boxVec.empty())
		MESHIO::addBox(boxVec, mesh);

	//********* Regularize mesh oritation *********
	if(resetOritation)
		MESHIO::resetOrientation(mesh);

	//********* modify facet orient ******
	if(reverseFacetOrient){
		MESHIO::reverseOrient(mesh.Topo);
	}

	//********* repair ********
	if(meshRepair){
		MESHIO::repair(mesh);
	}

	//********* export ********
    if(exportVTK) {
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.vtk";
        MESHIO::writeVTK(output_filename, mesh, "part");
    }
    if(exportMESH) {
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.mesh";
        MESHIO::writeMESH(output_filename, mesh);
    }
    if(exportPLY) { 
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.ply";
        MESHIO::writePLY(output_filename, mesh);
    }
    if(exportPLS){
        string output_filename = input_filename.substr(0, input_dotpos) + ".o.pls";
        MESHIO::writePLS(output_filename, mesh);
    }
	if(exportFacet) {
		string output_filename = input_filename.substr(0, input_dotpos) + ".o.facet";
		MESHIO::writeFacet(output_filename, mesh);
	}
	if(exportEpsVTK){// This is to generate AutoGrid to control local eps.
		string output_filename = input_filename.substr(0, input_dotpos) + ".eps.vtk";
		MESHIO::writeEpsVTK(output_filename, mesh, cou, mpd, mpi);
	}
	if(exportOBJ) {
		string output_filename = input_filename.substr(0, input_dotpos) + ".o.obj";
		MESHIO::writeOBJ(output_filename, mesh);
	}

    return 0;
}
