#include "meshIO.h"
#include "CLI11.hpp"
#include "fstream"

using namespace std;

int main(int argc, char** argv) {
    CLI::App app{"VTK to PLY"};
    string input_filename;    
    bool exportMESH = false;
    bool exportVTK = false;
    bool exportPLY = false;
    bool exportPLS = false;
	bool exportfacet = false;
	vector<double> rotateVec;
	app.add_option("-r", rotateVec, "input rotate param. Format is (start_x, start_y, start_z, end_x, end_y, end_z, angle) or (end_x, end_y, end_z, angle). angle value scale is (0, 2).");
    app.add_option("-i", input_filename, "input filename. (string, required)")->required();
    app.add_flag("-k", exportVTK, "Write mesh in VTK format.");
    app.add_flag("-m", exportMESH, "Write mesh in MESH/MEDIT format.");
    app.add_flag("-y", exportPLY, "Write mesh in PLY format.");
    app.add_flag("-s", exportPLS, "Write mesh in PLS format.");
	app.add_flag("-f", exportfacet, "Write mesh in facet format.");
    
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

    if(input_postfix == "vtk")
        MESHIO::readVTK(input_filename, V, F, M);
    else if(input_postfix == "mesh")
        MESHIO::readMESH(input_filename, V, F);
	else if (input_postfix == "pls")
		MESHIO::readPLS(input_filename, V, F, M);
    else {
        cout << "Unsupported input format - " << input_postfix << endl;
        return -1;
    }

	//********* Rotate *********
	double start_x = 0.0, start_y = 0.0, start_z = 0.0;
	double end_x = 0.0, end_y = 0.0, end_z = 0.0;
	double angle = 0.0;
	if(rotateVec.size() > 0){
		if(rotateVec.size() == 4) {
			end_x = rotateVec[0];
			end_y = rotateVec[1];
			end_z = rotateVec[2];
			angle = rotateVec[3];
		}
		else if(rotateVec.size() == 7){
			start_x = rotateVec[0]; start_y = rotateVec[1]; start_z = rotateVec[2];
			end_x = rotateVec[3]; end_y = rotateVec[4]; end_z = rotateVec[5];
			angle = rotateVec[6];
		}else{
			std::cout << "The format is Error.Format is (start_x, start_y, start_z, end_x, end_y, end_z, angle) or (end_x, end_y, end_z, angle). angle value scale is (0, 2).";
			return -1;
		}
	}

	if(rotateVec.size() == 7) {
		std::cout << "Rotating\n";
		end_x -= start_x;
		end_y -= start_y;
		end_z -= start_z;
		for (int i = 0; i < V.rows(); i++) {
			V(i, 0) -= start_x;
			V(i, 1) -= start_y;
			V(i, 2) -= start_z;
		}
	}
	if(rotateVec.size() > 0) {
		Eigen::AngleAxisd rotationVector(M_PI * angle, Eigen::Vector3d(end_x, end_y, end_z));
		Eigen::Matrix3d rotationMatrix = Eigen::Matrix3d::Identity();
		rotationMatrix = rotationVector.toRotationMatrix();

		std::cout << rotationMatrix << '\n';

		Eigen::MatrixXd tmp = rotationMatrix * V.transpose();
		V = tmp.transpose();
	}
	if(rotateVec.size() == 7){
		for(int i = 0; i < V.rows(); i++){
			V(i, 0) += start_x;
			V(i, 1) += start_y;
			V(i, 2) += start_z;
		}
		std::cout << "Rotated\n";
	}
	//******** Rotated ********

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
        MESHIO::writePLS(output_filename, V, F);
    }
	if (exportfacet) {
		string output_filename = input_filename.substr(0, input_dotpos) + ".o.facet";
		MESHIO::writeFacet(output_filename, V, F, M);
	}
    return 0;
}