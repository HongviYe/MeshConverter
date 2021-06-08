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
    app.add_option("-i", input_filename, "input filename. (string, required)")->required();
    app.add_flag("-k", exportVTK, "Write mesh in VTK format.");
    app.add_flag("-m", exportMESH, "Write mesh in MESH/MEDIT format.");
    app.add_flag("-y", exportPLY, "Write mesh in PLY format.");
    
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
    else {
        cout << "Unsupported input format - " << input_postfix << endl;
        return -1;
    }

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

    return 0;
}