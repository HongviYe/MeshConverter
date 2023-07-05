#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>
#include "MeshOrient.h"
using namespace std;

/**
 * @brief Seperate string origin by given a set of patterns.
 * 
 * @param origin 
 * @param patterns If it meets one of the patterns, delete the charactor and split it from this index.
 * @return std::vector<std::string> 
 */
std::vector<std::string> seperate_string(std::string origin) {
    std::vector<std::string> result;
    std::stringstream ss(origin);
    while (ss >> origin) result.push_back(origin);
    return result;
}

int read_VTK(const std::string& filename, vector<vector<double>> &V, vector<vector<int>> &T) {
    std::ifstream vtk_file;
    vtk_file.open(filename);
    if(!vtk_file.is_open()) {
        cout << "No such file." << endl;
        return 0;
    }
    std::string vtk_type_str = "POLYDATA ";
    char buffer[256];
    while(!vtk_file.eof()) {
        vtk_file.getline(buffer, 256);
        std::string line = (std::string)buffer;
        if(line.length() < 2 || buffer[0] == '#')
            continue;
        if(line.find("DATASET") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            if(words[1] == "POLYDATA") 
                vtk_type_str = "POLYGONS ";
            else if(words[1] == "UNSTRUCTURED_GRID")
                vtk_type_str = "CELLS ";
            else {
                cout << "WARN : The format of VTK file is illegal, No clear DATASET name." << endl;
            }
        }
        if(line.find("POINTS ") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            int nPoints = stoi(words[1]);
            for(int i = 0; i < nPoints; i++) {
                vtk_file.getline(buffer, 256);
                words = seperate_string(std::string(buffer));
                V.push_back( {stod(words[0]), stod(words[1]), stod(words[2])} );
            }
        }
        if(line.find(vtk_type_str) != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            int nFacets = stoi(words[1]);
            for(int i = 0; i < nFacets; i++) {
                vtk_file.getline(buffer, 256);
                words = seperate_string(std::string(buffer));
                T.push_back( {stoi(words[1]), stoi(words[2]), stoi(words[3])} );
            }
        }
    }
    vtk_file.close();
    return 1;
}

int write_VTK(const std::string& filename, const vector<vector<double>> &V, const vector<vector<int>> &T, const vector<int>& block_mark = vector<int>()) {
    std::ofstream f(filename);
    if(!f.is_open()) {
        cout << "Write VTK file failed. " << endl; 
        return 0;
    }
    f.precision(std::numeric_limits<double>::digits10 + 1);
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "TetWild Mesh" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET UNSTRUCTURED_GRID" << std::endl;
    f << "POINTS " << V.size() << " double" << std::endl;
    for(const auto & v : V)
        f << v[0] << " " << v[1] << " " << v[2] << std::endl;
    f << "CELLS " << T.size() << " " << T.size() * 4 << std::endl;
    for(const auto & t : T) {
        f << t.size() << " ";
        for(int i : t)
            f << i << " ";
        f << std::endl;
    }
    f << "CELL_TYPES " << T.size() << std::endl;
    int cellType = 5;
    for(int i = 0; i < T.size(); i++)
        f << cellType << std::endl;
    if(block_mark.size() == T.size()) {
        f << "CELL_DATA " << T.size() << std::endl;
        f << "SCALARS cell2bodyid int 1" << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for(int mark : block_mark) {
            f << mark << std::endl;
        }
    }
    f.close();
    return 1;
}

int main(int argc, char** argv) {
    if(argc < 2) {
        cout << "No file input." << endl;
        return -1;
    }
    string input_filename = argv[1];
    string output_filename = input_filename.substr(0, input_filename.find_last_of('.')) + ".o.vtk";
    vector<vector<double>> plist;
    vector<vector<int>> flist;
    int status = read_VTK(input_filename, plist, flist);
    if (status == 0) return -1;
    cout << "VTK file read succeed." << " - " << input_filename << endl;
    cout << "  " << plist.size() << " points" << endl;
    cout << "  " << flist.size() << " facets" << endl;
    vector<int> bMark;
    status = TIGER::resetOrientation(plist, flist, bMark);
    if(status == 0) return -1;
    cout << "Reset orientation finished." << endl;
    int blockNum = 0;
    for(auto &b : bMark)
        blockNum = std::max(b, blockNum);
    cout << "  " << ++blockNum << " blocks" << endl;
    cout << "Output file will be written to " << output_filename << endl;
    status = write_VTK(output_filename, plist, flist, bMark);
    return status;
}