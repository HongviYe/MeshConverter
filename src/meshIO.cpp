#include "meshIO.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <time.h>

#define BUFFER_LENGTH 256

using namespace std;

/**
 * @brief Seperate string origin by given a set of patterns.
 * 
 * @param origin 
 * @return std::vector<std::string> 
 */
std::vector<std::string> seperate_string(std::string origin) {
    std::vector<std::string> result;
    stringstream ss(origin);
    while(ss >> origin) result.push_back(origin);
    return result;
}

// The eps file format is 
/* eps 0.005          
 * id 1 5 7 9 100
 * eps 0.008
 * id 56 78 79
 * ...
 */
int MESHIO::readEPS(std::string filename, int& cou, std::map<int, double>& mpd, std::map<int, vector<int>>& mpi) {
	std::ifstream eps_file;
	eps_file.open(filename);
	if(!eps_file.is_open()){
		std::cout << "No Such file. - " << filename << std::endl;
		return -1;
	}
	string line;
	string _;
	stringstream input;
    // int cou = 0;
    // std::map<int, double> mpd;
    // std::map<int, vector<int>> mpi;
	while(getline(eps_file, line)){
		input.clear();
		input.str(line);
		input >> _;
		if(_[0] == '#') continue;
		else if(_ == "eps" || _ == "EPS"){
            cou++;
			input >> mpd[cou];
		}
		else if(_ == "id" || _ == "ID"){
			int tmpId;
			while (input >> tmpId)
			{
				mpi[cou].push_back(tmpId);
			}
		}
		else{
			continue;
		}
	}
	return 0;
}

int MESHIO::readVTK(std::string filename, Mesh& mesh, std::string mark_pattern) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    M.resize(1, 1);
    int nPoints = 0;
    int nFacets = 0;
    std::ifstream vtk_file;
    vtk_file.open(filename);
    if(!vtk_file.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    std::string vtk_type_str = "POLYDATA ";
    char buffer[BUFFER_LENGTH];
    while(!vtk_file.eof()) {
        vtk_file.getline(buffer, BUFFER_LENGTH);
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
                std::cout << "The format of VTK file is illegal, No clear DATASET name. - " << filename << std::endl;
            }
        }
        if(line.find("POINTS ") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nPoints = stoi(words[1]);
            V.resize(nPoints, 3);
            for(int i = 0; i < nPoints; i++) {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                V.row(i) << stod(words[0]), stod(words[1]), stod(words[2]);
            }
        }
        if(line.find(vtk_type_str) != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nFacets = stoi(words[1]);
            T.resize(nFacets, stoi(words[2]) / nFacets - 1);
            for(int i = 0; i < nFacets; i++) {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                for(int j = 0; j < stoi(words[0]); j++) 
                    T(i, j) = stoi(words[j + 1]);
            }
        }
        if(line.find("CELL_DATA ") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            if(stoi(words[1]) != nFacets) {
                std::cout << "The number of CELL_DATA is not equal to number of cells. -" << filename;
                std::cout << "Ignore CELL_DATA" << std::endl;
                return 0;
            }
            vtk_file.getline(buffer, BUFFER_LENGTH);
            std::string data_type = seperate_string(std::string(buffer))[1];
            if(data_type != mark_pattern) 
                continue;
            M.resize(nFacets, 1);
			for (int i = 0; i < nFacets; i++)
				M(i, 0) = 0;
            vtk_file.getline(buffer, BUFFER_LENGTH);
            for(int i = 0; i < nFacets; i++) {
                vtk_file.getline(buffer, BUFFER_LENGTH);
				int surface_id=stoi(std::string(buffer));
				M.row(i) << surface_id;
					
            }
        }
    }
    vtk_file.close();
    return 1;
}

int MESHIO::readOBJ(std::string filename, Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    ifstream objFile;
    objFile.open(filename);
    if(!objFile.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    char buffer[BUFFER_LENGTH];
    vector<vector<double>> plist;
    vector<vector<int>> flist;
    vector<int> mlist;
    int curMark = 0;
    while(!objFile.eof()) {
        objFile.getline(buffer, BUFFER_LENGTH);
        if(buffer[0] == 'v') {
            vector<string> words = seperate_string(string(buffer));
            if(words.size() < 4) {
                continue;
            }
            vector<double> coord;
            for(int ii = 1; ii < 4; ii++) {
                int wordLen = 0;
                while(wordLen < words[ii].length()) {
                    if(words[ii][wordLen] == '/') {
                        break;
                    }
                }
                coord.push_back(stod(words[ii].substr(0, wordLen)));
            }
            plist.push_back(coord);
        }
        if(buffer[0] == 'f') {
            vector<string> words = seperate_string(string(buffer));
            if(words.size() < 4) {
                continue;
            }
            vector<int> facet;
            for(int ii = 1; ii < 4; ii++) {
                int wordLen = 0;
                while(wordLen < words[ii].length()) {
                    if(words[ii][wordLen] == '/') {
                        break;
                    }
                }
                facet.push_back(stoi(words[ii].substr(0, wordLen)));
            }
            flist.push_back(facet);
            mlist.push_back(curMark);
        }
        if(buffer[0] == 'g') {
            ++curMark;
        }
    }

    V.resize(plist.size(), 3);
    for(int i = 0; i < plist.size(); i++) {
        for(int j = 0; j < 3; j++) {
            V(i, j) = plist[i][j];
        }
    }
    T.resize(flist.size(), 3);
    for(int i = 0; i < flist.size(); i++) {
        for(int j = 0; j < 3; j++) {
            T(i, j) = flist[i][j];
        }
    }
    M.resize(mlist.size(), 1);
    for(int i = 0; i < mlist.size(); i++) {
        M(i, 0) = mlist[i];
    }

    return 1;
}

int MESHIO::writeVTK(std::string filename, const Mesh &mesh, std::string mark_pattern) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    std::ofstream f(filename);
    if(!f.is_open()) {
        std::cout << "Write VTK file failed. - " << filename << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    f.precision(std::numeric_limits<double>::digits10 + 1);
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "TetWild Mesh" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET UNSTRUCTURED_GRID" << std::endl;
    f << "POINTS " << V.rows() << " double" << std::endl;
    for(int i = 0; i < V.rows(); i++)
        f << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    f << "CELLS " << T.rows() << " " << T.rows() * (T.cols() + 1) << std::endl;
    for(int i = 0; i < T.rows(); i++) {
        f << T.cols() << " ";
        for(int j = 0; j < T.cols(); j++)
            f <<  T(i, j) << " ";
        f << std::endl;
    }
    f << "CELL_TYPES " << T.rows() << std::endl;
    int cellType = 0;
    if(T.cols() == 2)
        cellType = 3;
    else if(T.cols() == 3)
        cellType = 5;
    else if(T.cols() == 4)
        cellType = 10;
    for(int i = 0; i < T.rows(); i++)
        f << cellType << std::endl;
    if(M.rows() != T.rows()) {
        f.close();
        return 1;
    }
    f << "CELL_DATA " << M.rows() << std::endl;
    f << "SCALARS " << mark_pattern << " int " << M.cols() << std::endl;
    f << "LOOKUP_TABLE default" << std::endl;
    for(int i = 0; i < M.rows(); i++) {
        for(int j = 0; j < M.cols(); j++)
            f << M(i, j);
        f << std::endl;
    }
    f << std::endl;
    f.close();
    return 1;
}

int MESHIO::writeEpsVTK(std::string filename, const Mesh &mesh, int& cou,  std::map<int, double> &mpd, std::map<int, vector<int>> &mpi, std::string mark_pattern) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	std::ofstream f(filename);
	if(!f.is_open()) {
		std::cout << "Write VTK file failed. - " << filename << std::endl;
		return -1;
	}
	std::cout << "Writing mesh to - " << filename << std::endl;
	f.precision(std::numeric_limits<double>::digits10 + 1);
	f << "# vtk DataFile Version 2.0" << std::endl;
	f << "TetWild Mesh" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;
	f << "POINTS " << V.rows() << " double" << std::endl;
	for(int i = 0; i < V.rows(); i++)
		f << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
	f << "CELLS " << V.rows() << " " << V.rows() * 2 << std::endl;
	for(int i = 0; i < V.rows(); i++) {
		f << 1 <<  " " << i ;
		f << std::endl;
	}
	f << "CELL_TYPES " << V.rows() << std::endl;
	int cellType = 1;
	for(int i = 0; i < V.rows(); i++)
		f << cellType << std::endl;

	f << "CELL_DATA " << V.rows() << std::endl;
	f << "SCALARS local_epsilon double 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;

	int curLoc = 0;

    vector<double> vec;
    vec.resize(V.rows(), -1);
    for(int i = 1; i <= cou; i++)
    {
        for(int id : mpi[i]){
            vec[id] = mpd[i];
        }
    }

	for(int i = 0; i < V.rows(); i++) {
		f << vec[i] << endl;
        if(vec[i] != -1){
            std::cout << i << " " << vec[i] << '\n';
        }
	}
	f << std::endl;
	f.close();
	return 1;
}

int MESHIO::readMESH(std::string filename, Mesh& mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    std::ifstream mesh_file;
    mesh_file.open(filename);
    if(!mesh_file.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    int dimension = 3;
    int nPoints;
    int nFacets;
    char buffer[BUFFER_LENGTH];
    while(!mesh_file.eof()) {
        mesh_file.getline(buffer, BUFFER_LENGTH);
        std::string line = (std::string)buffer;
        if(line.length() < 2)
            continue;
        std::vector<std::string> words = seperate_string(line);
        if(words.empty())
            continue;
        if(line.find("Dimension") != std::string::npos) {
            dimension = stoi(words[1]);
            std::cout << "Reading mesh dimension - " << dimension << std::endl;
        }
        if(line == "Vertices") {
            words = seperate_string(line);
            if(words.size() > 1 || words[0].length() > 8)
                continue;
            line.clear();
            while(line.empty()) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                line = (std::string)buffer;
            }
            words = seperate_string(line);
            nPoints = std::stoi(words[0]);
            std::cout << "Number of points : " << nPoints << std::endl;
            V.resize(nPoints, dimension);
            int i = 0;
            while(i < nPoints) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                if(words.size() != (dimension + 1)) continue;
                for(int j = 0; j < dimension; j++)
                    V(i, j) = std::stod(words[j]);
                i++;
            }
        }
        if(line.find("Triangles") != std::string::npos) {
            line.clear();
            while(line.empty()) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                line = (std::string)buffer;
            }
            words = seperate_string(line);
            nFacets = stoi(words[0]);
            std::cout << "Number of facets : " << nFacets << std::endl;
            T.resize(nFacets, 3);
            M.resize(nFacets, 1);
            int i = 0;
            while(i < nFacets) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                if(words.size() < 4)
                {
                    std::cout << "Warning : The number of triangles element is not equal 4.\n";
                }
                for(int j = 0; j < 3; j++)
                    T(i, j) = std::stoi(words[j]) - 1;
                M(i, 0) = std::stoi(words[3]) - 1;
                i++;
            }
        }
    }
    mesh_file.close();
    return 1;
}

int MESHIO::writeMESH(std::string filename, const Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    std::ofstream f(filename);
    if(!f.is_open()) {
        std::cout << "Write MESH file failed. - " << filename << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    f.precision(std::numeric_limits<double>::digits10 + 1);
    f << "MeshVersionFormatted 1" << std::endl;
    f << "Dimension " << V.cols() << std::endl;
    f << "Vertices" << std::endl;
    f << V.rows() << std::endl;
    for(int i = 0; i < V.rows(); i++) {
        for(int j = 0; j < V.cols(); j++)
            f << V(i, j) << " ";
        f << i + 1 << std::endl;
    }
    if(T.cols() == 3)
        f << "Triangles" << std::endl;
    else if(T.cols() == 4)
        f << "Tetrahedra" << std::endl;
    else {
        std::cout << "Unsupported format for .mesh file." << std::endl;
        return -1;
    }
    f << T.rows() << std::endl;
    for(int i = 0; i < T.rows(); i++) {
        for(int j = 0; j < T.cols(); j++)
            f << T(i, j) + 1 << " ";
        f << i + 1 << std::endl;
    }
    f.close();
    return 1;
}

int MESHIO::writePLY(std::string filename, const Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    if(T.cols() != 3) {
        std::cout << "Unsupported format for .ply file." << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    std::ofstream plyfile;

    plyfile.open(filename);
	plyfile.precision(std::numeric_limits<double>::digits10 + 1);
    plyfile << "ply" << std::endl;
    plyfile << "format ascii 1.0" << std::endl;
    plyfile << "comment VTK generated PLY File" << std::endl;
    plyfile << "obj_info vtkPolyData points and polygons: vtk4.0" << std::endl;
    plyfile << "element vertex " << V.rows() << std::endl;
    plyfile << "property float x" << std::endl;
    plyfile << "property float y" << std::endl;
    plyfile << "property float z" << std::endl;
    plyfile << "element face " << T.rows() << std::endl;
    plyfile << "property list uchar int vertex_indices" << std::endl;
    plyfile << "end_header" << std::endl;
    for(int i = 0; i < V.rows(); i++)
        plyfile << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    for(int i = 0; i < T.rows(); i++)
        plyfile << T.cols() << " " << T(i, 0) << " " << T(i, 1) << " " << T(i, 2) << std::endl;
    plyfile.close();
    return 1;
}

int MESHIO::writePLS(std::string filename, const Mesh &mesh)
{
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    if(T.cols() != 3)
    {
        std::cout << "Unsupported format for .pls file." << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    std::ofstream plsfile;
    plsfile.open(filename);
    plsfile << T.rows() << " " << V.rows() << " " << "0 0 0 0\n";
    for(int i = 0; i < V.rows(); i++)
        plsfile << i + 1 << " " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    for(int i = 0; i < T.rows(); i++)
        plsfile << i + 1 << " " << T(i, 0) + 1 << " " << T(i, 1) + 1 << " " << T(i, 2) + 1 << " " << M(i, 0) + 1 << std::endl;
        
    plsfile.close();
    std::cout << "Finish\n";
    return 1;
}

int MESHIO::readPLS(std::string filename, Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo; 
	std::ifstream pls_file;
	pls_file.open(filename);
	if (!pls_file.is_open()) {
		std::cout << "No such file. - " << filename << std::endl;
		return -1;
	}
	int nPoints = 0;
	int nFacets = 0;

	std::string str;
	std::getline(pls_file, str);
	std::stringstream ss(str);
	ss >> nFacets >> nPoints;
	V.resize(nPoints,3);
	for (int i = 0; i < nPoints; i++) {
		int index;
		double point_coordinate[3];
		pls_file >> index >> point_coordinate[0] >> point_coordinate[1] >> point_coordinate[2];
		for(int k=0;k<3;k++)
		V(index-1,k)=point_coordinate[k];
	}
	M.resize(nFacets,1);
	T.resize(nFacets, 3);
	for (int i = 0; i < nFacets; i++) {
		int index;
		int tri[3];
		pls_file >> index >> tri[0] >> tri[1] >> tri[2] >> M(i,0);
		tri[0]--; tri[1]--; tri[2]--;
		for (int k = 0; k < 3; k++)
			T(i, k) = tri[k];
	}

	return 1;
}

int MESHIO::writeFacet(std::string filename, const Mesh &mesh)
{
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	std::cout << "Writing mesh to - " << filename << std::endl;
	std::ofstream facetfile;
	facetfile.open(filename);
	facetfile.precision(std::numeric_limits<double>::digits10 + 1);
	facetfile << "FACET FILE V3.0  exported from Meshconverter http://10.12.220.71/tools/meshconverter " << std::endl;
	facetfile << 1 << std::endl;
	facetfile << "Grid" << std::endl;
	facetfile << "0, 0.00 0.00 0.00 0.00" << endl;
	facetfile << V.rows() << endl;
	for (int i = 0; i < V.rows(); i++)
		facetfile << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
	facetfile << 1 << endl;
	facetfile << "Triangles" << endl;
	facetfile << T.rows() << " 3" << endl;
	
	for (int i = 0; i < T.rows(); i++) {
		facetfile << " " << T(i, 0) + 1 << " " << T(i, 1) + 1 << " " << T(i, 2) + 1 << " 0 ";
		if (M.rows() >= T.rows())
			facetfile << M(i, 0);
		else
			facetfile << 0;
		facetfile<< " " << i + 1 << std::endl;
	}
	return 0;
}

int MESHIO::readFacet(string filename, Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo; 
	std::ifstream facetfile;
    facetfile.open(filename);
    if(!facetfile.is_open()) {
        return -1;
    }
    int nPoints = 0;
    int nFacets = 0;

    string line;
    while(!facetfile.eof()) {
        getline(facetfile, line);
        if(line[0] == '#') {
            continue;
        }
        if(line.find("Grid") != std::string::npos) {
            getline(facetfile, line);   // ignore line like "0, 0.00, 0.00, 0.00, 0.00"
            getline(facetfile, line);
            nPoints = stoi(line);
            V.resize(nPoints, 3);
            vector<string> words;
            for(int i = 0; i < nPoints; i++) {
                getline(facetfile, line);
                words = seperate_string(line);
                if(words.size() != 3) {
                    break;
                }
                for(int j = 0; j < 3; j++) {
                    V(i, j) = stod(words[j]);
                }
            }
        }
        if(line.find("Triangles") != std::string::npos) {
            getline(facetfile, line);
            vector<string> words = seperate_string(line);
            nFacets = stoi(words[0]);
            M.resize(nFacets, 1);
            T.resize(nFacets, 3);
            for(int i = 0; i < nFacets; i++) {
                getline(facetfile, line);
                words = seperate_string(line);
                for(int j = 0; j < 3; j++) {
                    T(i, j) = stoi(words[j]) - 1;
                }
                M(i, 0) = stoi(words[4]);
            }
        }
    }
    facetfile.close();
    return 1;
}

int MESHIO::writeOBJ(string filename, const Mesh &mesh) {
    // Facet group
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& F = mesh.Topo; 
	bool doGroup = (M.rows() == F.rows());
    vector<vector<int>> flist;
    for(int i = 0; i < F.rows(); i++) {
        vector<int> facet;
        facet.push_back(doGroup ? M(i, 0) : 0);
        for(int j = 0; j < F.cols(); j++) {
            facet.push_back(F(i, j));
        }
        flist.push_back(facet);
    }
    if(doGroup) {
        sort(flist.begin(), flist.end(), [](vector<int> A, vector<int> B){ return A[0] < B[0]; });
    }

	cout << "Writing mesh to - " << filename << endl;
	ofstream objFile;
	objFile.open(filename);
	objFile.precision(std::numeric_limits<double>::digits10 + 1);

    // Get current time.
    std::string export_time;
    char stime[BUFFER_LENGTH] = {0};
    time_t now_time;
    time(&now_time);
    strftime(stime,sizeof(stime),"%H:%M:%S",localtime(&now_time));
    export_time = stime + '\0';

    // Header 
    objFile << "# TIGER Mesh converter. (c) 2021." << endl;
    objFile << "# Created File: " << export_time << endl;
    objFile << "# " << endl;
    objFile << "# object default" << endl;
    objFile << "# " << endl;
    objFile << endl;

    // Write points
    for(int i = 0; i < V.rows(); i++) {
        objFile << "v";
        for(int j = 0; j < V.cols(); j++) {
            objFile << " " << V(i, j);
        }
        objFile << endl;
    }
    objFile << "# " << V.rows() << " vertices" << endl << endl;

    // Write facets with groups
    int curGroup = INT_MIN;
    for(int i = 0; i < flist.size(); i++) {
        if(flist[i][0] != curGroup) {
            curGroup = flist[i][0];
            objFile << "g " << curGroup << endl;
        }
        objFile << "f";
        for(int j = 1; j < flist[i].size(); j++) {
            objFile << " " << flist[i][j] + 1;
        }
        objFile << endl;
    }

    // Write facets 
    objFile << "# " << F.rows() << " faces" << endl << endl;

    objFile.close();

    return 1;
}
