#ifndef _VTK_RWER_H_ 
#define  _VTK_RWER_H_ 
#include "../include/MeshRWer.h"
#include <sstream>
#include <iostream>
class VTKRWer :public TriMeshRWer {
public:

	virtual std::string WhichFileType() { return "vtk"; };

	virtual void Read(std::basic_istream< char, std::char_traits<char >>& stream) {

		int nPoints = 0;
		int nFacets = 0;
		auto& vtk_file = stream;
		std::string vtk_type_str = "POLYDATA ";
		char buffer[256];
		while (!vtk_file.eof()) {
			vtk_file.getline(buffer, 256);
			std::string line = (std::string)buffer;
			if (line.length() < 2 || buffer[0] == '#')
				continue;
			if (line.find("DATASET") != std::string::npos) {
				std::vector<std::string> words = seperate_string(line);
				if (words[1] == "POLYDATA")
					vtk_type_str = "POLYGONS ";
				else if (words[1] == "UNSTRUCTURED_GRID")
					vtk_type_str = "CELLS ";
				else {
					std::cout << "WARN : The format of VTK file is illegal, No clear DATASET name." << endl;
				}
			}
			if (line.find("POINTS ") != std::string::npos) {
				std::vector<std::string> words = seperate_string(line);
				nPoints = stoi(words[1]);
				for (int i = 0; i < nPoints; i++) {
					vtk_file.getline(buffer, 256);
					words = seperate_string(std::string(buffer));
					V.push_back({ stod(words[0]), stod(words[1]), stod(words[2]) });
				}
			}
			if (line.find(vtk_type_str) != std::string::npos) {
				std::vector<std::string> words = seperate_string(line);
				nFacets = stoi(words[1]);
				for (int i = 0; i < nFacets; i++) {
					vtk_file.getline(buffer, 256);
					words = seperate_string(std::string(buffer));
					T.push_back({ stoi(words[1]), stoi(words[2]), stoi(words[3]) });
				}
			}
		}
	};
	virtual void Write(std::basic_ostream< char, std::char_traits<char >>& stream) {
		auto &f = stream;
		f.precision(std::numeric_limits<double>::digits10 + 1);
		f << "# vtk DataFile Version 2.0" << std::endl;
		f << "TetWild Mesh" << std::endl;
		f << "ASCII" << std::endl;
		f << "DATASET UNSTRUCTURED_GRID" << std::endl;
		f << "POINTS " << V.size() << " double" << std::endl;
		for (int i = 0; i < V.size(); i++)
			f << V[i][0] << " " << V[i][1] << " " << V[i][2] << std::endl;
		f << "CELLS " << T.size() << " " << T.size() * 4 << std::endl;
		for (int i = 0; i < T.size(); i++) {
			f << T[i].size() << " ";
			for (int j = 0; j < T[i].size(); j++)
				f << T[i][j] << " ";
			f << std::endl;
		}
		f << "CELL_TYPES " << T.size() << std::endl;
		int cellType = 5;
		for (int i = 0; i < T.size(); i++)
			f << cellType << std::endl;
		if (blockMark.size() == T.size()) {
			f << "CELL_DATA " << T.size() << std::endl;
			f << "SCALARS cell2bodyid int 1" << std::endl;
			f << "LOOKUP_TABLE default" << std::endl;
			for (int i = 0; i < blockMark.size(); i++)
				f << blockMark[i] << std::endl;
		}
	};

};
#endif // !_VTK_RWER_H_ 

