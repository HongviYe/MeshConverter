#ifndef _PLS_RWER_H_ 
#define  _PLS_RWER_H_ 
#include "../include/MeshRWer.h"
#include <sstream>
#include <iostream>
class PLSRWer :public TriMeshRWer {
public:
	PLSRWer() {}

	virtual std::string WhichFileType() { return "pls"; };

	virtual void Read(std::basic_istream< char, std::char_traits<char >>& stream) {


		int nPoints = 0;
		int nFacets = 0;

		auto& pls_file = stream;
		std::string str;
		std::getline(pls_file, str);
		std::stringstream ss(str);
		ss >> nFacets >> nPoints;
		for (int i = 0; i < nPoints; i++) {
			int index;
			vector<double> point_coordinate(3);
			pls_file >>index>> point_coordinate[0]>> point_coordinate[1]>> point_coordinate[2];
			V.push_back(point_coordinate);
		}
		blockMark.resize(nFacets);
		for (int i = 0; i < nFacets; i++) {
			int index;
			vector<int> tri(3);
			pls_file >> index >> tri[0] >> tri[1] >> tri[2] >> blockMark[i];
			tri[0]--; tri[1]--; tri[2]--;
			T.push_back(tri);
		}

	};
	virtual void Write(std::basic_ostream< char, std::char_traits<char >>& stream) {
		auto &f = stream;
		f.precision(std::numeric_limits<double>::digits10 + 1);
		f << T.size() << " " << V.size() << " 0 0 0 0" << std::endl;
		for (int i = 0; i < V.size(); i++) {
			f << i + 1 << " " << V[i][0] << " " << V[i][1] << " " << V[i][2] << std::endl;
		}
		for (int i = 0; i < T.size(); i++) {
			f << 3 << " " << T[i][0] + 1 << " " << T[i][1] + 1 << " " << T[i][2] + 1;
			if (!blockMark.empty())
				f << " " << blockMark[i];
			else
				f << " " << 0;
			f << std::endl;
		}
	};

};
#endif // !_VTK_RWER_H_ 
