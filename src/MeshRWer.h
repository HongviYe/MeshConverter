#ifndef _MESHRWER_H_
#define _MESHRWER_H_
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

using std::vector;
using std::cout;
using std::endl;
class TriMeshRWer {
public:
	TriMeshRWer();

	void Init(int num_t,int num_v) {
		V.resize(num_v);
		T.resize(num_t);

		std::fstream;
	}
	virtual std::string WhichFileType() = 0;

	virtual void Read(std::basic_istream< char, std::char_traits<char >>& stream) = 0;
	virtual void Write(std::basic_ostream< char, std::char_traits<char >>& stream) = 0;

	void readCoordinate(std::basic_istream< char, std::char_traits<char >>& stream, int index);
	void readConnection(std::basic_istream< char, std::char_traits<char >>& stream, int index);

	vector<int>& getConnection(int index);
	int getMark(int index);
	vector<double>& getCoordinate(int index);


	int getVertexSize();
	int getConnSize();
	int getblockMarkSize();
protected:
	vector<vector<double>> V;
	vector<vector<int>> T;
	vector<int> blockMark;
};


/**
 * @brief Seperate string origin by given a set of patterns.
 *
 * @param origin
 * @param patterns If it meets one of the patterns, delete the charactor and split it from this index.
 * @return std::vector<std::string>
 */
static std::vector<std::string> seperate_string(std::string origin, std::vector<std::string> patterns = { " ", "\t" }) {
	std::vector<std::string> result;
	if (origin.length() == 0) {
		return result;
	}
	origin += patterns[0];
	size_t startPos = origin.npos;
	for (auto pt = patterns.begin(); pt != patterns.end(); pt++) {
		startPos = (origin.find(*pt) < startPos) ? origin.find(*pt) : startPos;
	}
	size_t pos = startPos;
	while (pos != origin.npos) {
		std::string temp = origin.substr(0, pos);
		if (temp.length())
			result.push_back(temp);
		origin = origin.substr(pos + 1, origin.size());
		pos = origin.npos;
		for (auto pt = patterns.begin(); pt != patterns.end(); pt++) {
			pos = (origin.find(*pt) < pos) ? origin.find(*pt) : pos;
		}
	}
	return result;
}


#endif