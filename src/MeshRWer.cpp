#include "MeshRWer.h"
#include <exception>
TriMeshRWer::TriMeshRWer()
{
}
void TriMeshRWer::readCoordinate(std::basic_istream<char, std::char_traits<char>>& stream, int index)
{
	if (index >= V.size()||index<0)
		throw std::exception("exceed maximum index in reading coordinate");
	stream >> V[index][0]>> V[index][1]>> V[index][2];
}

void TriMeshRWer::readConnection(std::basic_istream<char, std::char_traits<char>>& stream, int index)
{
	if (index >= T.size()|| index < 0)
		throw std::exception("exceed maximum index in reading connection");
	stream >> T[index][0] >> T[index][1] >> T[index][2];
}

vector<int>& TriMeshRWer::getConnection(int index)
{
	if(index >= T.size() || index < 0)
		throw std::exception("exceed maximum index in reading connection");
	return T[index];
}

int TriMeshRWer::getMark(int index)
{
	if (index >= blockMark.size() || index < 0)
		throw std::exception("exceed maximum index in reading blockMark");
	return blockMark[index];
}

vector<double>& TriMeshRWer::getCoordinate(int index)
{
	if (index >= V.size() || index < 0)
		throw std::exception("exceed maximum index in reading coordinate");
	return V[index];
}

int TriMeshRWer::getVertexSize()
{
	return V.size();
}

int TriMeshRWer::getConnSize()
{
	return T.size();
}

int TriMeshRWer::getblockMarkSize()
{
	return blockMark.size();
}
