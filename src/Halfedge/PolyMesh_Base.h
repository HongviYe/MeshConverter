#pragma once


#include "assert.h"
#include "Eigen/Dense"

namespace MESHIO {
namespace polymesh{
	
class MVert;
class MEdge;
class MHalfedge;
class MPolyFace;

/**
* the class of Vert in Half_edge structure
* all index is begin with 0;
*/
class MVert
{
private:
	int index_;

	Eigen::Vector2d uv_;
	Eigen::Vector3d point_;
	Eigen::Vector3d normal_;

	MHalfedge* he_;
public:
	MVert() :
		index_(-1), uv_(0, 0), point_(0, 0, 0), normal_(0, 0, 0), he_(nullptr)
	{

	}

	MVert(double x, double y, double z, double u, double v) :
		index_(-1), uv_(u, v), point_(x, y, z), normal_(0, 0, 0), he_(nullptr)
	{
	}

	~MVert() { index_ = -1; }

public:
	MHalfedge* const halfEdge() { return he_; }
	const MHalfedge* const halfEdge() const { return he_; }

	void setHalfedge(MHalfedge* he) { he_ = he; }

	void setPosition(Eigen::Vector3d new_point) { point_ = new_point; }
	void setPosition(double x, double y, double z) { point_ = Eigen::Vector3d(x, y, z); }
	Eigen::Vector3d position() { return point_; }
	const Eigen::Vector3d& position() const { return point_; }

	double u() const { return uv_.x(); }
	double v() const { return uv_.y(); }

	double x() const { return point_.x(); }
	double y() const { return point_.y(); }
	double z() const { return point_.z(); }

	double nx() const { return normal_.x(); }
	double ny() const { return normal_.y(); }
	double nz() const { return normal_.z(); }

	void setNormal(Eigen::Vector3d& new_vec) { normal_ = new_vec; }
	void setNormal(double nx, double ny, double nz) { normal_ = Eigen::Vector3d(nx, ny, nz); }
	Eigen::Vector3d normal() { return normal_; }

	int index() const { return index_; }
	void set_index(int index) { index_ = index; }

	bool isIsolated() const { return he_ == nullptr; }

	//Try to set the half of the point as the boundary half. Be sure call it before adding faces by yourself.
	void adjustOutgoingHalfedge();
};

/**
* the class of Edge in Half_edge structure;
* the halfedge index is not from the edge index
* One can not select a halfedge or see a halfedge;
*/
class MHalfedge
{
private:
	int index_;

	MVert* v_;
	MEdge* e_;
	MPolyFace* poly_face_;

	MHalfedge* next_, * prev_;

	MHalfedge* pair_;

public:
	MHalfedge() : index_(-1), next_(nullptr), prev_(nullptr), pair_(nullptr), 
		v_(nullptr), e_(nullptr), poly_face_(nullptr) 
	{
	}
	MHalfedge(MHalfedge* next, MHalfedge* prev, MHalfedge* pair, MVert* v, MEdge* e, MPolyFace* p)
		: index_(-1), next_(next), prev_(prev), pair_(pair), v_(v), e_(e), poly_face_(p) 
	{
	}

	~MHalfedge() { index_ = -1; };
public:
	MHalfedge* const next() { return next_; }
	MHalfedge* const prev() { return prev_; }
	MHalfedge* const pair() { return pair_; }
	MVert* const fromVertex() { return v_; }
	MVert* const toVertex() { return next()->fromVertex(); }
	MEdge* const edge() { return e_; }
	MPolyFace* const polygon() { return poly_face_; }

	MHalfedge* const rotateNext() { return pair()->next(); }
	MHalfedge* const rotatePrev() { return prev()->pair(); }

	const MHalfedge* const next() const { return next_; }
	const MHalfedge* const prev() const { return prev_; }
	const MHalfedge* const pair() const { return pair_; }
	const MVert* const fromVertex() const { return v_; }
	const MVert* const toVertex() const { return next()->fromVertex(); }
	const MEdge* const edge() const { return e_; }
	const MPolyFace* const polygon() const { return poly_face_; }

	void setNext(MHalfedge* next) { next_ = next; }
	void setPrev(MHalfedge* prev) { prev_ = prev; }
	void setPair(MHalfedge* pair) { pair_ = pair; }
	void setVert(MVert* vert) { v_ = vert; }
	void setEdge(MEdge* edge) { e_ = edge; }
	void setPolygon(MPolyFace* poly) { poly_face_ = poly; }

	bool isBoundary() const { return poly_face_ == nullptr; }

	int index() { return index_; }
	//int edge_index() { return index_ / 2; }
	void set_index(int index) { index_ = index; }

	///get the direction of the edge, from fromVectex to toVertex;
	Eigen::Vector3d tangent()
	{
		Eigen::Vector3d t = toVertex()->position() - fromVertex()->position();
		t.normalize();
		return t;
	}


};

/**
* the class of Edge in Half_edge structure
*/
class MEdge
{
private:
	int index_;

	MVert* v1_; MVert* v2_;

	MHalfedge* he_;
public:
	MEdge() : index_(-1), v1_(nullptr), v2_(nullptr), he_(nullptr) {}
	MEdge(MVert* v1, MVert* v2) : index_(-1), v1_(v1), v2_(v2), he_(nullptr) {}
	MEdge(MVert* v1, MVert* v2, MHalfedge* he) : index_(-1), v1_(v1), v2_(v2), he_(he){}

	~MEdge() { index_ = -1; };

public:
	MHalfedge* const halfEdge() { return he_; }
	const MHalfedge* const halfEdge() const { return const_cast<MEdge*>(this)->halfEdge(); }

	void setHalfedge(MHalfedge* he) { he_ = he; }
	void setVert(MVert* v1, MVert* v2) { v1_ = v1; v2_ = v2; }
	void updateVert() { v1_ = he_->fromVertex(), v2_ = he_->toVertex(); }

	int index() const { return index_; }
	void set_index(int index) { index_ = index; }


	///get Vertex of the edge, the 0 is the first, the 1 is the second, the return is not orderd;
	MVert* getVert(int edge_v)
	{
		updateVert();
		if (edge_v == 0) return v1_;
		else if (edge_v == 1) return v2_;
		else return nullptr;
	}

	const MVert* getVert(int edge_v) const 
	{
		return const_cast<MEdge*>(this)->getVert(edge_v);
	}

	double length() 
	{
		updateVert();
		Eigen::Vector3d t = v1_->position() - v2_->position();
		return t.norm();
	}

	Eigen::Vector3d getCenter()
	{
		updateVert();
		return v1_->position() * 0.5 + v2_->position() * 0.5;
	}

};

/**
* the class of PolyFace in Half_edge structure;
* all index is begin with 0
*/
class MPolyFace
{
private:
	int index_;

	MHalfedge* he_begin_;

	Eigen::Vector3d normal_;		/* face normal */
	int polynum_;			/* number of vertices in the face */
	int material_index_;	/* material index, 暂时替代OpenMesh中的texindex */

public:
	MPolyFace() : 
		index_(-1), he_begin_(nullptr), normal_(0,0,0), polynum_(0), material_index_(0)
	{
	}
	MPolyFace(MHalfedge* he) : 
		index_(-1), he_begin_(he), normal_(0, 0, 0), polynum_(0), material_index_(0)
	{
	}
	~MPolyFace() { index_ = -1; }

public:
	MHalfedge* const halfEdge() { return he_begin_; }
	const MHalfedge* const halfEdge() const { return const_cast<MPolyFace*>(this)->halfEdge(); }

	void setHalfedge(MHalfedge* he) { he_begin_ = he; }

	unsigned int PolyNum() { 
		updatePolyNum();
		return polynum_; 
	}

	//be sure you have updatePoluNum
	unsigned int PolyNum() const { 
		return polynum_;
	}

	void setNormal(Eigen::Vector3d new_vec) { normal_ = new_vec; }
	void setNormal(double nx, double ny, double nz) { normal_ = Eigen::Vector3d(nx, ny, nz); }
	Eigen::Vector3d normal() { return normal_; }

	int index() const { return index_; }
	void set_index(int index) { index_ = index; }

	int updatePolyNum();
	void updataNormal();
	void update();				/* update polynum and normal */

	Eigen::Vector3d getFaceCenter();


};

}//polymesh
}//SmeshGen