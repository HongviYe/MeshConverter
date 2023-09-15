#ifndef FACET_CLASSIFICATION_H
#define FACET_CLASSIFICATION_H
#include <Eigen/Core>
#include <vector>
#include <queue>

#include "igl/triangle_triangle_adjacency.h"
#include "igl/per_face_normals.h"
namespace MESHALG
{
	// Compute connected components of facets based on edge-edge angle.
	//
	// For connected components on vertices see igl::vertex_components
	//
	// Inputs:
	//   V  #V by 3 coordinates
	//   F  #F by 3 list of triangle indices
	//   angle #the target angle between two facets
	// Outputs:
	//   C  #F list of connected component ids
	// example facet_classification(V,F,30,C);

	template<typename DerivedV, typename DerivedF, typename DerivedC>
	void facet_classification(const Eigen::MatrixBase<DerivedV>& V, const Eigen::MatrixBase<DerivedF>& F, double angle, Eigen::PlainObjectBase<DerivedC>& C)
	{
		double cos_th = cos(angle / 180.0 * 3.14159265358);
		std::vector<std::vector<std::vector<int>>> facet_adjacency;
		Eigen::MatrixXd N;
		igl::per_face_normals(V, F, N);
		igl::triangle_triangle_adjacency(F, facet_adjacency);
		int start_faceid = 0;
		C.resize(F.rows(), 1);
		C.setConstant(-1);
		for (int i = 0; i < F.rows(); i++) {
			if (C(i) >= 0)
				continue;
			int current_facet = start_faceid;
			start_faceid++;
			std::queue<int> Q;
			Q.push(i);
			C(i) = current_facet;
			while (!Q.empty()) {
				int f = Q.front();				
				Q.pop();
				for (auto j : facet_adjacency[f]) {
					for (auto k : j) {
						if (C(k) < 0)
							if (N.row(k).dot(N.row(f)) > cos_th) {
								Q.push(k);
								C(k) = current_facet;
							}
					}
				}
			}
		}
		std::cout << "there are " << start_faceid << " colors" << endl;
		return;
	}
}

#endif //FACET_CLASSIFICATION_H