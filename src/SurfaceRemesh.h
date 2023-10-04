#include "remesh.h"

namespace SurfaceRemesh
{
	void iso_remesh(Mesh& mesh);

	void iso_split_long_edges(PolyMesh* mesh, double high);

	void iso_collapse_short_edges(PolyMesh* mesh, double high, double low);

	void iso_equalize_valences(PolyMesh* mesh);

	void iso_tangential_relaxation(PolyMesh* mesh);

	double iso_calculateTargetEdgeLength(PolyMesh* mesh);

	double iso_calculate_target_lenth_shortest(PolyMesh* mesh);
}



