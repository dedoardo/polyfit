// polyvec
#include <polyvec/regularity/count_inflections_outside_subpath.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/core/log.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

int count_inflections_outside_subpath(
	const mat2x& P,
	const vecXi& V,
	const vec2i& B,
	bool circular
	) {

#if 0
	const Vertex ioff = circular ? 0 : 1;

	// which are the points that should be ignored?
	std::vector<bool> Vignore(V.size(), false);
	for (Vertex i = 0; i < V.size(); ++i) {
		if (!circular && (i == 0 || i == V.size() - 1)) {
			continue;
		}

		if (PathUtils::contains_closed(P, B(0), B(1), V[i], circular)) {
			Vignore[i] = true;
		}
	}

	int count = 0;

	for (Vertex i = ioff; i < V.size() - ioff; ++i) {
		const Vertex v0 = CircularAt(V, i - 1);
		const Vertex v1 = V[i];
		const Vertex v2 = CircularAt(V, i + 1);
		const Vertex v3 = CircularAt(V, i + 2);

		if (Vignore[v1] || Vignore[v2]) {
			continue;
		}

		count += AngleUtils::have_opposite_convexity(P.col(v0), P.col(v1), P.col(v2), P.col(v3));
	}

	return count;

#endif
	return 0;
}

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)