// polyvec
#include <polyvec/regularity/is_subpath_valid.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/matrix.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Symmetry)

bool is_subpath_valid(
	const mat2x& P,
	const vec2i S0,
	const vec2i S1,
	const vecXi V,
	const bool circular
) {

	bool is_symmetric = true;

#if 0

	// the two symmetric regions should span the same length, regardless if
	// they are disjoint or sequential
	PF_ASSERT(CircularDist(P, S0(0), S0(1), circular) ==
			  CircularDist(P, S1(0), S1(1), circular));

	// Checking whether potentially symmetric points have the same distance
	// relatively to the bounds of the symmetric region
	for (Vertex i = 0; is_symmetric && i < V.size(); ++i) {
		// ignoring this point if it's outside the first symmetric region
		if (!PathUtils::contains_closed(P, S0(0), S0(1), V[i], circular)) {
			continue;
		}

		const Vertex d0 = CircularDist(P, S0(0), V[i]);

		// we expect to find another vertex which has complementary distance
		bool is_i_symmetric = false;
		for (Vertex j = 0; j < V.size(); ++j) {
			if (PathUtils::contains_closed(P, S(0), S1(1), V[j])) {
				const Vertex d1 = CircularDist(P, V[j], S1(1)) ;

				if (d1 == d0) {
					is_i_symmetric = true;
					break;
				}
			}
		}

		if (!is_i_symmetric) {
			is_symmetric = false;
		}
	}

#endif

	return is_symmetric;
}

NAMESPACE_END(Symmetry)
NAMESPACE_END(polyfit)
