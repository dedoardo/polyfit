// Polyvec
#include <polyvec/regularity/collapse_opposite_short_segments.hpp>
#include <polyvec/utils/matrix.hpp>

using namespace Eigen;
using namespace polyfit;
using namespace polyfit::Regularity;

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(Regularity)

void collapse_opposite_short_segments(
	const polyfit::mat2x& B,
	polyfit::vecXi& P,
	const polyfit::Regularity::RegularityInformation& RE,
	const bool circular
) {
	mat3xi candidates;

	for(auto& r : RE.continuations()) {		
		const int i = r.v0;
		const int j = r.v1;

		candidates = mat3xi(3, 0);
		
		if (circular || i > 0) {
			MatrixUtils::append(candidates, vec3i(Circular(P, i - 1), i, j));
		}

		if (circular || i < P.size() - 1) {
			MatrixUtils::append(candidates, vec3i(i, Circular(P, i + 1), j));
		}

		if (circular || j > 0) {
			MatrixUtils::append(candidates, vec3i(Circular(P, j - 1), j, i));
		}

		if (circular || j < P.size() - 1) {
			MatrixUtils::append(candidates, vec3i(j, Circular(P, j + 1), i));
		}

		for (int k = 0; k < candidates.cols(); ++k) {

		}
	}
}

NAMESPACE_END(Regularity)
NAMESPACE_END(polyvec)