// polyvec
#include <polyvec/polygon-tracer/detect_perfect_lines_inner_points.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/matrix.hpp>

using namespace Eigen;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int detect_perfect_lines_inner_points(
	const mat2x& B,
	vecXi& L,
	const bool circular
) {
	// The logic below is correct for non-circular paths, but the failure conditions
	// check for ring-2 neighbors while some lines don't necessarily need them!!!
	assert(circular);

	int lines = 0;

	std::vector<int> C;
	PathUtils::compute_convexities(B, C);

	Index ring[5];
	int dist[5];

	for (Index i = 0; i < B.cols(); ++i) {
		if (C[i] == 0) {
			continue;
		}

		// ++readability
		ring[2] = i;
		dist[2] = 0;

		ring[1] = PathUtils::next_transition_vertex(C, ring[2], -1, circular);
		ring[3] = PathUtils::next_transition_vertex(C, ring[2], +1, circular);

		if (ring[1] == -1 || ring[3] == -1) {
			continue;
		}
		
		ring[0] = PathUtils::next_transition_vertex(C, ring[1], -1, circular);
		ring[4] = PathUtils::next_transition_vertex(C, ring[3], +1, circular);

		if (ring[0] == -1 || ring[4] == -1) {
			continue;
		}

		dist[0] = CircularDist(B, ring[0], ring[1]);
		dist[1] = CircularDist(B, ring[1], ring[2]);
		dist[2] = CircularDist(B, ring[2], ring[2]);
		dist[3] = CircularDist(B, ring[2], ring[3]);
		dist[4] = CircularDist(B, ring[3], ring[4]);

		for (int j = 0; j < 5; ++j) {
			PF_ASSERT(C[ring[j]] != 0);

			if (j != 2) {
				PF_ASSERT(dist[j] > 0);
			}
		}

		// this could be checked more elegantly (todo)
		// a lovely example of data-driven programming. The number of conditions and values
		// should be symmetrical, otherwise:
		// i) some conditions are overlapping (best case scenario)
		// ii) some cases are not handled
		if (C[ring[0]] == -1 && C[ring[1]] == +1 && C[ring[2]] == -1 && C[ring[3]] == +1 && C[ring[4]] == -1) {
			if (dist[0] == dist[3] && dist[1] == 1 && dist[4] == 1 && dist[2] == 0) {
				MatrixUtils::append(L, i);
			}
		}

		if (C[ring[0]] == +1 && C[ring[1]] == -1 && C[ring[2]] == +1 && C[ring[3]] == -1 && C[ring[4]] == +1) {
			if (dist[1] == dist[4] && dist[0] == 1 && dist[3] == 1 && dist[2] == 0) {
				MatrixUtils::append(L, i);
			}
		}
		
		if (C[ring[1]] == +1 && C[ring[2]] == -1 && C[ring[3]] == +1 && C[ring[4]] == -1) {
			if (dist[1] == dist[4] && dist[1] > 1 && dist[3] == 1) {
				MatrixUtils::append(L, i);
			}
		}

		if (C[ring[0]] == +1 && C[ring[1]] == -1 && C[ring[2]] == +1 && C[ring[3]] == -1) {
			if (dist[0] == dist[3] && dist[0] > 1 && dist[1] == 1) {
				MatrixUtils::append(L, i);
			}
		}

		if (C[ring[0]] == +1 && C[ring[1]] == -1 && C[ring[2]] == +1 && C[ring[3]] == -1 && C[ring[4]] == +1) {
			if (dist[0] == 1 && dist[1] == dist[4] && dist[3] == 1) {
				MatrixUtils::append(L, i);
			}
		}

		if (C[ring[1]] == -1 && C[ring[2]] == +1 && C[ring[3]] == -1 && C[ring[4]] == +1) {
			if (dist[1] == dist[4] && dist[1] > 1 && dist[3] == 1) {
				MatrixUtils::append(L, i);
			}
		}

		if (C[ring[0]] == -1 && C[ring[1]] == +1 && C[ring[2]] == -1 && C[ring[3]] == +1) {
			if (dist[0] == dist[3] && dist[0] > 1 && dist[1] == 1) {
				MatrixUtils::append(L, i);
			}
		}
	}

	return lines;
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)