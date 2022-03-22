// Polyvec
#include <polyvec/curve-tracer/refinement.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/utils/matrix.hpp>

// libc++
#include <algorithm>

using namespace Eigen;
using namespace polyfit;

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(CurveTracer)

int collapse_asymmetric_constant_steps(
	const polyfit::mat2x& B,
	polyfit::vecXi& P,
	const polyfit::Regularity::RegularityInformation& RE,
	std::vector<TangentFitType>& tangent_fits,
	const bool circular
) {
	PF_ASSERT(P.size() == tangent_fits.size());

	std::vector<int> C;
	PathUtils::compute_convexities(B, C);

	const double collapse_range = 1. + PF_EPS;

	int collapsed_tot = 0;
	int collapsed;

	do {
		collapsed = 0;

		for (size_t i = 0; i < tangent_fits.size(); ++i) {
			if (!circular && i == tangent_fits.size() - 1) {
				continue;
			}

			if (tangent_fits[i] != TANGENT_FIT_CONSTANT || CircularAt(tangent_fits, i + 1) != TANGENT_FIT_CONSTANT) {
				continue;
			}

			const Index j = Circular(P, i + 1);

			const Index vi = P(i);
			const Index vj = P(j);

			const vec2 pi = B.col(vi);
			const vec2 pj = B.col(vj);

			const int convexityi = AngleUtils::convexity(B.col(CircularAt(P, i - 1)), pi, pj, 1);
			const int convexityj = AngleUtils::convexity(pi, pj, B.col(CircularAt(P, j + 1)), 1);

			if (convexityi != convexityj) {
				continue;
			}

			const bool axis_aligned = (pi - pj).cwiseAbs().minCoeff() < PF_EPS;
			const bool within_range = (pi - pj).cwiseAbs().maxCoeff() < collapse_range;	

			if (!axis_aligned || !within_range) {
				continue;
			}

			//fprintf(stderr, "Collapse %d(%d) %d(%d)\n", i, vi, Circular(P, i + 1), vj);
			//fprintf(stderr, "dist %f %\n ", (pi - pj)(0), (pi - pj)(1));

			Index i_opp = -1, j_opp = -1;
			int i_opp_sym_size = 0, j_opp_sym_size = 0;
			int i_regularity = -1, j_regularity = -1;

			// Looking for the symmetric small edge
			for(auto& sym : RE.vertex_symmetries()) {			
				//printf("regularity: %d <-> %d size %d\n", sym.v0, sym.v1, sym.size);

				if (sym.v0 == i && sym.size > i_opp_sym_size) {
					i_opp = sym.v1;
					i_opp_sym_size = sym.size;
					i_regularity = sym.region;
				}

				if (sym.v0 == j && sym.size > j_opp_sym_size) {
					j_opp = sym.v1;
					j_opp_sym_size = sym.size;
					j_regularity = sym.region;
				}

				if (sym.v1 == i && sym.size > i_opp_sym_size) {
					i_opp = sym.v0;
					i_opp_sym_size = sym.size;
					i_regularity = sym.region;
				}

				if (sym.v1 == j && sym.size > j_opp_sym_size) {
					j_opp = sym.v0;
					j_opp_sym_size = sym.size;
					j_regularity = sym.region;
				}

				break;				
			}

			bool is_opp_valid = i_regularity == j_regularity && i_opp != -1 && j_opp != -1;

			if (i_opp == j || j_opp == i) {
				is_opp_valid = false;
			}

			const int vdist_prev = CircularDist(B, PathUtils::next_transition_vertex(C, vi, -1, circular), vi);
			const int vdist_next = CircularDist(B, vj, PathUtils::next_transition_vertex(C, vj, +1, circular));
			printf("dist to %d : %d\n", vi, vdist_prev);
			printf("dist from %d : %d\n", vj, vdist_next);

			Index vnew, vnew_opp = -1;
			if (vdist_prev > vdist_next) {
				vnew = vi;
				vnew_opp = is_opp_valid ? P(i_opp) : -1;
			}
			else if (vdist_prev < vdist_next) {
				vnew = vj;
				vnew_opp = is_opp_valid ? P(j_opp) : -1;
			}
			else {
				PF_ASSERT((vi - vj) % 2 == 0);
				PF_ASSERT((P(i_opp) - P(j_opp)) % 2 == 0);
				vnew = (vi + vj) / 2;
				vnew_opp = is_opp_valid ? ((P(i_opp) + P(j_opp)) / 2) : -1;
			}

			PF_ASSERT(vnew != -1);
			printf("relocating %d: %d -> %d\n", i, P(i), vnew);
			P(i) = vnew;
			tangent_fits[i] = TANGENT_FIT_CONSTANT;

			if (is_opp_valid) {
				PF_ASSERT(vnew_opp != -1);
				printf("relocating %d: %d -> %d\n", i_opp, P(i_opp), vnew_opp);
				P(i_opp) = vnew_opp;
				tangent_fits[i_opp] = TANGENT_FIT_CONSTANT;
			}

			std::vector<size_t> delete_list;
			delete_list.emplace_back((size_t)j);

			if (is_opp_valid) {
				delete_list.emplace_back((size_t)j_opp);
			}

			std::sort(delete_list.begin(), delete_list.end());
			EraseOrdered(P, delete_list);
			EraseOrdered(tangent_fits, delete_list);
			++collapsed;
			break;
		}

		collapsed_tot += collapsed;
	} while (collapsed);

	return collapsed_tot;
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyvec)