// polyvec
#include <polyvec/regularity/continuations_2018.hpp>
#include <polyvec/regularity/are_points_facing_inside.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/utils/matrix.hpp>

using namespace polyvec;
using namespace polyfit;
using namespace std;
using namespace Eigen;

NAMESPACE_BEGIN(polyvec)

std::vector<polyfit::Regularity::Continuation> find_continuation_candidates_2018(
    const Eigen::Matrix2Xd& B, 
    const Eigen::VectorXi& V,
    polyfit::Regularity::RegularityInformation& reg,
    const bool circular,
    const bool allow_move
) {
    Matrix2Xd P;
    BoundaryGraph::trace_to_points(B, V, P);

    vector<int> C;
    PathUtils::compute_convexities(P, C);

    vector<polyfit::Regularity::Continuation> continuations;

    vector<vector<double>> distances(V.size());
    for (Index i = 0; i < V.size(); ++i) {
        distances[i].resize(V.size(), INFINITY);
        
        for (Index ioff = 2; ioff < V.size() - 2; ++ioff) {
            const Index j = Circular(V, i + ioff);
            if (C[i] != -1 || C[j] != -1) {
                continue;
            }

            const double angle_at_i = AngleUtils::spanned_shortest(CircularAt(P, i - 1), P.col(i), CircularAt(P, i + 1));
            const double angle_at_j = AngleUtils::spanned_shortest(CircularAt(P, j - 1), P.col(j), CircularAt(P, j + 1));
            if (angle_at_i > PF_RAD(150) || angle_at_j > PF_RAD(150)) {
                continue;
            }

            distances[i][j] = (P.col(i) - P.col(j)).norm();
        }
    }

    unordered_set<int> matched;
    for (Index i = 0; i < V.size(); ++i) {
        int min_for_i = distance(distances[i].begin(), min_element(distances[i].begin(), distances[i].end()));
        
        if (distances[i][min_for_i] == INFINITY) {
            continue;
        }

        int min_for_j = distance(distances[min_for_i].begin(), min_element(distances[min_for_i].begin(), distances[min_for_i].end()));
        if (distances[min_for_i][min_for_j] == INFINITY || min_for_j != i) {
            continue;
        }

        auto j = min_for_i;

        const auto pi = B.col(V(i));
        const auto pj = B.col(V(j));

        // Test1: Inscribed circle
        const auto c = .5 * (pi + pj);
        const auto r = .5 * (pi - pj).norm() * .8;
        
        // test if a bunch of points around i j 
        bool is_inside = false;

        const int neighbors_to_test = 4;
        for (int nhb_off = -neighbors_to_test; !is_inside && nhb_off <= neighbors_to_test; ++nhb_off) {
            const int nhb = Circular(V, i + nhb_off);
            if (nhb == j) {
                break;
            }

            const double dist_sq = (B.col(V(nhb)) - c).norm();
            if (dist_sq + PF_EPS_MEDIUM < r) {
                is_inside = true;
                break;
            }
        }

        for (int nhb_off = -neighbors_to_test; !is_inside && nhb_off <= neighbors_to_test; ++nhb_off) {
            const int nhb = Circular(V, i + nhb_off);
            if (nhb == i) {
                break;
            }

            const double dist_sq = (B.col(V(nhb)) - c).norm();
            if (dist_sq + PF_EPS_MEDIUM < r) {
                is_inside = true;
                break;
            }
        }

        polyfit::Regularity::Continuation continuation;
        continuation.v0 = i;
        continuation.v1 = j;
        if (!is_inside) {
            
            // which one is the previous point? the most continuous one
            const double angle_i_prev = AngleUtils::spanned_shortest(CircularAt(P, i - 1), pi, pj);
            const double angle_i_next = AngleUtils::spanned_shortest(CircularAt(P, i + 1), pi, pj);
            const double angle_j_prev = AngleUtils::spanned_shortest(CircularAt(P, j - 1), pj, pi);
            const double angle_j_next = AngleUtils::spanned_shortest(CircularAt(P, j + 1), pj, pi);

            continuation.v0_prev = angle_i_prev < angle_i_next ? Circular(V, i + 1) : Circular(V, i - 1);
            continuation.v1_next = angle_j_prev < angle_j_next ? Circular(V, j + 1) : Circular(V, j - 1);

            continuations.emplace_back(continuation);
            matched.insert(i);
            matched.insert(j);
            continue; // not doing the second test
        }

        // test 1, same line continuation
        // need to test with the appropriate order
        const double angle_fw_i = AngleUtils::spanned_shortest(CircularAt(P, i - 1), pi, pj);
        const double angle_fw_j = AngleUtils::spanned_shortest(CircularAt(P, j + 1), pj, pi);
        if (angle_fw_i > PF_RAD(159) && angle_fw_j > PF_RAD(159)) {
            continuation.v0_prev = Circular(V, i - 1);
            continuation.v1_next = Circular(V, j + 1);
            continuations.emplace_back(continuation);
            matched.insert(i);
            matched.insert(j);
            continue;
        }

        const double angle_bw_i = AngleUtils::spanned_shortest(CircularAt(P, i + 1), pi, pj);
        const double angle_bw_j = AngleUtils::spanned_shortest(CircularAt(P, j - 1), pj, pi);
        if (angle_bw_i > PF_RAD(159) && angle_bw_j > PF_RAD(159)) {
            continuation.v0_prev = Circular(V, i + 1);
            continuation.v1_next = Circular(V, j - 1);
            continuations.emplace_back(continuation);
            matched.insert(i);
            matched.insert(j);
            continue;
        }
    }

    for (auto& continuation : continuations) {
        reg.add(move(continuation));
    }
    continuations.clear();
    return continuations;
}

NAMESPACE_END(polyvec)