// polyvec
#include <polyvec/regularity/add_90_degree_edges.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/utils/matrix.hpp>

// libc++
#include <unordered_map>
#include <algorithm>

using namespace Eigen;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

// Assumes that the graph construction never cuts off entire pixels.
class AxisAligned90DegreeRegularity
    : public RegularityAction
{
public:
    AxisAligned90DegreeRegularity(int v_aa_90) :
        v_aa_90(v_aa_90)
    { 
    }

    bool can_introduce_degenerate_edges() const { return false; }
    bool can_introduce_inflections() const { return false; }
    bool is_polygon_acceptable(const mat2x& raster, const vecXi& polygon_old, const vecXi& polygon_new, bool circular) const {
        // The corner that we were supposed to regularize was lost while trying to apply previous regularities.
        if (accuracy_old == INFINITY) {
            return false;
        }

        // Get the edge after the one which we have just regularized and check if the accuracy has improved
        PF_VERBOSE_F("is polygon acceptable %d", v_aa_90);
        for (Index i = 0; i < polygon_new.size(); ++i) {
            if (polygon_new(i) == v_aa_90) {
                int v_prev, v_next;
                if (direction < 0) {
                    v_prev = i - 2;
                    v_next = i - 1;
                }
                else {
                    v_prev = i + 1;
                    v_next = i + 2;
                }

                const vec2 p_prev = raster.col(CircularAt(polygon_new, v_prev));
                const vec2 p_next = raster.col(CircularAt(polygon_new, v_next));
                accuracy_new = GeomRaster::distance_bounds_from_points_with_slack(raster, p_prev, p_next, CircularAt(polygon_new, v_prev), CircularAt(polygon_new, v_next)).maxCoeff();
                PF_VERBOSE_F("accuracy_old %f accuracy_new %f direction %d vertex %d v_prev %d v_next %d", accuracy_old, accuracy_new, direction, v_aa_90, v_prev, v_next);
                return accuracy_new < accuracy_old + 0.1;
            }
        }

        return false;
    }

    std::vector<size_t> get_edges_to_delete(const mat2x& raster, const vecXi& polygon, const std::vector<BoundaryGraph::Edge>& E, bool circular) const
    {
        // Calculate the old accuracy
        for (Index i = 0; i < polygon.size(); ++i) {
            if (polygon(i) == v_aa_90) {
                const vec2 pp = raster.col(CircularAt(polygon, i - 1));
                const vec2 p = raster.col(polygon(i));
                const vec2 pn = raster.col(CircularAt(polygon, i + 1));

                if (!((p - pp).cwiseAbs().minCoeff() > PF_EPS_MEDIUM || (pn - p).cwiseAbs().minCoeff() > PF_EPS_MEDIUM)) {
                    return std::vector<size_t>();
                }

                if ((p - pp).cwiseAbs().minCoeff() < PF_EPS_MEDIUM) {
                    accuracy_old = GeomRaster::distance_bounds_from_points_with_slack(raster, p, pn, polygon(i), CircularAt(polygon, i + 1)).maxCoeff();
                    direction = 1;
                }else {
                    accuracy_old = GeomRaster::distance_bounds_from_points_with_slack(raster, pp, p, CircularAt(polygon, i -1), polygon(i)).maxCoeff();
                    direction = -1;
                }
                break;
            }
        }

        // Find the edges to be removed
        PF_ASSERT(v_aa_90 != -1);

        std::vector<size_t> result;
        PF_VERBOSE_F("get edges to delete for %d", v_aa_90);
        for (size_t i = 0; i < E.size(); ++i) {
            const auto& e = E[i];
            if (e.v0 == e.v1) {
                //PF_ABORT;
                continue;
            }

            if ((e.v1 == v_aa_90 || e.v0 == v_aa_90) && (raster.col(e.v0) - raster.col(e.v1)).cwiseAbs().minCoeff() > PF_EPS_MEDIUM) {
                PF_VERBOSE_F("remove edge %d %d", e.v0, e.v1);
                result.emplace_back(i);
            }
        }

        // Not regularizing a single pixel edge
        if (result.size() == 1) {
            if ((raster.col(E[result[0]].v0) - raster.col(E[result[0]].v1)).squaredNorm() < 1. + PF_EPS) {
                return std::vector<size_t>();
            }
        }
        
        return result;
    }

private:
    mutable int direction = 0;
    mutable double accuracy_old = INFINITY;
    mutable double accuracy_new = INFINITY;
    int v_aa_90 = -1;
};

struct RegularizeEdge {
    double length_sq;
    int vertex;
};

void add_90_degree_edges(
    const Eigen::Matrix2Xd& raster,
    const Eigen::VectorXi& polygon,
    const bool circular,
    std::vector<RegularityAction*>& regularity_actions
) {
    if (polygon.size() < 4) {
        return;
    }

    std::vector<int> C_raster;
    PathUtils::compute_convexities(raster, C_raster);

    for (Index i = 0; i < polygon.size(); ++i) {
        const vec2 pp = raster.col(CircularAt(polygon, i - 1));
        const vec2 p = raster.col(polygon(i));
        const vec2 pn = raster.col(CircularAt(polygon, i + 1));
        const vec2 pnn = raster.col(CircularAt(polygon, i + 2));

        PF_VERBOSE_F("polygon(i) = %d", polygon(i));

        // Skipping 1x1 pixels as the regularization will not be successful anyway
        const bool is_regularizable = abs(1. - (pn - p).cwiseAbs().minCoeff()) < PF_EPS && CircularDist(raster, polygon(i), CircularAt(polygon, i + 1)) > 2;
        if (!is_regularizable) {
            continue;
        }

        const bool is_prev_aa = (p - pp).cwiseAbs().minCoeff() < PF_EPS && (p - pp).squaredNorm() > 1. + PF_EPS;
        const bool is_next_aa = (pnn - pn).cwiseAbs().minCoeff() < PF_EPS && (pnn - pn).squaredNorm() > 1. + PF_EPS;
        
        if (!is_prev_aa && !is_next_aa) {
            continue;
        }

        // Handling conflicts when two edges want to create a 90 degree angle. We prefer the side with the longer pixel run as
        // it has the highest chance of succeeding
        bool make_prev_90;
        if (is_prev_aa && is_next_aa) {
            const int run_forward = PathUtils::distance_to_next_corner(C_raster, polygon(i), +1);
            const int run_backward = PathUtils::distance_to_next_corner(C_raster, CircularAt(polygon, i + 1), -1);
            
            // if run_forward == run_backward the regularization will fail as accuracy will get worse
            make_prev_90 = run_forward > run_backward;
        } else {
            make_prev_90 = is_prev_aa;
        }

        PF_VERBOSE_F("Regularizing %d %d %d %d make_prev %d", Circular(polygon, i - 1), i, Circular(polygon, i + 1), Circular(polygon, i + 2), make_prev_90);
        PF_VERBOSE_F("regularizing edge %d %d", polygon(i), CircularAt(polygon, i + 1));
        if (make_prev_90) {
            regularity_actions.emplace_back(new AxisAligned90DegreeRegularity(polygon(i)));
            PF_VERBOSE_F("Attempt to orthogonalize at vertex %d", polygon(i));
        } else {
            regularity_actions.emplace_back(new AxisAligned90DegreeRegularity(CircularAt(polygon, i + 1)));
            PF_VERBOSE_F("Attempt to orthogonalize at vertex %d", CircularAt(polygon, i + 1));
        }
    }
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)
