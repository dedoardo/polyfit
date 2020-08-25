// polyvec
#include <polyvec/polygon-tracer/error.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/utils/matrix.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

vec3 error_per_metric(const mat2x& B,
	const std::vector<BoundaryGraph::Edge>& E,
	const vecXi& P,
	const bool circular) {
	PF_ASSERT(circular); // todo

	// indices of each vertex
	std::vector<Vertex> VE;
	PF_ASSERT(BoundaryGraph::trace_to_edges(E, P, VE, circular));

	vec3 error = vec3::Zero();

	for (Vertex i = 0; i < (Vertex)VE.size(); ++i) {
		const BoundaryGraph::Edge& ep = E[CircularAt(VE, i - 1)];
		const BoundaryGraph::Edge& e = E[VE[i]];
		const BoundaryGraph::Edge& en = E[CircularAt(VE, i + 1)];

		const vec2 p[4] = {
			B.col(ep.v0),
			B.col(e.v0),
			B.col(e.v1),
			B.col(en.v1)
		};

		error += ErrorMetrics::calculate_inner_sep(e.d_min, e.d_max, p[0], p[1], p[2], p[3]);
	}

	return error;
}

double error_total(const mat2x& B,
	const std::vector<BoundaryGraph::Edge>& E,
	const vecXi& P,
	const bool circular) {
	return error_per_metric(B, E, P, circular).sum();
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)