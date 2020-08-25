// polyvec
#include <polyvec/io/draw.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/debug.hpp>

#define POINT_RADIUS 1

using namespace polyvec;

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( Draw )

void indices(const mat2x& B) {
	for (Vertex v = 0; v < B.cols(); ++v) {
		draw::text(B.col(v), std::to_string(v), draw::font_pdf, Style::text());
	}
}

void indices(const mat2x& B, const vecXi& P) {
	for (Vertex v = 0; v < P.size(); ++v) {
		draw::text(B.col(P(v)), std::to_string(v), draw::font_pdf, Style::text());
	}
}

void boundary ( const mat2x& B, bool circular ) {
    for ( Vertex v = 0; v < B.cols(); ++v ) {
        if ( !circular && v == B.cols() - 1 ) {
            break;
        }

        draw::line ( B.col ( v ), CircularAt ( B, v + 1 ), Style::outline ( colors::gray, 6. ) );
    }
}

void graph(const mat2x& B, const std::vector<BoundaryGraph::Edge>& E) {
	for (size_t i = 0; i < E.size(); ++i) {
		draw::line(B.col(E[i].v0), B.col(E[i].v1), Style::outline(colors::red, 1.));
	}
}

void polygon(const mat2x& B, const vecXi& P, bool circular) {
	for (Vertex v = 0; v < P.size(); ++v) {
		if (!circular && v == P.cols() - 1) {
			break;
		}

		const Eigen::Index v0 = P(v);
		const Eigen::Index v1 = CircularAt(P, v + 1);

		draw::line(B.col(v0), B.col(v1), Style::outline(colors::black, 9.));
	}
}

void polygon ( const mat2x& P, bool circular ) {
    for ( Vertex v = 0; v < P.cols(); ++v ) {
        if ( !circular && v == P.cols() - 1 ) {
            break;
        }

        draw::line ( P.col ( v ), CircularAt ( P, v + 1 ), Style::outline ( colors::black, 9. ) );
    }
}

void polygon_fill ( const mat2x& P, bool circular ) {
	std::vector<vec2> points;
	for (int i = 0; i < P.cols(); ++i) {
		points.emplace_back(P.col(i));
	}
	draw::polygon(points, Style::fill(colors::black));
}

void curve(polyvec::GlobFitCurve* curve) {
	draw::curve(curve, 5., colors::calm_blue);
}

void closure_same_line(const mat2x& B, const vecXi& P, Vertex i, Vertex j, bool circular) {
	const vec2 pi = B.col(P(i));
	const vec2 pj = B.col(P(j));

	draw::point(pi, 5., Style::fill(colors::green));
	draw::point(pj, 5., Style::fill(colors::green));
	draw::line(pi, pj, Style::outline(colors::red, 1.));
}

void closure_same_circle(const mat2x& B, const vecXi& P, Vertex i, Vertex j, bool circular) {
	const vec2 pi = B.col(P(i));
	const vec2 pj = B.col(P(j));

	draw::point(pi, 5., Style::fill(colors::green));
	draw::point(pj, 5., Style::fill(colors::green));
	draw::line(pi, pj, Style::outline(colors::red, 1.));
}

NAMESPACE_END ( Draw )
NAMESPACE_END ( polyfit )