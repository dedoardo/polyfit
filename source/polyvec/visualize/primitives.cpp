// polyvec
#include <polyvec/visualize/primitives.hpp>
#include <polyvec/debug.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/utils/matrix.hpp>

#define LABEL_FONT draw::font_pdf / 2.5, Style::text()
#define CURVE_STROKE_WIDTH 5.
#define CURVE_COLOR_ALT0 colors::talking_orange
#define CURVE_COLOR_ALT1 colors::forest_green
#define MIDPOINT_STYLE_ALT0 2., Style::fill(colors::talking_orange)
#define MIDPOINT_STYLE_ALT1 2., Style::fill(colors::forest_green)
#define NORMAL_STYLE_ALT0 Style::outline(colors::talking_orange, 2.)
#define NORMAL_STYLE_ALT1 Style::outline(colors::forest_green, 2.)

using namespace polyvec;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Visualize)

void primitives_indices(
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		const string label = StringUtils::fmt("%llu (%d)", i, prims[i].corner);
		draw::text(curve->get_curve()->pos(0.4), label, LABEL_FONT);
	}
}

void primitives_alternated(
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		draw::curve(curve->get_curve().get(), CURVE_STROKE_WIDTH, i % 2 ? CURVE_COLOR_ALT0 : CURVE_COLOR_ALT1);
	}
}

void primitives_alternated(
	const mat2x& P,
	bool circular
) {
	for (Vertex i = 0; i < P.cols(); ++i) {
		draw::line(P.col(i), CircularAt(P, i + 1), Style::outline(i % 2 ? CURVE_COLOR_ALT0 : CURVE_COLOR_ALT1, 2.));
	}
}

void primitives_midpoints(
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		for (size_t j = 0; j < prim.fitting_info.fit_midpoints.size(); ++j) {
			const vec2 p_mid = prim.fitting_info.fit_midpoints.at(j);

			if (i % 2) {
				draw::point(p_mid, MIDPOINT_STYLE_ALT0);
			}
			else {
				draw::point(p_mid, MIDPOINT_STYLE_ALT1);
			}
		}
	}
}

void primitives_midpoints_normals(
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		for (size_t j = 0; j < prim.fitting_info.fit_midpoints.size(); ++j) {
			const vec2 p_mid = prim.fitting_info.fit_midpoints.at(j);
			const vec2 n_mid = prim.fitting_info.fit_midpoint_normals.at(j).normalized();

			draw::line(p_mid, p_mid + n_mid, i % 2 ? NORMAL_STYLE_ALT0 : NORMAL_STYLE_ALT1);
		}
	}

}

NAMESPACE_END(Visualize)
NAMESPACE_END(polyfit)