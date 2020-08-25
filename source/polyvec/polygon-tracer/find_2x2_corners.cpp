// polyvec
#include <polyvec/polygon-tracer/find_2x2_corners.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/angle.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;
using namespace Eigen;

#define REQUIRE_CONTINUATION 1
#define ONLY_CONCAVE_CORNERS 0

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

void find_2x2_corners(
	const mat2x& B, // grid aligned points
	vecXi& C,       // corners
	const bool circular
) {
	C.resize(0);

	struct RasterSet
	{
		std::vector<int> points;
		bool sorted = false;
		int dimension;

		void assert_sorted(const mat2x& B)
		{
			if (!sorted)
				std::sort(points.begin(), points.end(), [&](int i, int j) { return B.coeff(dimension, i) < B.coeff(dimension, j); });
			sorted = true;
		}
	};

	auto get_id = [](const vec2& point, int dim) { return (int)std::round(2 * point(dim)); };

	// All 2x2 corners grouped by x/y coordinates
	std::map<int, RasterSet> raster_points_per_x;
	std::map<int, RasterSet> raster_points_per_y;

	std::vector<bool> is_corner(B.size(), false);
	std::vector<int> corners;

#if ONLY_CONCAVE_CORNERS
	std::vector<int> convexities;
	PathUtils::compute_convexities(B, convexities);
#endif

	for (size_t i = 0; i < B.cols(); ++i) {
		if (!circular && (i == 0 || i == B.size() - 1)) {
			continue;
		}

		raster_points_per_x[get_id(B.col(i), 0)].points.push_back(i);
		raster_points_per_y[get_id(B.col(i), 1)].points.push_back(i);

#if ONLY_CONCAVE_CORNERS
		if (convexities[i] != -1)
			continue;
#endif

		auto corner_size = GeomRaster::raster_corner_size(B, i);	

		if (corner_size >= 2) {
			is_corner[i] = true;
			corners.push_back(i);
		}
	}

	// Now check if the corners are raster continuations.
	for (auto i : corners)
	{
#if REQUIRE_CONTINUATION
		int x = get_id(B.col(i), 0);
		int y = get_id(B.col(i), 1);

		auto& x_raster = raster_points_per_x.at(x);
		x_raster.dimension = 1;
		x_raster.assert_sorted(B);

		auto& y_raster = raster_points_per_y.at(y);
		y_raster.dimension = 0;
		y_raster.assert_sorted(B);

		// find the current point in the raster sets
		auto in_x = std::lower_bound(x_raster.points.begin(), x_raster.points.end(), i, [&](int i, int j) { return B.coeff(1, i) < B.coeff(1, j); });
		auto in_y = std::lower_bound(y_raster.points.begin(), y_raster.points.end(), i, [&](int i, int j) { return B.coeff(0, i) < B.coeff(0, j); });

		bool is_continuation = false;
		if (std::next(in_x) != x_raster.points.end() && is_corner[*std::next(in_x)] && std::abs(*std::next(in_x) - i) > 1)
			is_continuation = true;
		if (in_x != x_raster.points.begin() && is_corner[*std::prev(in_x)] && std::abs(*std::prev(in_x) - i) > 1)
			is_continuation = true;

		if (std::next(in_y) != y_raster.points.end() && is_corner[*std::next(in_y)] && std::abs(*std::next(in_y) - i) > 1)
			is_continuation = true;
		if (in_y != y_raster.points.begin() && is_corner[*std::prev(in_y)] && std::abs(*std::prev(in_y) - i) > 1)
			is_continuation = true;		

		if(is_continuation)
#endif
			MatrixUtils::append(C, i);
	}
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)