#pragma once

#include <polyvec/api.hpp>

namespace polyvec {
	namespace util {

		//Calculates the normal direction for the given tangent direction. The normal
		//is the tangent rotated by 90° in clockwise direction. The magnitude of the
		//input vector is preserved.
		polyvec_inline Eigen::Vector2d normal_dir(const Eigen::Vector2d& vec) {
			return Eigen::Vector2d(-vec.y(), vec.x());
		}

	}
}