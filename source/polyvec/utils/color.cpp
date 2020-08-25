#include <polyvec/utils/color.hpp>
#include <polyvec/utils/num.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Color)

const vec3 COLOR_MIN(22. / 255, 64. / 255, 229. / 255);
const vec3 COLOR_MID(30. / 255, 255. / 255, 97. / 255);
const vec3 COLOR_MAX(255. / 255, 45. / 255, 30. / 255);

// Returns some kind of interpolated color
// <= t_min -> blue
// .5 (t_min + t_max) -> green
// > t_max -> red
// requires t_min < t_max
vec3 error(double t_min, double t_max, double t) {
	const double t_mid = polyvec::Num::lerp(t_min, t_max, .5);
	if (t <= t_mid) {
		return polyvec::Num::lerp(COLOR_MIN, COLOR_MID, t / t_mid);
	} else {
		return polyvec::Num::lerp(COLOR_MID, COLOR_MAX, polyvec::Num::saturate((t - t_mid) / t_mid));
	}
}

NAMESPACE_END(Color)
NAMESPACE_END(polyfit)