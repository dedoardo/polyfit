#ifndef POLYFIT_TYPEDEFS_H_
#define POLYFIT_TYPEDEFS_H_

// Eigen
#include <Eigen/Core>

// libc++
#include <cstdint>
#include <utility>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

// todo remove
#define PF_DEBUG
#ifdef _DEBUG 
#define PF_DEBUG
#endif

namespace polyfit {
	using vec2 = Eigen::Vector2d;
	using vec3 = Eigen::Vector3d;
	using vec4 = Eigen::Vector4d;

	using vec2i = Eigen::Vector2i;
	using vec3i = Eigen::Vector3i;
	using vec4i = Eigen::Vector4i;

	using mat2 = Eigen::Matrix<double, 2, 2>;
	using mat3 = Eigen::Matrix<double, 3, 3>;
	using mat4 = Eigen::Matrix<double, 4, 4>;

	using mat2x = Eigen::Matrix2Xd;
	using mat2xi = Eigen::Matrix2Xi;
	using mat3xi = Eigen::Matrix3Xi;
	using mat4xi = Eigen::Matrix4Xi;
	using mat4i = Eigen::Matrix4i;
	using mat23 = Eigen::Matrix<double, 2, 3>;
	using mat24 = Eigen::Matrix<double, 2, 4>;

	using mat24i = Eigen::Matrix<int, 2, 4>;

	using vecXd = Eigen::VectorXd;
	using vecXi = Eigen::VectorXi;

	using Vertex = Eigen::Index;

	using matXd = Eigen::MatrixXd;
	using matXi = Eigen::MatrixXi;

	// https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
	inline uint64_t hash_combine(uint64_t lhs, uint64_t rhs) {
		lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
		return lhs;
	}

	struct _StdHashVertexPair {
		uint64_t operator() (const std::pair<Vertex, Vertex>& t) const {
			return hash_combine((uint64_t) t.first, (uint64_t)t.second);
		}
	};
}

#endif // POLYFIT_TYPEDEFS_H_