/*
	Utilities for Eigen::Matrix and std::vector for circular access and 
	insertion/removal. Use libc++ standard algorithms when applicable.
	
	The indexing routines are in the polyvec namespace, operations in a
	nested namespace
*/
#pragma once

// polyvec
#include <polyvec/core/types.hpp>
#include <polyvec/core/macros.hpp>

#define KEEP_ERASED_ENTRIES 0

NAMESPACE_BEGIN(polyfit)

inline void Reverse(mat2x& m) {
	for (int i = 0; i < m.cols() / 2; ++i) {
		vec2 c = m.col(i);
		m.col(i) = m.col(m.cols() - i - 1);
		m.col(m.cols() - i - 1) = c;
	}
}

template <typename Scalar, int Rows>
bool IsContained(const Eigen::Matrix<Scalar, Rows, 1>& v, const Scalar& val) {
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        if (v(i) == val) {
            return true;
        }
    }

    return false;
}

template <typename T>
T Circular(const T size, const T i) {
	return (i + size * 2) % size;
}

template <typename Scalar, int Rows, int Cols>
Eigen::Index Circular(const Eigen::Matrix<Scalar, Rows, Cols>& m, const Eigen::Index i) {
	return Circular(m.cols(), i);
}

template <typename Scalar, int Rows>
Eigen::Index Circular(const Eigen::Matrix<Scalar, Rows, 1>& v, const Eigen::Index i) {
	return Circular(v.size(), i);
}

template <typename Scalar, int Rows, int Cols>
Eigen::VectorXd CircularAt(const Eigen::Matrix<Scalar, Rows, Cols>& m, const Eigen::Index i) {
	return m.col(Circular(m, i));
}

template <typename Scalar, int Rows>
const Scalar& CircularAt(const Eigen::Matrix<Scalar, Rows, 1>& v, const Eigen::Index i) {
	return v(Circular(v, i), 0);
}

static Eigen::Index CircularDist(int size, int from, int to) {
	return (to - from) + (from > to ? size : 0);
}

template <typename Scalar, int Rows, int Cols>
Eigen::Index CircularDist(const Eigen::Matrix<Scalar, Rows, Cols>& m, const Eigen::Index from, const Eigen::Index to) {
	return CircularDist(m.cols(), from, to);
}

template <typename Scalar, int Rows>
Eigen::Index CircularDist(const Eigen::Matrix<Scalar, Rows, 1>& v, const Eigen::Index from, const Eigen::Index to) {
	return CircularDist(v.size(), from, to);
}

template <typename Scalar, int Rows>
int DirectionBetween(const Eigen::Matrix<Scalar, Rows, 1>& v, const Eigen::Index from, const Eigen::Index to) {
	if (from == to) {
		return 0;
	}

	return CircularDist(v, from, to) < CircularDist(v, to, from) ? +1 : -1;
}

template <typename Scalar, int Rows, int Cols>
Eigen::Index ShortestDist(const Eigen::Matrix<Scalar, Rows, Cols>& v, const Eigen::Index from, const Eigen::Index to) {
	return ::std::min(
		::std::abs(to - from),
		(to - from + v.cols()) % v.cols()
	);
}

template <typename Scalar, int Rows>
Eigen::Index ShortestDist(const Eigen::Matrix<Scalar, Rows, 1>& v, const Eigen::Index from, const Eigen::Index to) {
	return std::min(
		::std::abs(to - from),
		(to - from + v.size()) % v.size()
	);
}

template <typename T>
Eigen::Index ShortestDist(const std::vector<T>& v, const Eigen::Index from, const Eigen::Index to) {
	return std::min(
		::std::abs(to - from),
		(Eigen::Index)((to - from + v.size()) % v.size())
	);
}


template <typename T>
int Circular(const std::vector<T>& vec, const size_t index) {
	return Circular(vec.size(), index);
}

template <typename T>
T& CircularAt(std::vector<T> & vec, const size_t index) {
	return vec[Circular(vec, index)];
}

template <typename T>
const T& CircularAt(const std::vector<T>& vec, size_t index) {
	return vec[Circular(vec, index)];
}

template <typename T>
Eigen::Index CircularDist(const std::vector<T>& vec, const Eigen::Index from, const Eigen::Index to) {
	return (to - from) + (from > to ? vec.size() : 0);
}

template <typename T>
void EraseDuplicatesOrdered(std::vector<T>& vec) {
	auto it = vec.begin();
	while (it != vec.end()) {
		if ((it + 1) != vec.end() && *it == *(it + 1)) {
			it = vec.erase(it);
		} else {
			++it;
		}
	}
}

#if KEEP_ERASED_ENTRIES
template <typename T>
struct ErasedEntries
{
	static std::vector<T> entries;
};
template <typename T> std::vector<T> ErasedEntries<T>::entries = { };
#endif

template <typename T>
void EraseOrdered(std::vector<T>& vec, const std::vector<size_t>& dlist) {
	int off = 0;

	for (const size_t idx : dlist) {
		// assert maybe?
#if KEEP_ERASED_ENTRIES
		ErasedEntries<T>::entries.push_back(vec[idx - off]);
#endif
		vec.erase(vec.begin() + (idx - off++));
	}
}

// why Ordered?
template <typename Scalar, int Rows>
void EraseOrdered(Eigen::Matrix<Scalar, Rows, 1>& vec, const Eigen::Index i_remove) {
	PF_ASSERT(i_remove < vec.size());

	for (Eigen::Index i = i_remove; i < vec.size() - 1; ++i) {
		vec(i) = vec(i + 1);
	}

	vec.conservativeResize(vec.size() - 1);
}

template <typename Scalar, int Rows>
void EraseOrdered(Eigen::Matrix<Scalar, Rows, 1>& vec, const std::vector<size_t>& dlist) {
	size_t off = 0;
	for (const size_t idx : dlist) {
		EraseOrdered(vec, idx - off++);
	}
}

template <typename Scalar, int Rows>
Eigen::Index FindIndex(const Eigen::Matrix<Scalar, Rows, 1>& vec, const Scalar& val) {
	for (Eigen::Index i = 0; i < vec.size(); ++i) {
		if (vec(i) == val) {
			return i;
		}
	}

	return -1;
}

template <typename T>
Eigen::Index FindIndex(const std::vector<T>& vec, const T& val) {
	for (size_t i = 0; i < vec.size(); ++i) {
		if (vec[i] == val) {
			return (Eigen::Index)i;
		}
	}

	return -1;
}

template <typename Scalar, int Rows>
Eigen::Index InsertAt(Eigen::Matrix<Scalar, Rows, 1>& v, Eigen::Index at, const Scalar& val) {
	at = Circular(v, at);
	v.conservativeResize(v.size() + 1);
	for (Eigen::Index j = v.size() - 2; j >= at; --j) {
		v(j + 1) = v(j);
	}
	v(at) = val;
	return v.size();
}

namespace MatrixUtils {
	inline void insert_at(Eigen::Matrix2Xd& points, Eigen::Index i, const Eigen::Vector2d& val) {
		i = Circular(points, i);
		points.conservativeResize(2, points.cols() + 1);
		for (Eigen::Index j = points.cols() - 2; j >= i; --j) {
			points.col(j + 1) = points.col(j);
		}
		points.col(i) = val;
	}

    inline void insert_at(Eigen::VectorXi& points, Eigen::Index i, const int& val) {
        i = Circular(points, i);
        points.conservativeResize(points.size() + 1);
        for (Eigen::Index j = points.size() - 2; j >= i; --j) {
            points(j + 1) = points(j);
        }
        points(i) = val;
    }

	inline void append(Eigen::Matrix2Xd& m, const Eigen::Vector2d& col) {
		m.conservativeResize(Eigen::NoChange, m.cols() + 1);
		m.col(m.cols() - 1) = col;
	}

	inline void append(Eigen::Matrix2Xi& m, const Eigen::Vector2i& col) {
		m.conservativeResize(Eigen::NoChange, m.cols() + 1);
		m.col(m.cols() - 1) = col;
	}

	inline void append(Eigen::Matrix3Xi& m, const Eigen::Vector3i& col) {
		m.conservativeResize(Eigen::NoChange, m.cols() + 1);
		m.col(m.cols() - 1) = col;
	}

	inline void append(Eigen::Matrix4Xi& m, const Eigen::Vector4i& col) {
		m.conservativeResize(Eigen::NoChange, m.cols() + 1);
		m.col(m.cols() - 1) = col;
	}

	inline void append(Eigen::VectorXi& v, const int c) {
		v.conservativeResize(v.size() + 1);
		v(v.size() - 1) = c;
	}

	template <typename T>
	inline void insert_at(std::vector<T>& vec, Eigen::Index i, const T& val) {
		i = Circular(vec, i);
		vec.insert(vec.begin() + i, val);
	}

	inline void remove_at(Eigen::Matrix2Xd& points, Eigen::Index i) {
		i = Circular(points, i);
		for (Eigen::Index j = i + 1; j < points.cols(); ++j)
		{
			points.col(j - 1) = points.col(j);
		}
		points.conservativeResize(2, points.cols() - 1);
	}
}

NAMESPACE_END(polyfit)