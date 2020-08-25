#include <Eigen/Core>

template <typename T, int R, int C>
struct ConstSerializationWrapper
{
	const Eigen::Matrix<T, R, C>& m;

	ConstSerializationWrapper(const Eigen::Matrix<T, R, C>& m)
		: m(m)
	{ }
};

template <typename T, int R, int C>
struct SerializationWrapper
{
	Eigen::Matrix<T, R, C>& m;

	SerializationWrapper(Eigen::Matrix<T, R, C>& m)
		: m(m)
	{ }	

	operator ConstSerializationWrapper<T, R, C>() const
	{
		return ConstSerializationWrapper<T, R, C>(m);
	}
};


template <typename T, int R, int C>
SerializationWrapper<T, R, C> Serialized(Eigen::Matrix<T, R, C>& m)
{
	return SerializationWrapper<T, R, C>(m);
}

template <typename T, int R, int C>
ConstSerializationWrapper<T, R, C> Serialized(const Eigen::Matrix<T, R, C>& m)
{
	return ConstSerializationWrapper<T, R, C>(m);
}

template <typename T, int R, int C>
std::ostream& operator<<(std::ostream& s, const ConstSerializationWrapper<T, R, C>& m)
{
	s << m.m.rows() << " " << m.m.cols();
	for (int i = 0; i < m.m.rows(); ++i)
		for (int j = 0; j < m.m.cols(); ++j)
			s << " " << m.m.coeff(i, j);
	return s;
}

template <typename T, int R, int C>
std::ostream& operator<<(std::ostream& s, const SerializationWrapper<T, R, C>& m)
{
	return s << static_cast<const ConstSerializationWrapper<T, R, C>>(m);
}

template <typename T, int R, int C>
std::istream& operator>>(std::istream& s, const SerializationWrapper<T, R, C>& m)
{
	int rows, cols;
	s >> rows >> cols;
	m.m.resize(rows, cols);
	for (int i = 0; i < m.m.rows(); ++i)
		for (int j = 0; j < m.m.cols(); ++j)
			s >> m.m.coeffRef(i, j);
	return s;
}