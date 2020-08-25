#pragma once

#include <ostream>
#include <Eigen/Core>

namespace nse
{
	namespace util
	{

		//Wraps an Eigen matrix. If this object is printed with an operator<<, the it will be formatted in Mathematica's syntax.
		template <typename EigenMatrix>
		class MathematicaFormatter
		{
		public:
			MathematicaFormatter(const EigenMatrix& matrix)
				: matrix(matrix)
			{ }

		private:
			const EigenMatrix& matrix;

			friend std::ostream& (operator<<)(std::ostream& stream, const MathematicaFormatter<EigenMatrix>& fmt)
			{
				stream << "{ ";
				for (int i = 0; i < fmt.matrix.rows(); ++i)
				{
					stream << "{ ";
					for (int j = 0; j < fmt.matrix.cols(); ++j)
					{
						stream << fmt.matrix.coeff(i, j);
						if (j != fmt.matrix.cols() - 1)
							stream << ", ";
					}
					stream << " }";
					if (i != fmt.matrix.rows() - 1)
						stream << ", ";
				}
				stream << " }";
				return stream;
			}
		};

		//Specialization of the MathematicaFormatter for sparse matrices
		template <typename Scalar, int Options>
		class MathematicaFormatter<Eigen::SparseMatrix<Scalar, Options>>
		{
		public:
			MathematicaFormatter(const Eigen::SparseMatrix<Scalar, Options>& matrix)
				: matrix(matrix)
			{ }

		private:
			const Eigen::SparseMatrix<Scalar, Options>& matrix;

			friend std::ostream& (operator<<)(std::ostream& stream, const MathematicaFormatter<Eigen::SparseMatrix<Scalar, Options>>& fmt)
			{
				bool first = true;
				stream << "SparseArray[{ ";
				for (int i = 0; i < fmt.matrix.rows(); ++i)
				{
					for (typename Eigen::SparseMatrix<Scalar, Options>::InnerIterator it(fmt.matrix, i); it; ++it)
					{
						if (!first)
						{
							stream << ", ";
						}
						first = false;
						int row, col;
						if (Options & Eigen::RowMajor)
						{
							row = i;
							col = it.index();
						}
						else
						{
							row = it.index();
							col = i;
						}
						stream << "{" << (row + 1) << ", " << (col + 1) << "}->" << it.value();
					}
				}
				stream << " }, {" << fmt.matrix.rows() << ", " << fmt.matrix.cols() << "}]";
				return stream;
			}
		};

		template <typename EigenMatrix>
		MathematicaFormatter<EigenMatrix> FormatMathematica(const EigenMatrix& matrix)
		{
			return MathematicaFormatter<EigenMatrix>(matrix);
		}
	}
}