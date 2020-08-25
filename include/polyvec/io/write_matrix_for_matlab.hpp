#ifndef polyvec_write_matrix_for_matlab_
#define polyvec_write_matrix_for_matlab_

#include <ostream>
#include <vector>
#include <Eigen/Core>

namespace polyvec {

    template<typename T>
    void
    _write_vector_for_matlab_impl ( std::ostream& output,
                                    T& a,
                                    const char* variable_name,
                                    bool column_vector = true,
                                    int significant_digits = 18 ) {
        output << variable_name << "=[";
        std::streamsize old_precision = output.precision();
        output.precision ( significant_digits );

        for ( unsigned int i = 0; i < a.size(); ++i ) {
            output << a[i] << " ";
        }

        output << "]";

        if ( column_vector ) {
            output << "'";
        }

        output << ";" << std::endl;
        output.precision ( old_precision );
    }

    inline void
    write_vector_for_matlab ( std::ostream& output,
                              const Eigen::VectorXd& a,
                              const char* variable_name,
                              bool column_vector = true,
                              int significant_digits = 18 ) {
        _write_vector_for_matlab_impl ( output, a, variable_name, column_vector, significant_digits );
    }

    inline void
    write_vector_for_matlab ( std::ostream& output,
                              const std::vector<double>& a,
                              const char* variable_name,
                              bool column_vector = true,
                              int significant_digits = 18 ) {
        _write_vector_for_matlab_impl ( output, a, variable_name, column_vector, significant_digits );
    }


    inline void
    write_full_matrix_for_matlab ( std::ostream& output,
                                   const Eigen::MatrixXd& a,
                                   const char* variable_name,
                                   int significant_digits = 18 ) {
        output << variable_name << "=[";
        std::streamsize old_precision = output.precision();
        output.precision ( significant_digits );

        for ( unsigned int i = 0; i < a.rows(); ++i ) {
            for ( unsigned int j = 0; j < a.cols(); ++j ) {
                output << a ( i, j ) << " ";
            }

            if ( i != a.rows() - 1 ) {
                output << "; ... \n";
            }
        }

        output << "]";
        output << ";" << std::endl;
        output.precision ( old_precision );
    }

}

#endif // polyvec_polygon_h_
