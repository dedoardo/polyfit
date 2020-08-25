#ifndef PIXVEC_VTK_CURVE_WIRTER_IS_INCLUDED
#define PIXVEC_VTK_CURVE_WIRTER_IS_INCLUDED

#include <Eigen/Core>
#include <string>
#include <vector>

namespace polyvec {

class VtkCurveWriter {
public:

    void add_line ( const Eigen::Vector2d&, const Eigen::Vector2d& );
    void add_polyline ( const Eigen::Matrix2Xd&, const bool is_close=false );
    void add_point ( const Eigen::Vector2d& );
    void dump ( const std::string );
    void clear();

private:
    template<typename DataType>
    using EigenVector = std::vector<DataType, Eigen::aligned_allocator<DataType>>;

    EigenVector<Eigen::Matrix2Xd> _lines;
    EigenVector<Eigen::Vector2d> _points;
};

} // end of polyvec

#endif /* PIXVEC_VTK_CURVE_WIRTER_IS_INCLUDED */
