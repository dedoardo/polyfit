// Header
#include <polyvec/geom.hpp>

// Polyvec
#include <polyvec/misc.hpp>

// Wykobi
#include <wykobi/hdr/wykobi.hpp>
#include <wykobi/hdr/wykobi_math.hpp>
#include <wykobi/hdr/wykobi_algorithm.hpp>

// Eigen
#include <Eigen/Eigenvalues>

using namespace wykobi;
using namespace std;

namespace polyvec {

const Curve::Type Curve::types[3] = { Type::LINE, Type::ARC, Type::BEZIER };

//
// Clean up this mess of types below please...
//

// Internal types
using wray = wykobi::ray<double, 2>;
using wpoint = wykobi::point2d<double>;
using wsegment = wykobi::segment<double, 2>;
using wline = wykobi::line<double, 2>;

namespace geom {
// To/From backend (wykobi) types
namespace {
    polyvec_inline wray construct ( const real2& orig, const real2& dir ) {
        return make_ray<double> ( orig.x(), orig.y(), dir.x(), dir.y() );
    }

    polyvec_inline wpoint construct ( const real2& point ) {
        return make_point<double> ( point.x(), point.y() );
    }

    polyvec_inline wsegment construct ( const geom::segment& segment ) {
        return  make_segment <double> ( construct ( segment.src ), construct ( segment.dst ) );
    }

    polyvec_inline wline construct_axis ( const geom::segment& segment ) {
        return make_line<double> ( construct ( segment.src ), construct ( segment.dst ) );
    }

    polyvec_inline real2 destruct ( const wpoint& point ) {
        return real2 ( point.x, point.y );
    }
}


    real2 project ( const real2& pt, const segment& s ) {
        wpoint p = construct ( pt );
        wline axis = construct_axis ( s );
        wpoint p_proj = project_onto_axis ( p, axis ) [1];
        return destruct ( p_proj );
    }

    bool intersect ( const segment& lhs, const segment& rhs, real2* pt ) {
        wsegment a = construct ( lhs );
        wsegment b = construct ( rhs );

        if ( intersect ( a, b ) ) {
            if ( pt ) {
                wpoint p = intersection_point ( a, b );
                *pt = destruct ( p );
            }

            return true;
        }

        return false;
    }

    bool intersect_any ( const segment& l, const vector<segment>& others, real2* out ) {
        segment l_eps = l;
        l_eps.src += constants::SmallEps * l_eps.dir();
        l_eps.dst -= constants::SmallEps * l_eps.dir();

        for ( size_t i = 0; i < others.size(); ++i ) {
            if ( intersect ( l_eps, others[i], out ) ) {
                return true;
            }
        }

        return false;
    }

    bool ray_intersect ( const real2& orig, const real2& dir, const segment& l, double& ray_d ) {
        ray_d = 0.;
        wline ray = construct_axis ( geom::segment ( orig, orig + dir ) );
        wline segment = construct_axis ( l );

        if ( wykobi::intersect ( ray, segment ) ) {
            real2 p = destruct ( wykobi::intersection_point ( ray, segment ) );

            // Is the point along the segment ?
            if ( !point_on_line ( l, p ) ) {
                return false;
            }

            // ray distance
            real2 p_d = ( p - orig );
            ray_d = ( p_d ).norm();

            if ( abs ( p_d.x() ) > constants::SmallEps && misc::sign ( dir.x() ) != misc::sign ( p_d.x() ) ) {
                ray_d *= -1;
            } else if ( abs ( p_d.y() ) > constants::SmallEps && misc::sign ( dir.y() ) != misc::sign ( p_d.y() ) ) {
                ray_d *= -1;
            }

            return true;
        }

        return false;
    }

    double range_overlap ( double src0, double dst0, double src1, double dst1 ) {
        if ( src0 >= dst1 || dst0 <= src1 ) {
            return 0.;
        }

        return src0 < src1 ? ( dst0 - src1 ) : ( dst1 - src0 );
    }

    bool point_on_line ( const segment& segment, const real2& pt ) {
#if 0
        const double len = segment.arc_length() + constants::Eps;
        const real2 src2p = pt - segment.src;
        const real2 dst2p = pt - segment.dst;
        return src2p.norm() < len && dst2p.norm() > len;
#else
        return wykobi::point_on_segment ( construct ( pt ), construct ( segment ) );
#endif
    }

    double overlap ( const segment& lhs, const segment& rhs ) {
        real2 begp = project ( lhs.src, rhs );
        real2 endp = project ( lhs.dst, rhs );
        real2 dir = rhs.dst - rhs.src;

        begp = ( begp - rhs.src );
        endp = ( endp - rhs.src );

        if ( dir.x() > 0 ) {
            begp.x() /= dir.x();
            endp.x() /= dir.x();
        }

        if ( dir.y() > 0 ) {
            begp.y() /= dir.y();
            endp.y() /= dir.y();
        }

        double begs = begp.x() > 0 ? begp.norm() : -begp.norm();
        double ends = endp.x() > 0 ? endp.norm() : -endp.norm();

        if ( begs < ends ) {
            return misc::max0 ( ends ) - misc::max0 ( begs );
        }

        return misc::max0 ( begs ) - misc::max0 ( ends );
    }
}

}