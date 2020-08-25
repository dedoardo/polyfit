#include <polyvec/utils/summable_to_one_sampler.hpp>


namespace polyvec {

    SummableToOneSampler::SummableToOneSampler ( const index n_variables_in, const index n_grid_points_in ) {
        _n_variables = n_variables_in;
        _n_grid_points = n_grid_points_in;
        assert_break ( n_variables_in >= 2 );
        assert_break ( n_grid_points_in >= 3 );
        _delta = 1. / ( n_grid_points_in-1 );
        _iteration_no = -1;
        _positions_on_grid = std::vector<index> ( n_variables_in, -1 );
        _has_finished = false;
    }

    void SummableToOneSampler::begin() {
        std::fill ( _positions_on_grid.begin(), _positions_on_grid.end(), 0 );
        _positions_on_grid.back() = _n_grid_points-1;
        _has_finished = false;
        _iteration_no = 0;
    }

    void SummableToOneSampler::step() {
        if ( has_began() && ( !has_finished() ) ) {

            auto increment_jth_variable =[this] ( const int j ) ->void {
                assert_break ( j != this->_n_variables-1 );
                assert_break ( this->_positions_on_grid[j] != this->_n_grid_points - 1 );
                ++ this->_positions_on_grid[j];

                for ( int i = this->_n_variables-2 ; i > j ; --i ) {
                    assert_break ( this->_positions_on_grid[i] == this->_n_grid_points - 1 );
                    this->_positions_on_grid[i] = this->_positions_on_grid[j];
                }
            };

            bool could_increment = false;

            for ( int i = this->_n_variables-2 ; i >= 0 ; --i ) {
                if ( _positions_on_grid[i] != _n_grid_points -1 ) {
                    increment_jth_variable ( i );
                    could_increment = true;
                    break;
                }
            }

            if ( !could_increment ) {
                _has_finished = true;
            } else {
                ++_iteration_no;
            }
        } // end of has began and not finished
    }

    bool SummableToOneSampler::has_finished() {
        return _has_finished;
    }

    bool SummableToOneSampler::has_began() {
        return _iteration_no >= 0;
    }

    std::vector<double> SummableToOneSampler::get_current_iteration() {
        check_validity();

        if ( has_began() && ( !has_finished() ) ) {
            std::vector<double> ans ( _n_variables );

            ans[0] = _positions_on_grid[0] * _delta;

            for ( index i = 1 ; i < _n_variables ; ++i ) {
                ans[i] = ( _positions_on_grid[i] - _positions_on_grid[i-1] ) * _delta;
            }

            return ans;
        } else {
            return std::vector<double>();
        }
    }

    index SummableToOneSampler::n_variables() {
        return _n_variables;
    }

    index SummableToOneSampler::n_grid_points() {
        return _n_grid_points;
    }

    index SummableToOneSampler::delta() {
        return _delta;
    }

    void SummableToOneSampler::check_validity() {
        if ( has_began() && ( !has_finished() ) ) {
            for ( index i = 0 ; i < _n_variables ; ++i ) {
                assert_break ( _positions_on_grid[i] >= 0 );
                assert_break ( ( i==0 ) || ( _positions_on_grid[i] >=  _positions_on_grid[i-1] ) );
            }

            assert_break ( _positions_on_grid.back() ==  _n_grid_points-1 );
        }
    }

    index SummableToOneSampler::iteration_no() {
        return _iteration_no;
    }

    std::vector<index> SummableToOneSampler::get_positions_on_grid() {
        return _positions_on_grid;
    }

    // bool SummableToOneSampler::has_finished_helper() {

    //     if ( !has_began() ) {
    //         return false;
    //     }

    //     for ( index i = 0 ; i < _n_variables ; ++i ) {
    //         if ( _positions_on_grid[i] != _n_grid_points - 1 ) {
    //             return false;
    //         }
    //     }

    //     return true;
    // }


} // end of polyvec