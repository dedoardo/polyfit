// Sample a bunch of values from the [0 1] such that their sum is 1
// Author: Shayan Hoshyari

#include <polyvec/api.hpp>

namespace polyvec {
    class SummableToOneSampler {
    public:
        // n_variables_in: number of variables
        // n_grid_points: number of values that each variable can attain in 0 1
        SummableToOneSampler ( const index n_variables_in, const index n_grid_points_in ) ;

        // Start iterating
        void begin();

        // Step the iteration
        void step();

        // Have we called begin()?
        bool has_began();

        // Have we reached the last iterate
        bool has_finished();

        // Get the values for the current iteration
        std::vector<double> get_current_iteration();


        // Querry
        index n_variables();
        index n_grid_points();
        index delta();
        index iteration_no();
        std::vector<index> get_positions_on_grid();

    private:
        // General info
        index _n_variables;
        index _n_grid_points;
        double _delta;
        bool _has_finished;

        // Current iteration
        index _iteration_no;
        std::vector<index> _positions_on_grid;

        void check_validity();
        // bool has_finished_helper();
    };
}