#include <polyvec/utils/summable_to_one_sampler.hpp> // Optional

using namespace polyvec;


namespace polyvectest {
namespace util {

int
test_summable_to_one_sampler ( int /*argc*/, char** /*argv*/ ) {

    const polyvec::index n_variables = 4;
    const polyvec::index n_values = 5;
    SummableToOneSampler sampler(n_variables, n_values);

    printf("n_variables = %d \n", (int)n_variables);
    printf("n_values = %d \n", (int)n_values);
    
    sampler.begin();
    while(!sampler.has_finished()){
      std::vector<double> iterate = sampler.get_current_iteration();
      std::vector<polyvec::index> iterate_pos = sampler.get_positions_on_grid();
      printf("%05d: ", (int)sampler.iteration_no());
      for(polyvec::index i = 0 ; i < sampler.n_variables() ; ++i)
      {
        printf("%4.1f ", iterate[i]*(n_values-1));
      }
      printf( "|" );
      for(polyvec::index i = 0 ; i < sampler.n_variables() ; ++i)
      {
        printf("%4d ", (int)iterate_pos[i]);
      }
      printf("\n");
      sampler.step();
    }

    return 0;
}

} // globfitter
} // polyvectest