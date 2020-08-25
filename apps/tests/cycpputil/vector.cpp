#ifdef POLYVEC_BUILD_PYTHON

#include <cstdio>

#include <polyvec/cycpputil.hpp>


namespace cycpputil_test
{

void
vector(FILE * fl)
{
  const double eps = 1e-10;

  CYCPPUTIL_DEFAULT_TRY;

    // Create an integer and make sure that it is equal to itself
    {
      const std::vector<int> value1 = {2341234, 2341, 1234, 624562456};
      cycpputil::Py_unique_ptr py_value(cycpputil::vectori_to_pyobject(value1));
      const std::vector<int> value2 = cycpputil::pyobject_to_vectori(py_value.get());
      CYCPPUTIL_ASSERT(value1 == value2, "");
      if(fl)
      {
        for(unsigned i = 0; i < value1.size(); ++i)
          fprintf(fl, "%d == %d, ", value1[i], value2[i]);
        fprintf(fl, "\n");
      }
    }
    // Create a double and make sure that it is equal to itself
    {
      const std::vector<double> value1 = {554.123, 2314.564, 54454.234, 1234.54};
      cycpputil::Py_unique_ptr py_value(cycpputil::vectord_to_pyobject(value1));
      const std::vector<double> value2 = cycpputil::pyobject_to_vectord(py_value.get());
      for(unsigned i = 0; i < value1.size(); ++i)
        CYCPPUTIL_ASSERT(std::abs(value1[i] - value2[i]) < eps, "");
      if(fl)
      {
        for(unsigned i = 0; i < value1.size(); ++i)
          fprintf(fl, "%lf == %lf, ", value1[i], value2[i]);
        fprintf(fl, "\n");
      };
    }
    // Create a string and make sure that it is equal to itself
    {
      const std::vector<std::string> value1 = {"wqeqwe", "xxx", "ZdxQ", "123dX"};
      cycpputil::Py_unique_ptr py_value(cycpputil::vectorstr_to_pyobject(value1));
      const std::vector<std::string> value2 = cycpputil::pyobject_to_vectorstr(py_value.get());
      CYCPPUTIL_ASSERT(value1 == value2, "");
      if(fl)
      {
        for(unsigned i = 0; i < value1.size(); ++i)
          fprintf(fl, "%s == %s, ", value1[i].c_str(), value2[i].c_str());
        fprintf(fl, "\n");
      }
    }
    CYCPPUTIL_DEFAULT_CATCH();


} // end of vector()
} // end of cycpputil_test

#endif