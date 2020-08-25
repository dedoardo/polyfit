#ifdef POLYVEC_BUILD_PYTHON

#include <cstdio>

#include <polyvec/cycpputil.hpp>

namespace cycpputil_test
{
void
basic_types(FILE * fl)
{
  const double eps = 1e-10;

  CYCPPUTIL_DEFAULT_TRY;

    // Create an integer and make sure that it is equal to itself
    {
      const int value1 = 2341234;
      cycpputil::Py_unique_ptr py_value(cycpputil::int_to_pyobject(value1));
      const int value2 = cycpputil::pyobject_to_int(py_value.get());
      CYCPPUTIL_ASSERT(value1 == value2, "");
      if(fl) fprintf(fl, "%d == %d \n", value1, value2);
    }
    // Create a double and make sure that it is equal to itself
    {
      const double value1 = 2341234.2341829347;
      cycpputil::Py_unique_ptr py_value(cycpputil::double_to_pyobject(value1));
      const double value2 = cycpputil::pyobject_to_double(py_value.get());
      CYCPPUTIL_ASSERT(std::abs(value1 - value2) < eps, "");
      if(fl) fprintf(fl, "%lf == %lf \n", value1, value2);
    }
    // Create a string and make sure that it is equal to itself
    {
      const std::string value1 = "klsadjfklasjd89012374981jasdfkl;jka2903";
      cycpputil::Py_unique_ptr py_value(cycpputil::string_to_pyobject(value1));
      const std::string value2 = cycpputil::pyobject_to_string(py_value.get());
      CYCPPUTIL_ASSERT(value1 == value2, "");
      if(fl) fprintf(fl, "%s == %s \n", value1.c_str(), value2.c_str());
    }

  CYCPPUTIL_DEFAULT_CATCH();
}

} //  end of cycpputil_test
#endif