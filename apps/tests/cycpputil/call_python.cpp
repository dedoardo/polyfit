#ifdef POLYVEC_BUILD_PYTHON

#include <cstdio>

#include <polyvec/cycpputil.hpp>


namespace cycpputil_test
{

void
call_python_basic(FILE * fl)
{

  const double eps = 1e-10;
  CYCPPUTIL_UNUSED(eps);

      CYCPPUTIL_DEFAULT_TRY;

    fprintf(fl, "BASIC MODULE \n");

    // Load the module
    cycpputil::Py_unique_ptr py_module(cycpputil::pycall_load_module("basic_module"));
#define LOAD_MODULE(mname) \
  cycpputil::Py_unique_ptr py_##mname(cycpputil::pycall_load_function(py_module.get(), #mname));
    LOAD_MODULE(sum_two_numbers_as_str);
#undef LOAD_MODULE

    const int n0 = 12;
    const int n1 = 20;
    cycpputil::Py_unique_ptr py_kwargs(PyDict_New());
    cycpputil::pydict_absorb_item(py_kwargs.get(), "n0", cycpputil::double_to_pyobject(n0));
    cycpputil::pydict_absorb_item(py_kwargs.get(), "n1", cycpputil::double_to_pyobject(n1));
    cycpputil::Py_unique_ptr py_ans(
        cycpputil::pycall_call_function(py_sum_two_numbers_as_str.get(), NULL, py_kwargs.get()));
    const std::string ans = cycpputil::pyobject_to_string(py_ans.get());
    int int_ans;
    sscanf(ans.c_str(), "%d", &int_ans);
    CYCPPUTIL_ASSERT(n0 + n1 == int_ans, n0 + n1 << "!=" << int_ans);
    fprintf(fl, " Success in summing: %d+%d == %d \n", n0, n1, int_ans);

     CYCPPUTIL_DEFAULT_CATCH();

}

void
call_python_advanced(FILE * fl)
{

  const double eps = 1e-10;
  CYCPPUTIL_UNUSED(eps);

      CYCPPUTIL_DEFAULT_TRY;

    fprintf(fl, "ADVANCED MODULE \n");

    // Load the module
    cycpputil::Py_unique_ptr py_module(cycpputil::pycall_load_module("advanced_module"));
#define LOAD_MODULE(mname) \
  cycpputil::Py_unique_ptr py_##mname(cycpputil::pycall_load_function(py_module.get(), #mname));
    LOAD_MODULE(flatten_array);
    LOAD_MODULE(generate_array);
    LOAD_MODULE(fprint_object);
    LOAD_MODULE(print_object);
    LOAD_MODULE(pickle_object);
    LOAD_MODULE(unpickle_object);
#undef LOAD_MODULE

    // Create the matrix [1,2,3 \\ 4,5,6] in numpy
    std::vector<double> ref_data = {1, 2, 3, 4, 5, 6};
    std::vector<int> ref_dims = {2, 3};

     // Print the data
    fprintf(fl, "ref_data: ");
    for(auto p : ref_data)
      fprintf(fl, "%g ", p);
    fprintf(fl, "\n");
    fprintf(fl,"ref_dims: ");
    for(auto d : ref_dims)
      fprintf(fl, "%d ", d);
    fprintf(fl, "\n");

    // Convert to numpy matrix
    cycpputil::Py_unique_ptr py_numpy_matrix;
    {
      cycpputil::Py_unique_ptr py_data(cycpputil::vectord_to_pyobject(ref_data));
      cycpputil::Py_unique_ptr py_dims(cycpputil::vectori_to_pyobject(ref_dims));
      cycpputil::Py_unique_ptr py_kwargs(PyDict_New());
      cycpputil::pydict_absorb_item(py_kwargs.get(), "dims", py_dims.get());
      cycpputil::pydict_absorb_item(py_kwargs.get(), "data", py_data.get());
      py_numpy_matrix.reset(cycpputil::pycall_call_function(py_generate_array.get(), nullptr, py_kwargs.get()));
    }

    // Print the numpy array
    fprintf(fl, "Python is going to print the matrix (stdout only)\n");
    {
      cycpputil::Py_unique_ptr py_kwargs(PyDict_New());
      cycpputil::pydict_set_item(py_kwargs.get(), "obj", py_numpy_matrix.get());
      cycpputil::pycall_call_function(py_print_object.get(), nullptr, py_kwargs.get());
    }

    // Now pickle the numpy matrix
    fprintf(fl, "pickling to _dummy.pickle\n");
    {
      cycpputil::Py_unique_ptr py_kwargs(PyDict_New());
      cycpputil::pydict_set_item(py_kwargs.get(), "obj", py_numpy_matrix.get());
      cycpputil::pydict_absorb_item(py_kwargs.get(), "fname", cycpputil::string_to_pyobject("_dummy.pickle"));
      cycpputil::pycall_call_function(py_pickle_object.get(), nullptr, py_kwargs.get());
    }

    // Now unpickle it
    fprintf(fl, "unpickling from _dummy.pickle\n");
    cycpputil::Py_unique_ptr py_unpickled_matrix;
    {
      cycpputil::Py_unique_ptr py_kwargs(PyDict_New());
      cycpputil::pydict_absorb_item(py_kwargs.get(), "fname", cycpputil::string_to_pyobject("_dummy.pickle"));
      py_unpickled_matrix.reset(cycpputil::pycall_call_function(py_unpickle_object.get(), nullptr, py_kwargs.get()));
    }

    // Print the unpickled matrix
    fprintf(fl, "Python is going to print the unpickled object (stdout only)\n");
    {
      cycpputil::Py_unique_ptr py_kwargs(PyDict_New());
      cycpputil::pydict_set_item(py_kwargs.get(), "obj", py_unpickled_matrix.get());
      cycpputil::pycall_call_function(py_print_object.get(), nullptr, py_kwargs.get());
    }

    // Convert back to array
    fprintf(fl, "Converting back array to c++ \n");
    std::vector<double> data;
    std::vector<int> dims;
    {
      cycpputil::Py_unique_ptr py_ans_tuple;
      cycpputil::Py_unique_ptr py_kwargs(PyDict_New());
      cycpputil::pydict_set_item(py_kwargs.get(), "array", py_unpickled_matrix.get());
      py_ans_tuple.reset(cycpputil::pycall_call_function(py_flatten_array.get(), nullptr, py_kwargs.get()));

      // borrow ref
      PyObject * py_dims = PyTuple_GetItem(py_ans_tuple.get(), 0);
      PyObject * py_data = PyTuple_GetItem(py_ans_tuple.get(), 1);
      data = cycpputil::pyobject_to_vectord(py_data);
      dims = cycpputil::pyobject_to_vectori(py_dims);
    }

    // Print the data
    fprintf(fl, "data: ");
    for(auto p : data)
      fprintf(fl, "%g ", p);
    fprintf(fl, "\n");
    fprintf(fl,"dims: ");
    for(auto d : dims)
      fprintf(fl, "%d ", d);
    fprintf(fl, "\n");

    // make sure it is the same as the initial ones
    CYCPPUTIL_ASSERT( ref_data == data, "data does not match");
    CYCPPUTIL_ASSERT( ref_dims == dims, "dims don't match" );
    
      CYCPPUTIL_DEFAULT_CATCH();

}
}

#endif