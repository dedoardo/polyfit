#include <limits>

#ifdef POLYVEC_BUILD_PYTHON

#include <polyvec/cycpputil.hpp>


namespace cycpputil
{


// ============================================================
//                              UTIL
// ============================================================

namespace
{


template<typename Data_type>
std::vector<Data_type>
pyobject_to_vector_impl(PyObject * py_array, Data_type (*convert_function)(PyObject * obj))
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyList_Check(py_array), "Object is not an array");
  const unsigned n_items = PyList_Size(py_array);
  std::vector<Data_type> ans(n_items);
  for(unsigned i = 0; i < n_items; ++i)
  {
    PyObject * py_val = PyList_GetItem(py_array, (Py_ssize_t)i); // borrow ref
    ans[i] = convert_function(py_val);
  }
  return ans;

  CYCPPUTIL_DEFAULT_CATCH(std::vector<Data_type>());
}

template<typename Data_type>
PyObject *
vector_to_pyobject_impl(const std::vector<Data_type> & array, PyObject * (*convert_function)(const Data_type obj))
{
  CYCPPUTIL_DEFAULT_TRY;

  PyObject * py_array = PyList_New((Py_ssize_t)array.size());
  for(unsigned i = 0; i < array.size(); ++i)
  {
    PyObject * entry = convert_function(array[i]);
    PyList_SET_ITEM(py_array, (Py_ssize_t)i, entry); // steals the reference
  }
  return py_array;

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

} // end of anonymus namespace

// ============================================================
//                              ERROR
// ============================================================

void
Error::report(std::ostream & ostr)
{
  const char * c_red = "\033[1;31m";
  const char * c_reset = "\033[0m";

  ostr << c_red << "================= ERR BEGIN =================" << c_reset << "\n";
  ostr << this->msg;
  for(unsigned i = 0; i < this->funcs.size(); ++i)
  {
    ostr << i << this->funcs[i] << "()\n";
  }
  ostr << c_red << "================= ERR END   =================" << c_reset << "\n";
}

void
Error::report(FILE * fl)
{
  const char * c_red = "\033[1;31m";
  const char * c_reset = "\033[0m";

  fprintf(fl, " %s ================= ERR BEGIN ================= %s \n", c_red, c_reset);
  fprintf(fl, "%s \n", this->msg.c_str());
  for(unsigned i = 0; i < this->funcs.size(); ++i)
  {
    fprintf(fl, "%s() \n", this->funcs[i].c_str());
  }
  fprintf(fl, " %s ================= ERR END ================= %s \n", c_red, c_reset);
}


// ============================================================
//                              Py_unique_ptr
// ============================================================


Py_unique_ptr::Py_unique_ptr(PyObject * val)
{
  _val = val;
}

void
Py_unique_ptr::reset(PyObject * val)
{
  if(_val != nullptr)
    Py_DECREF(_val);
  _val = val;
}

Py_unique_ptr::~Py_unique_ptr()
{
  if(_val != nullptr)
    Py_DECREF(_val);
}

PyObject *
Py_unique_ptr::get()
{
  return _val;
}

// ============================================================
//                              Py_init
// ============================================================

Py_init::Py_init()
{
  ::Py_Initialize();
}

Py_init::~Py_init()
{
  ::Py_Finalize();
}

// ============================================================
//                              VECTOR
// ============================================================

std::vector<int>
pyobject_to_vectori(PyObject * py_array)
{
  CYCPPUTIL_DEFAULT_TRY;

  return pyobject_to_vector_impl<int>(py_array, pyobject_to_int);

  CYCPPUTIL_DEFAULT_CATCH(std::vector<int>());
}

PyObject *
vectori_to_pyobject(const std::vector<int> & array)
{
  CYCPPUTIL_DEFAULT_TRY;

  return vector_to_pyobject_impl<int>(array, int_to_pyobject);

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

std::vector<double>
pyobject_to_vectord(PyObject * py_array)
{
  CYCPPUTIL_DEFAULT_TRY;

  return pyobject_to_vector_impl<double>(py_array, pyobject_to_double);

  CYCPPUTIL_DEFAULT_CATCH(std::vector<double>());
}

PyObject *
vectord_to_pyobject(const std::vector<double> & array)
{
  CYCPPUTIL_DEFAULT_TRY;

  return vector_to_pyobject_impl<double>(array, double_to_pyobject);

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

std::vector<std::string>
pyobject_to_vectorstr(PyObject * py_array)
{
  CYCPPUTIL_DEFAULT_TRY;

  return pyobject_to_vector_impl<std::string>(py_array, pyobject_to_string);

  CYCPPUTIL_DEFAULT_CATCH(std::vector<std::string>());
}

PyObject *
vectorstr_to_pyobject(const std::vector<std::string> & array)
{
  CYCPPUTIL_DEFAULT_TRY;

  return vector_to_pyobject_impl<std::string>(array, string_to_pyobject);

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

// ============================================================
//                              BASIC TYPES
// ============================================================

int
pyobject_to_int(PyObject * obj)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyLong_Check(obj), "obj is not of type int (python long)");
  return (int)PyLong_AsLong(obj);

  CYCPPUTIL_DEFAULT_CATCH(std::numeric_limits<int>::max());
}


PyObject *
int_to_pyobject(const int inp)
{
  CYCPPUTIL_DEFAULT_TRY;

  return PyLong_FromLong((long)inp);

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}


double
pyobject_to_double(PyObject * obj)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyFloat_Check(obj), "obj is not of type float ");
  return PyFloat_AsDouble(obj);

  CYCPPUTIL_DEFAULT_CATCH(std::numeric_limits<double>::max());
}

PyObject *
double_to_pyobject(const double inp)
{
  CYCPPUTIL_DEFAULT_TRY;

  return PyFloat_FromDouble(inp);

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}


std::string
pyobject_to_string(PyObject * obj)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyUnicode_Check(obj), "Object is not string");
  return std::string(PyUnicode_AsUTF8AndSize(obj, NULL));

  CYCPPUTIL_DEFAULT_CATCH("");
}

PyObject *
string_to_pyobject(const std::string inp)
{
  CYCPPUTIL_DEFAULT_TRY;
  return PyUnicode_FromString(inp.c_str());

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

// ============================================================
//                              DICTIONARY
// ============================================================


bool
pydict_has_item(PyObject * py_dict, const std::string & key)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyDict_Check(py_dict), "Object is not python DICT");
  Py_unique_ptr py_key(string_to_pyobject(key));
  return (bool)PyDict_Contains(py_dict, py_key.get());

  CYCPPUTIL_DEFAULT_CATCH(false);
}

PyObject *
pydict_get_item(PyObject * py_dict, const std::string & key)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(pydict_has_item(py_dict, key), "Object does not contain key: " << key);
  return PyDict_GetItemString(py_dict, key.c_str());

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

bool
pydict_set_item(PyObject * py_dict, const std::string & key, PyObject * item)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyDict_Check(py_dict), "Object is not dict. ");
  return PyDict_SetItemString(py_dict, key.c_str(), item); // no steal -- bump ref count

  CYCPPUTIL_DEFAULT_CATCH(false);
}

bool
pydict_absorb_item(PyObject * py_dict, const std::string & key, PyObject * item)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyDict_Check(py_dict), "Object is not dict. ");
  Py_unique_ptr py_item(item); // absorb
  return PyDict_SetItemString(py_dict, key.c_str(), py_item.get()); // no steal

  CYCPPUTIL_DEFAULT_CATCH(false);
}

// ============================================================
//                         CALLING PYTHON
// ============================================================


PyObject *
pycall_load_module(const std::string module_name)
{
  CYCPPUTIL_DEFAULT_TRY;

  Py_unique_ptr py_name(string_to_pyobject(module_name));
  PyObject * py_module = PyImport_Import(py_name.get());
  // PyObject * py_module = PyImport_ImportModule(module_name.c_str());
  if(PyErr_Occurred())
  {
    PyErr_Print();
    CYCPPUTIL_ASSERT(0, "Error occured inside module " << module_name);
  }
  CYCPPUTIL_ASSERT(py_module != NULL, "Module ``" << module_name << "'' was not found");
  return py_module;

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

PyObject *
pycall_load_function(PyObject * py_module, const std::string function_name)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyModule_Check(py_module), "Object is not a module");
  PyObject * py_func = PyObject_GetAttrString(py_module, function_name.c_str());
  if(PyErr_Occurred())
  {
    PyErr_Print();
    CYCPPUTIL_ASSERT(0, "Error occured ");
  }
  CYCPPUTIL_ASSERT(PyCallable_Check(py_func), function_name << " is not a callable attribute");
  return py_func;

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

PyObject *
pycall_call_function(PyObject * py_function, PyObject * py_args, PyObject * py_kwargs)
{
  CYCPPUTIL_DEFAULT_TRY;

  CYCPPUTIL_ASSERT(PyCallable_Check(py_function), "Object is not callable");
  PyObject * ans;
  if(py_args == nullptr)
  {
    Py_unique_ptr dummy_args(PyTuple_New(0));
    ans = PyObject_Call(py_function, dummy_args.get(), py_kwargs);
  }
  else
  {
    ans = PyObject_Call(py_function, py_args, py_kwargs);
  }
  if(PyErr_Occurred())
  {
    PyErr_Print();
    CYCPPUTIL_ASSERT(0, "Error occured ");
  }
  return ans;

  CYCPPUTIL_DEFAULT_CATCH(NULL);
}

void
pypath_append(const std::string path)
{
  CYCPPUTIL_DEFAULT_TRY;

  PyObject * py_sys_path = PySys_GetObject("path"); // borrowed reference
  Py_unique_ptr py_path(string_to_pyobject(path));
  // New reference: https://stackoverflow.com/questions/3512414/does-this-pylist-appendlist-py-buildvalue-leak
  PyList_Append(py_sys_path, py_path.get());

  CYCPPUTIL_DEFAULT_CATCH();
}


} // end of cycpputil
#endif
