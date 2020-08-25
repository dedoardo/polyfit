#ifdef POLYVEC_BUILD_PYTHON

#include <cstdio>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

// Sadly forward declaration is not possible
#include <Python.h>

namespace cycpputil
{

struct Error
{
  std::vector<std::string> funcs = std::vector<std::string>();
  std::string msg = "";
  void report(std::ostream &);
  void report(FILE *);
};

struct Py_init
{
  Py_init();
  ~Py_init();
};

// Use for temp objects that have to be DecRef'd at the
// end
struct Py_unique_ptr
{
  Py_unique_ptr(PyObject * val_in = nullptr);
  ~Py_unique_ptr();
  void reset(PyObject * val_in = nullptr);
  PyObject * get();

private:
  PyObject * _val;
};


// ===================================================================
// Vectors
// ===================================================================

std::vector<int>
pyobject_to_vectori(PyObject *);

// Ref Info: new reference
PyObject *
vectori_to_pyobject(const std::vector<int> &);

std::vector<double>
pyobject_to_vectord(PyObject *);


// Ref Info: new reference
PyObject *
vectord_to_pyobject(const std::vector<double> &);


std::vector<std::string>
pyobject_to_vectorstr(PyObject *);


// Ref Info: new reference
PyObject *
vectorstr_to_pyobject(const std::vector<std::string> &);

//
// TODO
// PyObject *
// vectorpyobject_to_pyobject(const std::vector<std::string> &);


// ===================================================================
// Int and double
// ===================================================================

int
pyobject_to_int(PyObject * obj);

// Ref Info: new reference
PyObject *
int_to_pyobject(const int);

double
pyobject_to_double(PyObject * obj);

// Ref Info: new reference
PyObject *
double_to_pyobject(const double);

std::string
pyobject_to_string(PyObject * obj);

// Ref Info: new reference
PyObject *
string_to_pyobject(const std::string);

// ===================================================================
// Dictionary
// ===================================================================

bool
pydict_has_item(PyObject * py_dict, const std::string & key);

// Ref Info: Stolen reference
PyObject *
pydict_get_item(PyObject * py_dict, const std::string & key);

// Ref Info: we will create a new reference of the item in the dict
bool
pydict_set_item(PyObject * py_dict, const std::string & key, PyObject * item);

// Ref Info: the item will be absorbed and set to NULL
bool
pydict_absorb_item(PyObject * py_dict, const std::string & key, PyObject * item);

// ===================================================================
// Callable objects
// NOTE: PYTHON PATH MUST BE SET CORRECTLY.
// ===================================================================

// Ref Info: returns New reference
PyObject *
pycall_load_module(const std::string module_name);

// Ref Info: returns New reference
PyObject *
pycall_load_function(PyObject * py_module, const std::string function_name);

// Ref Info: returns New reference
PyObject *
pycall_call_function(PyObject * py_function, PyObject * py_args, PyObject * py_kwargs);

// ===================================================================
// Append a string to python system path
// ===================================================================

void
pypath_append(const std::string path);

} // end of cycpputil


// ===================================================================
// Macros for error reporting
// ===================================================================

#define CYCPPUTIL_STR(stuff) static_cast<std::ostringstream &>(std::ostringstream().flush() << stuff).str()

#define CYCPPUTIL_C_STR(stuff) (CYCPPUTIL_STR(stuff).c_str())

#define CYCPPUTIL_ASSERT(__cond, __msg)                                                                             \
  if(!(__cond))                                                                                                     \
  {                                                                                                                 \
    ::cycpputil::Error err;                                                                                         \
    err.msg = CYCPPUTIL_STR(                                                                                        \
        "CONDITION: " << #__cond << "\n MSG:" << __msg << "\n " << __FILE__ << " " << __func__ << " " << __LINE__); \
    assert(0);                                                                                                      \
    throw err;                                                                                                      \
  }

#define CYCPPUTIL_PROPAGATE_ERR(err)                                   \
  do                                                                   \
  {                                                                    \
    err.funcs.push_back(CYCPPUTIL_C_STR(__FILE__ << " " << __func__)); \
    throw err;                                                         \
  } while(0)

#define CYCPPUTIL_TERMINATE_ERR(err, ostr)                             \
  do                                                                   \
  {                                                                    \
    err.funcs.push_back(CYCPPUTIL_C_STR(__FILE__ << " " << __func__)); \
    err.report(ostr);                                                  \
  } while(0)

#define CYCPPUTIL_DEFAULT_TRY \
  try                     \
  {

#define CYCPPUTIL_DEFAULT_CATCH(DEFAULT_RETURN) \
  }                                         \
  catch(::cycpputil::Error ierr)            \
  {                                         \
    CYCPPUTIL_PROPAGATE_ERR(ierr);          \
    return DEFAULT_RETURN;                  \
  }


#define CYCPPUTIL_UNUSED(_X) (void)_X;

#endif