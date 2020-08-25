#ifdef POLYVEC_BUILD_PYTHON

#include <cstdio>

#include <polyvec/cycpputil.hpp>


namespace cycpputil_test
{
void
dict(FILE *fl)
{

  const double eps = 1e-10;
  CYCPPUTIL_UNUSED( eps );

    CYCPPUTIL_DEFAULT_TRY;

    std::vector<int> item1 = {1,2,3,4,5,6};
    std::vector<std::string> item2 = {"hello", "dear", "amigos"};
    int item3 = 324;
    double item4 = 1000;
    std::string item5 = "why writing tests takes so long!!";

    cycpputil::Py_unique_ptr py_dict( PyDict_New() );

    cycpputil::pydict_absorb_item( py_dict.get(), "key_item1", cycpputil::vectori_to_pyobject( item1 ) );
    cycpputil::pydict_absorb_item( py_dict.get(), "key_item2", cycpputil::vectorstr_to_pyobject( item2 ) );
    cycpputil::pydict_absorb_item( py_dict.get(), "key_item3", cycpputil::int_to_pyobject( item3 ) );
    cycpputil::pydict_absorb_item( py_dict.get(), "key_item4", cycpputil::double_to_pyobject( item4 ) );
    cycpputil::pydict_absorb_item( py_dict.get(), "key_item5", cycpputil::string_to_pyobject( item5 ) );

    // Borrowed stuff no unique_ptr
    PyObject *py_item1 = cycpputil::pydict_get_item(py_dict.get(), "key_item1");
    PyObject *py_item2 = cycpputil::pydict_get_item(py_dict.get(), "key_item2");
    PyObject *py_item3 = cycpputil::pydict_get_item(py_dict.get(), "key_item3");
    PyObject *py_item4 = cycpputil::pydict_get_item(py_dict.get(), "key_item4");
    PyObject *py_item5 = cycpputil::pydict_get_item(py_dict.get(), "key_item5");

    CYCPPUTIL_ASSERT(  Py_REFCNT( py_item1 ) == 1 , "Py_REFCNT( py_item1 ) = " << Py_REFCNT( py_item1 ));
    CYCPPUTIL_ASSERT(  Py_REFCNT( py_item2 ) == 1 , "Py_REFCNT( py_item2 ) = " << Py_REFCNT( py_item2 ));
    CYCPPUTIL_ASSERT(  Py_REFCNT( py_item3 ) == 1 , "Py_REFCNT( py_item3 ) = " << Py_REFCNT( py_item3 ));
    CYCPPUTIL_ASSERT(  Py_REFCNT( py_item4 ) == 1 , "Py_REFCNT( py_item4 ) = " << Py_REFCNT( py_item4 ));
    CYCPPUTIL_ASSERT(  Py_REFCNT( py_item5 ) == 1 , "Py_REFCNT( py_item5 ) = " << Py_REFCNT( py_item5 ));

    // convert back to c type
    std::vector<int> c_py_item1 = cycpputil::pyobject_to_vectori( py_item1 );
    std::vector<std::string> c_py_item2 =  cycpputil::pyobject_to_vectorstr( py_item2 );
    int c_py_item3 =  cycpputil::pyobject_to_int( py_item3 );
    double c_py_item4 =  cycpputil::pyobject_to_double( py_item4 );
    std::string c_py_item5 =  cycpputil::pyobject_to_string( py_item5 );

    CYCPPUTIL_ASSERT( c_py_item1 == item1, "");
    CYCPPUTIL_ASSERT( c_py_item2 == item2, "");
    CYCPPUTIL_ASSERT( c_py_item3 == item3, "");
    CYCPPUTIL_ASSERT( std::abs( c_py_item4 - item4 ) < eps, "");
    CYCPPUTIL_ASSERT( c_py_item5 == item5, "");

    CYCPPUTIL_DEFAULT_CATCH();
}
}

#endif