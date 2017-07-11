#include <Python.h>

/*
 * Implements an example function.
 */
PyDoc_STRVAR(IsingCore_example_doc, "example(obj, number)\
\
Example function");

PyObject *IsingCore_example(PyObject *self, PyObject *args, PyObject *kwargs) {
    /* Shared references that do not need Py_DECREF before returning. */
    PyObject *obj = NULL;
    int number = 0;

    /* Parse positional and keyword arguments */
    static char* keywords[] = { "obj", "number", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", keywords, &obj, &number)) {
        return NULL;
    }

    /* Function implementation starts here */

    if (number < 0) {
        PyErr_SetObject(PyExc_ValueError, obj);
        return NULL;    /* return NULL indicates error */
    }

    Py_RETURN_NONE;
}

/*
 * List of functions to add to IsingCore in exec_IsingCore().
 */
static PyMethodDef IsingCore_functions[] = {
    { "example", (PyCFunction)IsingCore_example, METH_VARARGS | METH_KEYWORDS, IsingCore_example_doc },
    { NULL, NULL, 0, NULL } /* marks end of array */
};

/*
 * Initialize IsingCore. May be called multiple times, so avoid
 * using static state.
 */
int exec_IsingCore(PyObject *module) {
    PyModule_AddFunctions(module, IsingCore_functions);

    PyModule_AddStringConstant(module, "__author__", "Megan");
    PyModule_AddStringConstant(module, "__version__", "1.0.0");
    PyModule_AddIntConstant(module, "year", 2017);

    return 0; /* success */
}

/*
 * Documentation for IsingCore.
 */
PyDoc_STRVAR(IsingCore_doc, "The IsingCore module");


static PyModuleDef_Slot IsingCore_slots[] = {
    { Py_mod_exec, exec_IsingCore },
    { 0, NULL }
};

static PyModuleDef IsingCore_def = {
    PyModuleDef_HEAD_INIT,
    "IsingCore",
    IsingCore_doc,
    0,              /* m_size */
    NULL,           /* m_methods */
    IsingCore_slots,
    NULL,           /* m_traverse */
    NULL,           /* m_clear */
    NULL,           /* m_free */
};

PyMODINIT_FUNC PyInit_IsingCore() {
    return PyModuleDef_Init(&IsingCore_def);
}
