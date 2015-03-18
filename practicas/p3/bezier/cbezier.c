#include <Python.h>
#include <stdlib.h>
#include <stdio.h>

// Configuration constants

#define CACHE_POWERS true
#define CACHE_POWERS_SIZE 4
typedef double  FLOAT;

//
// Prototype declarations
//

static PyObject * cbezier_system(PyObject *self, PyObject *args);
static PyObject * cbezier_bernstein(PyObject *self, PyObject *args);



//
// Method table and globals
//
static PyMethodDef CBezierMethods[] = {
 
    {"system",  cbezier_system, METH_VARARGS,
     "Execute a shell command."},

    {"bernstein",  cbezier_bernstein, METH_VARARGS,
     "compute bernstein"},

 
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

#if CACHE_POWERS
static struct {FLOAT* list; size_t points; size_t power;} 
    cache_powers[CACHE_POWERS_SIZE];
#endif // CACHE_POWERS

// static binom


//
// Module initialization
//
PyMODINIT_FUNC initcbezier(void) {
    PyObject *m;

    m = Py_InitModule("cbezier", CBezierMethods);
    if (m == NULL)
        return;

    /* SpamError = PyErr_NewException("spam.error", NULL, NULL); */
    /* Py_INCREF(SpamError); */
    /* PyModule_AddObject(m, "error", SpamError); */
}



//
// Method implementation
//

PyObject * cbezier_system(PyObject *self, PyObject *args) {
    const char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return Py_BuildValue("i", sts);
}



PyObject * cbezier_bernstein(PyObject *self, PyObject *args) {
    int i;
    unsigned npoints;
    struct PyObject* polyarg, *temp;

    // Parse args
    if (!PyArg_ParseTuple(args, "OI", &polyarg, &npoints)) return NULL;

    temp = PySequence_Fast(polyarg, "argument must be iterable");
    size_t length = PySequence_Fast_GET_SIZE(temp);
    FLOAT* array1 = malloc((length << 1) * sizeof(FLOAT));

    // Parse polynomial to array1
    // printf("arg1: %uL; l: %u\n", length, npoints);
    for(i=0; i<length; i++) {
	struct PyObject* aux = PySequence_Fast(PySequence_Fast_GET_ITEM(temp, i), 
					       "points of list must be iterables");
	array1[i << 1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(aux,0));
	array1[(i << 1) + 1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(aux,1));
	Py_DECREF(aux);
    }



    // Use array1 to compute bernstein

    // do magic with array1!!!
    // array1[i << 1], array1[(i << 1) + 1]


    // Create output list of points as python object

    PyObject *list = PyList_New(length);
    Py_INCREF(list);

    for(i=0; i<length; i++) {
	PyList_SetItem(list, i, Py_BuildValue("[dd]", array1[i << 1], array1[(i << 1) + 1]));
    }

    Py_DECREF(temp);
    free(array1);

    //Py_RETURN_NONE;
    return list;
}

/*
temp = PySequence_Fast(arg1, "argument must be iterable");
if (!temp) {
    return NULL;
}
length1 = PySequence_Fast_GET_SIZE(temp);
array1 = malloc(length1 * sizeof(long);
if (!array1) {
    Py_DECREF(temp);
    return PyErr_NoMemory();
}
for (i=0; i < length1; i++) {
    t = PyLong_AsLong(PySequence_Fast_GET_ITEM(temp, i));
    if (t == -1 && PyErr_Occurred()) {
        Py_DECREF(temp);
        free(array1);
        return NULL;
    }
    array1[i] = temp;
}

# Do magic stuff.
# Assume array1 contents have changed and need to be returned.

result1 = PyTuple_New(length1);
if (!result1) {
    Py_DECREF(temp);
    free(array1);
    return NULL;
}
for (i=0; i < length1; i++) {
    z = PyLong_FromLong(array1[i]);
    if (!z)  {
        Py_DECREF(temp);
        free(array1);
        return NULL;
    }
    PyTuple_SET_ITEM(result1, i, z);
}
*/
