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

static PyObject * p3cbezier_system(PyObject *self, PyObject *args);
static PyObject * p3cbezier_bernstein(PyObject *self, PyObject *args);



//
// Method table and globals
//
static PyMethodDef CBezierMethods[] = {
 
    {"system",  p3cbezier_system, METH_VARARGS,
     "Execute a shell command."},

    {"bernstein",  p3cbezier_bernstein, METH_VARARGS,
     "compute bernstein"},

 
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

#if CACHE_POWERS
static struct {FLOAT* list; size_t points; size_t power;} 
    cache_powers[CACHE_POWERS_SIZE];
#endif // CACHE_POWERS


#include "binom/cbinom.data"
#define BINOM(n,m) (binomials[(((n)*((n)+1)) >> 1) + (m)])

#include "exponentials.data"

static double pool[1<<15];
static unsigned points_list[1<<15];


//
// Module initialization
//
PyMODINIT_FUNC initp3cbezier(void) {
    PyObject *m;
    m = Py_InitModule("p3cbezier", CBezierMethods);
    if (m == NULL) return;
}



//
// Method implementation
//

PyObject * p3cbezier_system(PyObject *self, PyObject *args) {
    const char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return Py_BuildValue("i", sts);
}

#define PX(i) ((i) << 1)
#define PY(i) (((i) << 1) + 1)

PyObject * p3cbezier_bernstein(PyObject *self, PyObject *args) {
    int i;
    unsigned npoints;
    struct PyObject* polyarg, *temp;

    // Parse args
    if (!PyArg_ParseTuple(args, "OI", &polyarg, &npoints)) return NULL;

    temp = PySequence_Fast(polyarg, "argument must be iterable");
    size_t length = PySequence_Fast_GET_SIZE(temp);
    FLOAT* polynom = malloc((length << 1) * sizeof(FLOAT));

    // Parse polynomial to polynom
    // printf("arg1: %uL; l: %u\n", length, npoints);
    for(i=0; i<length; i++) {
	struct PyObject* aux = PySequence_Fast(PySequence_Fast_GET_ITEM(temp, i), 
					       "points of list must be iterables");
	polynom[i << 1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(aux,0));
	polynom[(i << 1) + 1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(aux,1));
	Py_DECREF(aux);
    }
    // Py_DECREF(polyarg);



    // Use polynom to compute bernstein

    // do magic with polynom!!!
    // polynom[i << 1], polynom[(i << 1) + 1]


        /*     m = len(t) - 1  */
        /* d, r = cls.limit / m , cls.limit % m */
        /* def gen(end, step, off, points):  */
        /*     i,r,no = 0, 0, 0 */
        /*     while r + no <= end: yield r + no; i += 1; r += step ; no = (i*off)/points */
        /* return cls.precalc[:n+1,list(gen(cls.limit, d, r, m))] */

    PyObject *list = PyList_New(npoints);
    Py_INCREF(list);


    unsigned degree = length - 1, m = npoints - 1, d = (jlimit-1)/m, r = (jlimit-1) % m;
    unsigned j=0, k=0, jj=0, no = 0;

    for(; j < npoints; ++j, k+=d, (no=(j*r)/m), (jj=k+no)) points_list[j] = jj;

    for(j=0,jj=0, k=points_list[npoints-1]; j < npoints; ++j) {
	// printf("index %u ->", j); 
	jj=points_list[j]; k=points_list[npoints-j-1];
	// printf(" point %u\n", jj);

	double rx = 0, ry = 0;
	for(i=0; i<= degree; i++) {
	    // printf("\t%f %f %f\n", BINOM(degree,i), exponentials[i][jj], exponentials[degree - i][k]);

	    double aux = BINOM(degree,i) * exponentials[i][jj] * exponentials[degree - i][k];
	    rx += polynom[PX(i)]* aux;
	    ry += polynom[PY(i)]* aux;
	}
	PyList_SetItem(list, j, Py_BuildValue("[dd]", rx, ry));
	// printf("[%f,%f]\n", rx,ry);
    }

    // printf("asdfdsfa\n");

    /* binomials[0]; */
    /* BINOM(4,2); */
     
    /* ilimit; jlimit; exponentials; */



    // Create output list of points as python object

    /* PyObject *list = PyList_New(length); */
    /* Py_INCREF(list); */

    /* for(i=0; i<length; i++) { */
    /* 	PyList_SetItem(list, i, Py_BuildValue("[dd]", polynom[i << 1], polynom[(i << 1) + 1])); */
    /* } */


    // printf("pre-returning\n");

    Py_DECREF(temp);
    // free(polynom);

    // printf("returning\n");

    //Py_RETURN_NONE;
    return list;
}

/*
temp = PySequence_Fast(arg1, "argument must be iterable");
if (!temp) {
    return NULL;
}
length1 = PySequence_Fast_GET_SIZE(temp);
polynom = malloc(length1 * sizeof(long);
if (!polynom) {
    Py_DECREF(temp);
    return PyErr_NoMemory();
}
for (i=0; i < length1; i++) {
    t = PyLong_AsLong(PySequence_Fast_GET_ITEM(temp, i));
    if (t == -1 && PyErr_Occurred()) {
        Py_DECREF(temp);
        free(polynom);
        return NULL;
    }
    polynom[i] = temp;
}

# Do magic stuff.
# Assume polynom contents have changed and need to be returned.

result1 = PyTuple_New(length1);
if (!result1) {
    Py_DECREF(temp);
    free(polynom);
    return NULL;
}
for (i=0; i < length1; i++) {
    z = PyLong_FromLong(polynom[i]);
    if (!z)  {
        Py_DECREF(temp);
        free(polynom);
        return NULL;
    }
    PyTuple_SET_ITEM(result1, i, z);
}
*/
