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

static PyObject * p3cbezier_bernstein(PyObject *self, PyObject *args);



//
// Method table and globals
//
static PyMethodDef CBezierMethods[] = {
 
    {"bernstein",  p3cbezier_bernstein, METH_VARARGS,
     "compute bernstein"},

 
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

#if CACHE_POWERS
static struct {FLOAT* list; size_t points; size_t power;} 
    cache_powers[CACHE_POWERS_SIZE];
#endif // CACHE_POWERS


#include "data/cbinom.data"
#define BINOM(n,m) (binomials[(((n)*((n)+1)) >> 1) + (m)])

static double pool[1<<22];



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

#define PX(i) ((i) << 1)
#define PY(i) (((i) << 1) + 1)
#define EXP(i,n) (npoints*(i)+(n))

static FLOAT polynom[1000];
PyObject * p3cbezier_bernstein(PyObject *self, PyObject *args) {
    unsigned i,j, npoints, length, degree;
    PyObject* polyarg, *temp, *list;
    FLOAT *pidx = polynom, *bidx = pool, k=0.0, step;

    // Parse args
    if (!PyArg_ParseTuple(args, "OI", &polyarg, &npoints)) return NULL;

    // Parse polynomial to polynom
    temp = PySequence_Fast(polyarg, "argument must be iterable");
    degree = (length = PySequence_Fast_GET_SIZE(temp)) - 1;
    for(i=0; i<length; i++) {
	list = PySequence_Fast(PySequence_Fast_GET_ITEM(temp, i), 
					       "points of list must be iterables");
	*pidx++ = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(list,0));
	*pidx++ = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(list,1));
	Py_DECREF(list);
    }
    Py_DECREF(temp);
    // Py_DECREF(polyarg);


    ////////////////////////////
    //  Compute exponentials  //
    ////////////////////////////



    step = 1.0/(npoints-1); 
    pidx = pool + npoints;
    for(; k<= 1.0; k+=step){ *pidx++ = k; *bidx++= 1.0;  }
    for(j=2; j<length; j++) for(k=0; k<=1; k+=step) { *pidx++ = (*bidx++)*k; }

    memset(pidx, 0, sizeof(double)*2*npoints); 
    for(i=0, bidx=pidx; i<= degree; i++, pidx = bidx) { 
	for(j=0; j < npoints; ++j) {
	    double aux = BINOM(degree,i) * pool[EXP(i,j)] * pool[EXP(degree - i,npoints-j-1)];
	    *pidx++ += polynom[PX(i)]* aux;
	    *pidx++ += polynom[PY(i)]* aux;
	    //printf("idx (%u,%u): [%fl,%fl]\n", i, j, *(pidx-2), *(pidx-1));
	}
    }

    list = PyList_New(npoints); Py_INCREF(list);
    for(j=0; j<npoints; ++j,bidx+=2){
	//double a1 = *bidx++, a2=*bidx++;
	PyList_SetItem(list, j, Py_BuildValue("[dd]", *bidx, *(bidx+1)));
	//printf("built point [%fl,%fl]\n", a1, a2);
    }

    //printf("returnning\n");

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
