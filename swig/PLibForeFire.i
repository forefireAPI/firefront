/* PLibForeFire */

%module forefire
%{
#define SWIG_FILE_WITH_INIT
#include "PLibForeFire.h"

%}

%fragment("NumPy_Backward_Compatibility", "header")
{

}
%typemap(out) std::string {
    $result = PyString_FromString($1.c_str());
}

%pythoncode %{

%}


%include <numpy.i>

%init %{
import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {( int mdimx,  double* meshx), ( int mdimy,  double* meshy), ( int mdimz,  double* zgrid)};
%apply (int DIM1, int DIM2, int DIM3, double* IN_ARRAY3){(int nnx, int nny, int nnz,  double* values)};
%apply (int DIM1, int DIM2, int DIM3, int* IN_ARRAY3){(int nnx, int nny, int nnz,  int* values)};
%apply (double ARGOUT_ARRAY3[ANY][ANY][ANY]) {(double lower[2][2][2])};
%apply (double** ARGOUTVIEW_ARRAY3, int* DIM1, int* DIM2, int* DIM3){(double** outA, int* outNI, int* outNJ, int* outNK)};

// Include the header file to be wrapped

%include "PLibForeFire.h"

