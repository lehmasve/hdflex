#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _hdflex_dsc_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_stsc_loop_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_stsc_loop_par_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_tvc_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_hdflex_dsc_",           (DL_FUNC) &_hdflex_dsc_,           12},
    {"_hdflex_stsc_loop_",     (DL_FUNC) &_hdflex_stsc_loop_,     16},
    {"_hdflex_stsc_loop_par_", (DL_FUNC) &_hdflex_stsc_loop_par_, 17},
    {"_hdflex_tvc_",           (DL_FUNC) &_hdflex_tvc_,            7},
    {NULL, NULL, 0}
};

void R_init_hdflex(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

