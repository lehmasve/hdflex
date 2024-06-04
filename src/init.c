#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _hdflex_active_models_dsc(void *, void *);
extern SEXP _hdflex_agg_density_dsc(void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_dsc_loop(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_forget_dsc(void *, void *);
extern SEXP _hdflex_init_dsc(void *);
extern SEXP _hdflex_matrix_subset_idx(void *, void *, void *);
extern SEXP _hdflex_stsc_loop_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_stsc_loop_par_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_tvc_(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_update_dsc(void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_hdflex_active_models_dsc", (DL_FUNC) &_hdflex_active_models_dsc,  2},
    {"_hdflex_agg_density_dsc",   (DL_FUNC) &_hdflex_agg_density_dsc,    6},
    {"_hdflex_dsc_loop",          (DL_FUNC) &_hdflex_dsc_loop,           9},
    {"_hdflex_forget_dsc",        (DL_FUNC) &_hdflex_forget_dsc,         2},
    {"_hdflex_init_dsc",          (DL_FUNC) &_hdflex_init_dsc,           1},
    {"_hdflex_matrix_subset_idx", (DL_FUNC) &_hdflex_matrix_subset_idx,  3},
    {"_hdflex_stsc_loop_",        (DL_FUNC) &_hdflex_stsc_loop_,        18},
    {"_hdflex_stsc_loop_par_",    (DL_FUNC) &_hdflex_stsc_loop_par_,    19},
    {"_hdflex_tvc_",              (DL_FUNC) &_hdflex_tvc_,               7},
    {"_hdflex_update_dsc",        (DL_FUNC) &_hdflex_update_dsc,         6},
    {NULL, NULL, 0}
};

void R_init_hdflex(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
