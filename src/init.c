#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _hdflex_active_models_dsc(SEXP, SEXP, SEXP);
extern SEXP _hdflex_agg_density_dsc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_dsc_loop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_forget_dsc(SEXP, SEXP);
extern SEXP _hdflex_init_dsc(SEXP);
extern SEXP _hdflex_init_tvc(SEXP, SEXP, SEXP);
extern SEXP _hdflex_init_tvc_forecast(SEXP, SEXP);
extern SEXP _hdflex_matrix_subset_idx(SEXP, SEXP, SEXP);
extern SEXP _hdflex_tvc_model(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_tvc_model_forecasts(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_tvc_model_loop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_tvc_model_loop_forecasts(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _hdflex_update_dsc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_hdflex_active_models_dsc",        (DL_FUNC) &_hdflex_active_models_dsc,         3},
    {"_hdflex_agg_density_dsc",          (DL_FUNC) &_hdflex_agg_density_dsc,           6},
    {"_hdflex_dsc_loop",                 (DL_FUNC) &_hdflex_dsc_loop,                 10},
    {"_hdflex_forget_dsc",               (DL_FUNC) &_hdflex_forget_dsc,                2},
    {"_hdflex_init_dsc",                 (DL_FUNC) &_hdflex_init_dsc,                  1},
    {"_hdflex_init_tvc",                 (DL_FUNC) &_hdflex_init_tvc,                  3},
    {"_hdflex_init_tvc_forecast",        (DL_FUNC) &_hdflex_init_tvc_forecast,         2},
    {"_hdflex_matrix_subset_idx",        (DL_FUNC) &_hdflex_matrix_subset_idx,         3},
    {"_hdflex_tvc_model",                (DL_FUNC) &_hdflex_tvc_model,                 8},
    {"_hdflex_tvc_model_forecasts",      (DL_FUNC) &_hdflex_tvc_model_forecasts,       7},
    {"_hdflex_tvc_model_loop",           (DL_FUNC) &_hdflex_tvc_model_loop,           10},
    {"_hdflex_tvc_model_loop_forecasts", (DL_FUNC) &_hdflex_tvc_model_loop_forecasts,  9},
    {"_hdflex_update_dsc",               (DL_FUNC) &_hdflex_update_dsc,                6},
    {NULL, NULL, 0}
};

void R_init_hdflex(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
