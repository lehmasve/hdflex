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
extern SEXP _hdflex_dsc_active_models_(void *, void *);
extern SEXP _hdflex_dsc_agg_density_(void *, void *, void *, void *);
extern SEXP _hdflex_dsc_dpll_cands_(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_dsc_dpll_comb_(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_dsc_init_(void *, void *, void *, void *);
extern SEXP _hdflex_dsc_loop(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_dsc_loop_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_forget_dsc(void *, void *);
extern SEXP _hdflex_init_dsc(void *);
extern SEXP _hdflex_init_tvc(void *, void *, void *);
extern SEXP _hdflex_init_tvc_(void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_init_tvc_forecast(void *, void *, void *);
extern SEXP _hdflex_matrix_subset_idx(void *, void *, void *);
extern SEXP _hdflex_rank_comb_(void *, void *, void *);
extern SEXP _hdflex_stsc_loop(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_tvc_model(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_tvc_model_(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_tvc_model_cand_(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_tvc_model_loop(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _hdflex_update_dsc(void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_hdflex_active_models_dsc",  (DL_FUNC) &_hdflex_active_models_dsc,   2},
    {"_hdflex_agg_density_dsc",    (DL_FUNC) &_hdflex_agg_density_dsc,     6},
    {"_hdflex_dsc_active_models_", (DL_FUNC) &_hdflex_dsc_active_models_,  2},
    {"_hdflex_dsc_agg_density_",   (DL_FUNC) &_hdflex_dsc_agg_density_,    4},
    {"_hdflex_dsc_dpll_cands_",    (DL_FUNC) &_hdflex_dsc_dpll_cands_,     9},
    {"_hdflex_dsc_dpll_comb_",     (DL_FUNC) &_hdflex_dsc_dpll_comb_,      9},
    {"_hdflex_dsc_init_",          (DL_FUNC) &_hdflex_dsc_init_,           4},
    {"_hdflex_dsc_loop",           (DL_FUNC) &_hdflex_dsc_loop,            9},
    {"_hdflex_dsc_loop_",          (DL_FUNC) &_hdflex_dsc_loop_,          13},
    {"_hdflex_forget_dsc",         (DL_FUNC) &_hdflex_forget_dsc,          2},
    {"_hdflex_init_dsc",           (DL_FUNC) &_hdflex_init_dsc,            1},
    {"_hdflex_init_tvc",           (DL_FUNC) &_hdflex_init_tvc,            3},
    {"_hdflex_init_tvc_",          (DL_FUNC) &_hdflex_init_tvc_,           6},
    {"_hdflex_init_tvc_forecast",  (DL_FUNC) &_hdflex_init_tvc_forecast,   3},
    {"_hdflex_matrix_subset_idx",  (DL_FUNC) &_hdflex_matrix_subset_idx,   3},
    {"_hdflex_rank_comb_",         (DL_FUNC) &_hdflex_rank_comb_,          3},
    {"_hdflex_stsc_loop",          (DL_FUNC) &_hdflex_stsc_loop,          16},
    {"_hdflex_tvc_model",          (DL_FUNC) &_hdflex_tvc_model,           8},
    {"_hdflex_tvc_model_",         (DL_FUNC) &_hdflex_tvc_model_,          8},
    {"_hdflex_tvc_model_cand_",    (DL_FUNC) &_hdflex_tvc_model_cand_,     8},
    {"_hdflex_tvc_model_loop",     (DL_FUNC) &_hdflex_tvc_model_loop,     10},
    {"_hdflex_update_dsc",         (DL_FUNC) &_hdflex_update_dsc,          6},
    {NULL, NULL, 0}
};

void R_init_hdflex(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
