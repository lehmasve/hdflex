#include <RcppArmadillo.h>
#include <thread>
#include <RcppThread.h>
using namespace Rcpp;
using namespace RcppThread;

// 1) TV-C Regressions
// Function I - Initialize TV-C-Parameter
List init_tvc_par_(const arma::vec& y, 
                   const arma::mat& S,
                   int n_raw_sig,
                   int sample_length,
                   arma::vec lambda_grid,
                   arma::vec kappa_grid,
                   int n_threads) {
  
  // Get Dimensions   
  const int n_signal = S.n_cols;
  const int n_cands  = n_signal * lambda_grid.n_elem * kappa_grid.n_elem;
  
  // Define Variables 
  arma::cube   theta_cube(2, 1, n_cands);
  arma::cube   cov_mat_cube(2, 2, n_cands);
  arma::rowvec h_vec(n_cands);
  arma::vec    y_sample, cm_na_idx(n_cands); cm_na_idx.fill(arma::datum::nan); // x 
  List         ret_all(4);  
  
  // Define Counter Cube
  arma::cube ctr_cube(lambda_grid.n_elem, kappa_grid.n_elem, n_signal);
  unsigned int ctr = 0;
  // fill ctr cube according to loop
  for (unsigned int l = 0; l < lambda_grid.n_elem; l++) {
    for (unsigned int k = 0; k < kappa_grid.n_elem; k++) {
      for (unsigned int j = 0; j < n_signal; j++) {
        ctr_cube(l, k, j) = ctr;
        ctr++;
      }
    }
  }
  
  // Parallel Loop over lambdas
  parallelFor(0, lambda_grid.n_elem, [&y, &S, &n_cands, &sample_length, 
              &kappa_grid, &n_signal, &n_raw_sig,
              &theta_cube, &cov_mat_cube, &h_vec,
              &ctr_cube] (unsigned int l) {
              
              // Define variables
              arma::vec     x;
                arma::uvec    init_idx, non_finite;
                arma::vec     y_sample, cm_na_idx(n_cands); cm_na_idx.fill(arma::datum::nan); 
                arma::mat     x_sample_one, theta, cov_mat;
                arma::colvec  coef;
                double        intercept, var_y, var_x, h;
                
                // Loop over Kappas and Sinals
                for (unsigned int k = 0; k < kappa_grid.n_elem; k++) {
                  for (unsigned int j = 0; j < n_signal; j++) {
                                        
                    // Select Signal
                    x = S.col(j);
                    
                    // Check and Count for NA-Values
                    non_finite = arma::find_nonfinite(x);
                    int na_ctr = non_finite.n_elem;
                    if (na_ctr > 0) {
                      // get lth lement from ctr
                      cm_na_idx(ctr_cube(l, k, j)) = ctr_cube(l, k, j);
                    }
                    
                    // Index for subsetting the Initialisation Sample
                    init_idx = arma::regspace<arma::uvec>(0 + na_ctr, na_ctr + sample_length - 1);
                    
                    // Define and prepare matrices for regression  
                    y_sample = y.elem(init_idx);  
                    arma::mat x_sample(x.elem(init_idx));   
                    x_sample_one = arma::ones<arma::mat>(x_sample.n_rows, 1);
                    x_sample.insert_cols(0, x_sample_one);
                                        
                    // Initialize - Theta
                    theta = arma::zeros<arma::mat>(2,1);
                    
                    // Initialize - System Covariance
                    coef      = solve(x_sample, y_sample);
                    intercept = coef(0);
                    var_y     = arma::var(y_sample);         
                    var_x     = arma::var(x_sample.col(1));  
                    cov_mat   = arma::zeros<arma::mat>(2, 2); 
                    cov_mat(0, 0) = pow(intercept, 2) + var_y; // Set to Zero for Constant Intercept
                    
                    // Distinguish between raw and processed signals
                    if(j < n_raw_sig) {
                      if(arma::is_finite(var_y / var_x)) {
                        cov_mat(1, 1) =  var_y / var_x; 
                      } else {
                        cov_mat(1, 1) = 0;
                      }
                    } else {
                      theta(1, 0) = 1;
                      cov_mat(1, 1) = 0;
                    }
                    
                    // Initialize - Observational Variance
                    h = var_y;
                    
                    // Fill Cubes  
                    unsigned int index = ctr_cube(l, k, j);
                    theta_cube.slice(index)   = theta;
                    cov_mat_cube.slice(index) = cov_mat; 
                    h_vec(index)              = h;
                    //ctr++;
                  }
                }
              }, n_threads);
  
  // Fill Return-List  
  ret_all[0] = theta_cube;
  ret_all[1] = cov_mat_cube; 
  ret_all[2] = h_vec;
  ret_all[3] = arma::conv_to<arma::uvec>::from(cm_na_idx.elem(arma::find_finite(cm_na_idx)));
  
  // Return Results
  return ret_all;
}
// ----------


// 1) TV-C Regressions
// Function II - Predictive and Update Step for Signal j and time t
arma::field<double> tvc_model_par_(double y_t, 
                                   double s_t_j, 
                                   double s_pred_j,
                                   double lambda, 
                                   double kappa, 
                                   arma::mat& theta, 
                                   arma::mat& cov_mat, 
                                   double& h) { 
  
  // Define Variables
  arma::field<double> ret(2);
  //arma::mat r_upt;
  //double mu, variance, inv_tvar, e_t;
  
  // Get Predictor for Time t and Predicting t + 1
  const arma::mat z_t = {1.0, s_t_j};
  const arma::mat z_pred = {1.0, s_pred_j};
  
  // Add noise to uncertainty of Coefficients in time t (Equation 5)
  const arma::mat r_upt = cov_mat / lambda;
  
  // Calculate (OOS) Forecast Error for time t (see Equation 7)
  const double e_t = arma::as_scalar(y_t - z_t * theta);
  
  // Update Observational Variance in time t (Equation 10 and 11)
  h = arma::as_scalar(kappa * h + (1 - kappa) * pow(e_t, 2));
  
  // Update Coefficients in time t (Equation 7)
  const double inv_tvar = arma::as_scalar(1.0 / (h + z_t * r_upt * z_t.t()));
  theta = theta + r_upt * z_t.t() * inv_tvar * e_t;
  
  // Update Uncertainty of Coefficients in time t (Equation 8)
  cov_mat = r_upt - r_upt * z_t.t() * inv_tvar * (z_t * r_upt);
  
  // Get Predictive Density for Predicting t + 1 (Equation 9)
  const double       mu = arma::as_scalar(z_pred * theta);
  const double variance = arma::as_scalar(h + z_pred * ((1.0 / lambda) * cov_mat) * z_pred.t());
  
  // Fill Return-List  
  ret(0) = mu;
  ret(1) = variance;
  
  // Return List  
  return ret;
}


// 1) TV-C Regressions
// Function III - Loop over Signals, Lambda and Kappa
arma::field<arma::rowvec> tvc_model_cand_par_(double y_t, 
                                              const arma::rowvec& s_t,
                                              const arma::rowvec& s_pred, 
                                              arma::vec lambda_grid, 
                                              arma::vec kappa_grid, 
                                              arma::cube& theta_cube, 
                                              arma::cube& cov_mat_cube, 
                                              arma::rowvec& h_vec,
                                              int n_threads) {
  
  // Get Dimensions
  const int n_cands = s_t.n_elem * lambda_grid.n_elem * kappa_grid.n_elem;
  
  // Define Variables
  arma::field<arma::rowvec> ret(2);
  arma::rowvec mu_vec(n_cands);
  arma::rowvec variance_vec(n_cands);
  
  // Loop over Lambda
  // Define Counter Cube
  arma::cube ctr_cube(lambda_grid.n_elem, kappa_grid.n_elem, s_t.n_elem);
  unsigned int ctr = 0;
  // fill ctr cube according to loop
  for (unsigned int l = 0; l < lambda_grid.n_elem; l++) {
    for (unsigned int k = 0; k < kappa_grid.n_elem; k++) {
      for (unsigned int j = 0; j < s_t.n_elem; j++) {
        ctr_cube(l, k, j) = ctr;
        ctr++;
      }
    }
  }
    
  // Parallel Loop over lambdas
  parallelFor(0, lambda_grid.n_elem, [&y_t, &ctr_cube, &lambda_grid, &kappa_grid,
                                      &s_t, &s_pred, &theta_cube,
                                      &cov_mat_cube, &h_vec,
                                      &mu_vec, &variance_vec] (unsigned int l) {

    // Set Lambda
    const double lambda = lambda_grid(l);

    // Loop over Kappa
    for (unsigned int k = 0; k < kappa_grid.n_elem; k++) {

      // Set Kappa
      const double kappa =  kappa_grid(k);

      // Loop over all candidates
      for (unsigned int j = 0; j < s_t.n_elem; j++) {

        // Get Counter
        unsigned int counter = ctr_cube(l, k, j);

        // Set Signals
        const double    s_t_j = s_t(j);
        const double s_pred_j = s_pred(j);

        // Check if signal is NA or not
        const bool is_na = !arma::is_finite(s_t_j);
        if(!is_na) {

          // Apply TV-C-Function
          const arma::field<double> tvc_results = tvc_model_par_(y_t,
                                                                 s_t_j,
                                                                 s_pred_j,
                                                                 lambda,
                                                                 kappa,
                                                                 theta_cube.slice(counter),
                                                                 cov_mat_cube.slice(counter),
                                                                 h_vec(counter));
          // Assign TV-C-Model-Results
          mu_vec(counter)             = tvc_results(0);
          variance_vec(counter)       = tvc_results(1);

        } else {

          // Assign TV-C-Model-Results
          mu_vec(counter)             = arma::datum::nan;
          variance_vec(counter)       = arma::datum::nan;
        }

      }
    }
  }, n_threads);
  
  // Fill list 
  ret(0) = mu_vec;
  ret(1) = variance_vec;
  
  // Return list
  return ret;
}

// ###################################################################################
// ###################################################################################
// ###################################################################################

// 2) Dynamic Subset Combination
// Function I - Initialize DSC-Parameter
arma::field<arma::field<arma::rowvec>> dsc_init_par_(int n_cands,
                                                     int n_combs,
                                                     int n_gamma,
                                                     arma::uvec na_idx) {
  
  // Define Variables
  arma::field<arma::field<arma::rowvec>> ret(2); 
  
  // Initialize Vector for Discounted-Log-Likelihood-Score (Subset Combinations)  
  arma::rowvec dpll_combs(n_combs, arma::fill::zeros);
  
  // Initialize Vector for Discounted-Log-Likelihood-Score (Candidate Models)  
  arma::rowvec vec(n_cands, arma::fill::zeros);
  vec.elem(na_idx).fill(arma::datum::nan);
  
  // Fill Field for Candidate Models
  arma::field<arma::rowvec> dpll_cands(n_gamma);
  for (unsigned int i = 0; i < n_gamma; ++i) {
    dpll_cands(i) = vec;
  }
  
  // Fill Return-List 
  ret(0) = dpll_cands;
  ret(1) = arma::field<arma::rowvec>(1);
  ret(1)(0) = dpll_combs;
  
  // Return Vector
  return ret;
}

// Function II - Rank and Set Active Model Set (Active Models)
arma::field<arma::uvec> dsc_active_models_par_(const arma::rowvec& dpll_cands_gamma,  
                                               int psi) {
  
  // Get indices of Non-NaN values
  const arma::uvec non_na = arma::find_finite(dpll_cands_gamma);
  const int non_na_ctr = non_na.n_elem;
  const int psi_ = std::min(non_na_ctr, psi);
  
  // Check if all elements are equal (without NAs)
  const arma::vec vec = dpll_cands_gamma.elem(non_na);
  const bool isequal = arma::all(vec == vec(0));
  
  // If all elements are equal, return the indices from left to right
  arma::uvec idx;
  if (isequal) {
    idx = non_na.head(psi_);
  } else {
    // Get psi-highest values (-> indices)
    const arma::uvec sorted_idx = arma::sort_index(vec, "descend");
    idx = non_na.elem(sorted_idx.head(psi_));
  }
  
  // Create field to return
  arma::field<arma::uvec> ret(2);
  ret(0) = idx;
  ret(1) = arma::uvec(1).fill(non_na_ctr);
  
  // Return
  return ret;
}

// Function III - Compute Aggregated Predictive Distribution
arma::field<double> dsc_agg_density_par_(const arma::rowvec& active_weights, 
                                         const arma::rowvec& forecast_tvc_t,
                                         const arma::rowvec& variance_tvc_t, 
                                         const arma::uvec& idx_sub) {
  
  // Define Variables
  arma::field<double> ret(2);
  
  // Subset Matrices (select only active models)
  const arma::rowvec oos_forecast_tvp_sub = arma::conv_to<arma::rowvec>::from(forecast_tvc_t(idx_sub));
  const arma::rowvec oos_variance_tvp_sub = arma::conv_to<arma::rowvec>::from(variance_tvc_t(idx_sub));
  
  // Calculate Combined Predictive Density (Logarithmic Combination Rule)
  const double variance_comb = 1.0 / accu(active_weights / oos_variance_tvp_sub);
  const double       mu_comb = accu(active_weights % oos_forecast_tvp_sub / oos_variance_tvp_sub) *
    variance_comb;
  
  // Fill list  
  ret(0) = mu_comb;
  ret(1) = variance_comb;
  
  // Return list
  return ret;
}


// Function IV - Calculate (exponentially down-weighted) Log-Predictive-Likelihoods for Candidate Forecasts
void dsc_dpll_cands_par_(arma::rowvec& dpll_cands_gamma, 
                         double y_t, 
                         const arma::rowvec& forecast_tvc_t,             
                         const arma::rowvec& variance_tvc_t,
                         double gamma,
                         int method,
                         double risk_aversion,
                         double min_weight,
                         double max_weight) {
  
  // Define Variables
  const int n_cands = dpll_cands_gamma.n_elem;
  arma::rowvec performance_score(n_cands); performance_score.fill(arma::datum::nan);
  
  // Calculate Performance
  for (unsigned int i=0; i<n_cands; i++) {
    
    // Check for NA value
    if (arma::is_finite(forecast_tvc_t(i))) {  
      
      // Optimization-Method
      switch (method) {
      case 1: { // Log-Likelihood
    performance_score(i) = arma::log_normpdf(y_t,   
                      forecast_tvc_t(i),
                      pow(variance_tvc_t(i), 0.5));
    break;
  }
      case 2: { // MSE
        performance_score(i) = -pow(y_t - forecast_tvc_t(i), 2.0);
        break;
      }
      case 3: { // AE
        performance_score(i) = -std::abs(y_t - forecast_tvc_t(i));
        break;
      }
      case 4: { // Compounded Returns
        // Calculate Market Weight
        double w = (1.0 / risk_aversion) * (forecast_tvc_t(i) / variance_tvc_t(i));
        
        // Restrict Market Weight
        double weight = std::min(std::max(w, min_weight), max_weight);
        
        // Returns
        if (weight * y_t <= -1.0) {
          performance_score(i) = -10000;
        } else {
          performance_score(i) = log(1 + weight * y_t);
        }
        break;
      }
      default:
        throw std::invalid_argument("Error: Method not available");
      }
    }
  }
  
  // Exponentially down-weighted past predictive log-likelihoods
  dpll_cands_gamma = dpll_cands_gamma * gamma + performance_score * gamma;
}

// Function V - Rank and Select Aggregated Forecast
arma::field<double> rank_comb_par_(const arma::rowvec& dpll_combs, 
                                   const arma::rowvec& mu_comb_vec,
                                   const arma::rowvec& variance_comb_vec) {
  
  // Define Variables
  arma::field<double> ret(3);
  
  // Get Index of highest Rank
  const arma::uword best_idx = index_max(dpll_combs);
  
  // Select STSC-Forecast
  auto it_mu  = mu_comb_vec.begin() + best_idx;
  auto it_var = variance_comb_vec.begin() + best_idx;
  const double forecast_stsc = *it_mu;
  const double variance_stsc = *it_var;
  
  // Fill list  
  ret(0) = forecast_stsc;
  ret(1) = variance_stsc;
  ret(2) = best_idx;
  
  // Return list
  return ret;
}


// Function VI - Calculate (exponentially down-weighted) Log-Predictive-Likelihoods for Combinations (Aggregated Predictive Distributions)
void dsc_dpll_comb_par_(arma::rowvec& dpll_combs, 
                        double y_t, 
                        const arma::rowvec& forecasts_comb,             
                        const arma::rowvec& variances_comb,
                        double delta,
                        int method,
                        double risk_aversion,
                        double min_weight,
                        double max_weight) {
  
  // Define Variables
  const int n_combs = dpll_combs.n_cols;
  arma::rowvec performance_score(n_combs);
  
  // Calculate Performance
  for (unsigned int i=0; i<n_combs; i++) {
    
    // Optimization-Method
    switch(method) {
    case 1: { // Log-Likelihood
    performance_score(i) = arma::log_normpdf(y_t,   
                      forecasts_comb(i),
                      pow(variances_comb(i), 0.5));
    break;
  }
    case 2: { // SE
      performance_score(i) = -pow(y_t - forecasts_comb(i), 2.0);
      break;
    }
    case 3: { // AE
      performance_score(i) = -std::abs(y_t - forecasts_comb(i));
      break;
    }
    case 4: { // Compounded Returns
      // Calculate Market Weight
      double w = (1.0 / risk_aversion) * (forecasts_comb(i) / variances_comb(i));
      
      // Restrict Market Weight
      double weight = std::min(std::max(w, min_weight), max_weight); 
      
      // Returns
      if (weight * y_t <= -1.0) {
        performance_score(i) = -10000;
      } else {
        performance_score(i) = log(1 + weight * y_t);
      }
      break;
    }
    default:
      throw std::invalid_argument("Error: Method not available");
    }
  }
  
  // Exponentially down-weighted Past Predictive Log-Likelihoods
  dpll_combs = delta * dpll_combs + performance_score * delta;
}


// Function VII - Loop over Gamma and Psi
arma::field<arma::rowvec> dsc_loop_par_(arma::field<arma::rowvec>& dpll_cands,
                                        arma::rowvec& dpll_combs, 
                                        arma::rowvec gamma_grid, 
                                        arma::irowvec psi_grid, 
                                        double y_t,
                                        const arma::rowvec& forecast_tvc_t,             
                                        const arma::rowvec& variance_tvc_t,
                                        double delta,
                                        int method,
                                        bool equal_weight,
                                        double risk_aversion,
                                        double min_weight,
                                        double max_weight,
                                        int n_threads) { 
  
  // Define Variables
  const int n_combs = dpll_combs.n_cols;
  arma::rowvec forecasts_comb(n_combs); forecasts_comb.fill(arma::datum::nan);
  arma::rowvec variances_comb(n_combs); variances_comb.fill(arma::datum::nan);
  arma::field<arma::rowvec> chosen_cands(n_combs); 
  
  arma::field<double> stsc_results(3);
  arma::field<arma::rowvec> ret(4); 
  
  // ctr cube
  int ctr = 0;
  arma::mat ctr_mat(gamma_grid.n_elem, psi_grid.n_elem);
  for (unsigned int g=0; g<gamma_grid.n_elem; g++) {
    for (unsigned int p=0; p<psi_grid.n_elem; p++) {
      ctr_mat(g, p) = ctr;
      ctr++;
    }
  }
  
  // write comment
  const int psi_max = max(psi_grid);
  
  // for (unsigned int g=0; g<gamma_grid.n_elem; g++) {
  // Parallel Loop over gammas
  parallelFor(0, gamma_grid.n_elem, [&dpll_cands, &psi_max, &psi_grid, &ctr_mat,
                                     &gamma_grid, &chosen_cands, &equal_weight,
                                     &forecasts_comb, &variances_comb,
                                     &forecast_tvc_t, &variance_tvc_t,
                                     &y_t, &method, &risk_aversion, 
                                     &min_weight, &max_weight] (unsigned int g) {
    
    // Set Gamma
    const double gamma = gamma_grid(g);
    unsigned int ctr; 
                                       
    // Define variables for local scope                                    
    arma::field<arma::uvec> active_models(2);
    arma::field<double> agg_density(2);
    
    // Compute Active Set of Candidate Models
    active_models = dsc_active_models_par_(dpll_cands(g), psi_max);
    const arma::uvec active_idx = active_models(0);
    const int non_na_ctr = active_models(1)(0);
    
    // Loop over Psi
    for (unsigned int p=0; p<psi_grid.n_elem; p++) {
      
      // Get ctr
        ctr = ctr_mat(g, p);
      
        // Set Psi
        const int psi = std::min(non_na_ctr, psi_grid(p));
        
        // Get the Active Set of Candidate Models
        const arma::uvec seq_psi = arma::regspace<arma::uvec>(0, psi - 1);
        const arma::uvec active_idx_uvec = arma::conv_to<arma::uvec>::from(active_idx.elem(seq_psi));
        chosen_cands(ctr) = arma::conv_to<arma::rowvec>::from(active_idx_uvec);
        arma::rowvec active_weights(psi);
        
        // Set Active Weights
        if (equal_weight) {
          active_weights.fill(1.0 / psi);
        } else {
          const arma::rowvec raw_weights = arma::conv_to<arma::rowvec>::from(dpll_cands(g).elem(active_idx.elem(seq_psi))); 
          const arma::rowvec exp_raw_weights = exp(raw_weights);
          active_weights = exp_raw_weights / accu(exp_raw_weights); 
        }
        
        // Calculate Aggregated Predictive Density
        agg_density = dsc_agg_density_par_(active_weights,
                                           forecast_tvc_t,
                                           variance_tvc_t,
                                           active_idx_uvec); 
        
        // Assign Results
        forecasts_comb(ctr) = agg_density(0);
        variances_comb(ctr) = agg_density(1);
        //  ctr++;
    }
    
    // Update DPLL for Candidate Models
    dsc_dpll_cands_par_(dpll_cands(g), 
                        y_t, 
                        forecast_tvc_t,             
                        variance_tvc_t,
                        gamma,
                        method,
                        risk_aversion,
                        min_weight,
                        max_weight);
  }, n_threads); 
  
  // Select Aggregated Forecast
  stsc_results = rank_comb_par_(dpll_combs, 
                                forecasts_comb,
                                variances_comb);
  
  // Assign Results
  const double stsc_forecast = stsc_results(0);
  const double stsc_variance = stsc_results(1);
  const int stsc_idx = stsc_results(2);           
  
  // Update DPLL for Combinations (Aggregated Predictive Distributions)
  dsc_dpll_comb_par_(dpll_combs, 
                     y_t, 
                     forecasts_comb,             
                     variances_comb,
                     delta,
                     method,
                     risk_aversion,
                     min_weight,
                     max_weight);
  
  // Fill field   
  ret(0) = stsc_forecast;
  ret(1) = stsc_variance;
  ret(2) = stsc_idx;
  ret(3) = chosen_cands(stsc_idx); 
  
  // Return list
  return ret;
}

// ###################################################################################
// ###################################################################################
// ###################################################################################

// Helper Function -- Replication of Rcpp::setdiff()
arma::uvec my_setdiff_par(arma::uvec y, arma::uvec x) {
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  arma::uvec diff(x.n_elem);
  arma::uvec::iterator it = std::set_difference(x.begin(), x.end(), y.begin(), y.end(), diff.begin());
  diff.resize(it - diff.begin());
  return diff;
}

// 3.) STSC 
// Function I - Loop over t
// [[Rcpp::export]]
List stsc_loop_par_(const arma::vec& y, 
                    Nullable<const NumericMatrix&> X_, 
                    Nullable<const NumericMatrix&> Ext_F_, 
                    int sample_length,
                    arma::vec lambda_grid,
                    arma::vec kappa_grid,
                    int burn_in_tvc,
                    arma::rowvec gamma_grid, 
                    arma::irowvec psi_grid,
                    double delta,
                    int burn_in_dsc,
                    int method,
                    bool equal_weight,
                    int n_threads,
                    Nullable<NumericVector> risk_aversion_ = R_NilValue,
                    Nullable<NumericVector> min_weight_ = R_NilValue,
                    Nullable<NumericVector> max_weight_ = R_NilValue) { 
  
  
  // Check whether Simple Signals and / or Point Forecasts are provided and create combined Signal-Matrix
  arma::mat S;
  int n_raw_sig = 0, n_point_f = 0;
  
  if (X_.isNull() && Ext_F_.isNull()) {
    throw std::invalid_argument("Error: No signals provided");
  } else if (X_.isNotNull() && Ext_F_.isNotNull()) {
    n_raw_sig = as<arma::mat>(X_.get()).n_cols;
    n_point_f = as<arma::mat>(Ext_F_.get()).n_cols;
    S = arma::join_rows(as<arma::mat>(X_.get()), as<arma::mat>(Ext_F_.get()));
  } else if (Ext_F_.isNull()) {
    n_raw_sig = as<arma::mat>(X_.get()).n_cols;
    S = as<arma::mat>(X_.get());
  } else {
    n_point_f = as<arma::mat>(Ext_F_.get()).n_cols;
    S = as<arma::mat>(Ext_F_.get());
  }
  
  // Check Nullable Objects for method 4
  if (method == 4 && (risk_aversion_.isNull() || min_weight_.isNull() || max_weight_.isNull())) {
    throw std::invalid_argument("Error: Relevant parameter not provided!");
  }
  
  // Cast Nullable Objects
  double risk_aversion = arma::datum::nan;
  double min_weight    = arma::datum::nan;
  double max_weight    = arma::datum::nan;
  if (method == 4) {
    // Cast to double
    risk_aversion = as<double>(risk_aversion_.get());
    min_weight    = as<double>(min_weight_.get());
    max_weight    = as<double>(max_weight_.get());
  }
  
  //// Check Parameter
  //if (equal_weight != true && equal_weight != false) {
  //  throw std::invalid_argument("Error: Equal Weight argument wrong");
  //}
  
  // Number of Candiate Models and Signals
  const int n_signal = n_raw_sig + n_point_f;
  const int n_cands  = n_signal * lambda_grid.n_elem * kappa_grid.n_elem;
  
  // Define Variables for TV-C-Models
  const int tlength = y.n_elem;
  arma::rowvec forecast_tvc_pred(n_cands); forecast_tvc_pred.fill(arma::datum::nan);
  arma::rowvec variance_tvc_pred(n_cands); variance_tvc_pred.fill(arma::datum::nan);
  arma::cube   theta_cube(2, 1, n_cands), cov_mat_cube(2, 2, n_cands);
  arma::rowvec h_vec(n_cands);
  arma::uvec   current_na_cm, new_na_cm;
    
  // Define Variables for Dynamic Subset Combinations
  const int n_combs = gamma_grid.n_elem * psi_grid.n_elem;
  arma::vec  stsc_forecast(tlength); stsc_forecast.fill(arma::datum::nan);
  arma::vec  stsc_variance(tlength); stsc_variance.fill(arma::datum::nan);
  arma::vec stsc_idx(tlength); stsc_idx.fill(arma::datum::nan);  
  arma::field<arma::rowvec> dpll_cands(gamma_grid.n_elem);
  arma::rowvec dpll_combs(n_combs);
  List chosen_cands(tlength);
  
  // --- 
  // Apply TV-C-Init-Function
  List init_tvc_results(1);
  init_tvc_results = init_tvc_par_(y, 
                                   S,
                                   n_raw_sig,
                                   sample_length,
                                   lambda_grid,
                                   kappa_grid,
                                   n_threads);
  
  // Assign Results
  theta_cube    = as<arma::cube>(init_tvc_results(0));
  cov_mat_cube  = as<arma::cube>(init_tvc_results(1));
  h_vec         = as<arma::rowvec>(init_tvc_results(2));
  current_na_cm = as<arma::uvec>(init_tvc_results(3));
  
  // ---
  // Apply DSC-Init-Function
  arma::field<arma::field<arma::rowvec>> init_dsc_results(2);
  init_dsc_results = dsc_init_par_(n_cands,
                                   n_combs,
                                   gamma_grid.n_elem,
                                   current_na_cm); 
  
  // Assign Results
  dpll_cands = init_dsc_results(0);
  dpll_combs = init_dsc_results(1)(0);
  
  // --- 
  // Loop over t
  for (unsigned int t=0; t<(tlength-1); t++ ) {
    
    // Subset Data
    const double           y_t = y[t];
    const double        y_pred = y[t+1];
    const arma::rowvec     s_t = S.row(t);
    const arma::rowvec  s_pred = S.row(t+1);
    
    // Check for NA-Values in Candidate Models in t 
    new_na_cm.clear();
    int ctr = 0;
    for (unsigned int l = 0; l < lambda_grid.n_elem; l++) {
      for (unsigned int k = 0; k < kappa_grid.n_elem; k++) {
        for (unsigned int j = 0; j < S.n_cols; j++) {
          
          // Check and Count for NA-Values
          if (!arma::is_finite(s_t(j))) {
            new_na_cm.insert_rows(new_na_cm.n_rows, 1);
            new_na_cm(new_na_cm.n_rows-1) = ctr;                        
          }
          ctr++;
        }
      } 
    }
    
    // Identify Candidate Models that went to Non-Na
    if (new_na_cm.n_elem < current_na_cm.n_elem) {
      // Get the Index for the Signals that are not NA anymore
      arma::uvec vec_diff = my_setdiff_par(new_na_cm, current_na_cm);
      current_na_cm = new_na_cm;    
      for (unsigned int g=0; g<gamma_grid.n_elem; g++ ) {
        for (auto i : vec_diff) {
          dpll_cands(g)(i) = 0.0;
        }  
      }
    }
    
    // Check for Burn-In-Period
    if (t == (burn_in_tvc-1)) {
      
      arma::field<arma::field<arma::rowvec>> init_dsc_results_after_burn_in = dsc_init_par_(n_cands, n_combs, gamma_grid.n_elem, current_na_cm);
      dpll_cands = init_dsc_results_after_burn_in(0); 
      dpll_combs = init_dsc_results_after_burn_in(1)(0); 
      stsc_forecast.fill(arma::datum::nan); 
      stsc_variance.fill(arma::datum::nan); 
      stsc_idx.fill(arma::datum::nan);   
      IntegerVector idx = Range(0, t);
      chosen_cands[idx] = NA_INTEGER;
    }
    
    if (t == (burn_in_dsc-1)) {
      
      arma::field<arma::field<arma::rowvec>> init_dsc_results_after_burn_in = dsc_init_par_(n_cands, n_combs, gamma_grid.n_elem, current_na_cm);
      dpll_combs = init_dsc_results_after_burn_in(1)(0); 
      stsc_forecast.fill(arma::datum::nan); 
      stsc_variance.fill(arma::datum::nan); 
      stsc_idx.fill(arma::datum::nan);   
      IntegerVector idx = Range(0, t);
      chosen_cands[idx] = NA_INTEGER;
    }
    
    // Apply TV-C-Function
    arma::field<arma::rowvec> tvc_model_cand_results(2);
    tvc_model_cand_results = tvc_model_cand_par_(y_t, 
                                                 s_t,
                                                 s_pred, 
                                                 lambda_grid, 
                                                 kappa_grid, 
                                                 theta_cube, 
                                                 cov_mat_cube, 
                                                 h_vec,
                                                 n_threads); 
    
    // Assign Results
    forecast_tvc_pred = tvc_model_cand_results(0); 
    variance_tvc_pred = tvc_model_cand_results(1); 
        
    // Apply DSC-Function
    arma::field<arma::rowvec> dsc_results(4);
    dsc_results = dsc_loop_par_(dpll_cands,
                                dpll_combs, 
                                gamma_grid, 
                                psi_grid, 
                                y_pred,
                                forecast_tvc_pred,
                                variance_tvc_pred,
                                delta,
                                method,
                                equal_weight,
                                risk_aversion,
                                min_weight,
                                max_weight,
                                n_threads); 
    
    // Assign Results
    stsc_forecast(t+1) = dsc_results(0)(0);
    stsc_variance(t+1) = dsc_results(1)(0);
    stsc_idx(t+1)      = dsc_results(2)(0);
    chosen_cands(t+1)  = dsc_results(3);
  }
  
  // Fill list    
  List ret(4);
  ret[0] = stsc_forecast;
  ret[1] = stsc_variance;
  ret[2] = stsc_idx;
  ret[3] = chosen_cands;
  
  // Return list
  return ret;
}
