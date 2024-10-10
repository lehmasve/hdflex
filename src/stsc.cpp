#include <RcppArmadillo.h>
using namespace Rcpp;

// 1) TVC Regressions
// Function I - Initialize TVC-Parameter
   List init_tvc_(const arma::vec& y, 
                  const arma::mat& S,
                  int n_raw_sig,
                  int init,
                  arma::vec lambda_grid,
                  arma::vec kappa_grid,
                  bool bias) {

   // Get Dimensions   
      int n_signal = S.n_cols;
      int n_cands  = n_signal * lambda_grid.n_elem * kappa_grid.n_elem;

   // Define Variables 
      arma::cube theta_cube(2, 1, n_cands);
      arma::cube cov_mat_cube(2, 2, n_cands);
      arma::rowvec h_vec(n_cands);
      arma::vec y_sample, x, cm_na_idx(n_cands); cm_na_idx.fill(arma::datum::nan);
      arma::uvec non_finite, init_idx;
      arma::mat x_sample_one, theta, cov_mat;
      arma::colvec coef;
      double intercept, var_y, var_x, h;
      List ret_all(4);  

   // Loop over all candidates models
      int ctr = 0;
         for (unsigned int l = 0; l < lambda_grid.n_elem; l++) {
           for (unsigned int k = 0; k < kappa_grid.n_elem; k++) {
               for (int j = 0; j < n_signal; j++) {

               // Select Signal
                  x = S.col(j);

               // Check and Count NA-Values
                  non_finite = arma::find_nonfinite(x);
                  int na_ctr = non_finite.n_elem;
                  if (na_ctr > 0) {
                     cm_na_idx(ctr) = ctr;
                  }
   
               // Index for subsetting the initialisation sample
                  init_idx = arma::regspace<arma::uvec>(0 + na_ctr, na_ctr + init - 1);
   
               // Define and prepare matrices for regression  
                  y_sample = y.elem(init_idx);  
                  arma::mat x_sample(x.elem(init_idx));   
                  x_sample_one = arma::ones<arma::mat>(x_sample.n_rows, 1);
                  x_sample.insert_cols(0, x_sample_one);

               // Initialize - Theta
                  theta = arma::zeros<arma::mat>(2,1);

               // Initialize - System Covariance
                  coef = solve(x_sample, y_sample);
                  intercept = coef(0);
                  var_y = arma::var(y_sample);         
                  var_x = arma::var(x_sample.col(1));  
                  cov_mat = arma::zeros<arma::mat>(2, 2); 
               
               // Set Intercept Variance
                  cov_mat(0, 0) = pow(intercept, 2) + var_y; // Set to Zero for Constant Intercept

               // Distinguish between P- and F-Signals
                  if(j < n_raw_sig) {
                    cov_mat(1, 1) = (var_x != 0.0) ? var_y / var_x : var_y;
                  } else {
                     theta(1, 0) = 1.0;   // -> Slope Coefficient 1.0
                     cov_mat(1, 1) = 0.0; // -> Constant Slope Coefficient
                     if(!bias) {
                        cov_mat(0, 0) = 0.0; // -> Constant Intercept
                     }
                  }

               // Initialize - Observational Variance
                  h = var_y;

               // Fill Cubes  
                  theta_cube.slice(ctr) = theta;
                  cov_mat_cube.slice(ctr) = cov_mat; 
                  h_vec(ctr) = h;
                  ctr++;
          }
        }
      }
   
      // Fill Return-List  
         ret_all[0] = theta_cube;
         ret_all[1] = cov_mat_cube; 
         ret_all[2] = h_vec;
         ret_all[3] = arma::conv_to<arma::uvec>::from(cm_na_idx.elem(arma::find_finite(cm_na_idx)));
     
      // Return Results
         return ret_all;
   }
// ----------


// 1) TVC Regressions
// Function II - Predictive and Update Step for Signal j and time t
   arma::field<double> tvc_model_(double y_t, 
                                  double s_t_j, 
                                  double s_pred_j,
                                  double lambda, 
                                  double kappa, 
                                  arma::mat& theta, 
                                  arma::mat& cov_mat, 
                                  double& h) { 
  
   // Define Variables
      arma::field<double> ret(2);

   // Get Signal for time t and t + 1
      const arma::mat z_t = {1.0, s_t_j};
      const arma::mat z_pred = {1.0, s_pred_j};
   
   // Add noise to uncertainty of coefficients in time t (see Equation 5)
      const arma::mat r_upt = cov_mat / lambda;

   // Calculate (OOS) Forecast Error for time t (see Equation 7)
      const double e_t = arma::as_scalar(y_t - z_t * theta);

   // Update Observational Variance in time t (see Equation 10 and 11)
      h = arma::as_scalar(kappa * h + (1 - kappa) * pow(e_t, 2));
    
   // Update Coefficients in time t (see Equation 7)
      const double inv_tvar = arma::as_scalar(1.0 / (h + z_t * r_upt * z_t.t()));
      theta = theta + r_upt * z_t.t() * inv_tvar * e_t;
     
   // Update Uncertainty of Coefficients in time t (see Equation 8)
      cov_mat = r_upt - r_upt * z_t.t() * inv_tvar * (z_t * r_upt);

   // Get Predictive Density for Predicting t + 1 (see Equation 9)
      const double mu = arma::as_scalar(z_pred * theta);
      const double variance = arma::as_scalar(h + z_pred * ((1.0 / lambda) * cov_mat) * z_pred.t());

   // Fill Return-Field  
      ret(0) = mu;
      ret(1) = variance;
   
   // Return  
      return ret;
}


// 1) TVC Regressions
// Function III - Loop over Signals, Lambda and Kappa
   arma::field<arma::rowvec> tvc_model_cand_(double y_t, 
                                             const arma::rowvec& s_t,
                                             const arma::rowvec& s_pred, 
                                             arma::vec lambda_grid, 
                                             arma::vec kappa_grid, 
                                             arma::cube& theta_cube, 
                                             arma::cube& cov_mat_cube, 
                                             arma::rowvec& h_vec) {      
  
   // Get Dimensions
      int n_cands = s_t.n_elem * lambda_grid.n_elem * kappa_grid.n_elem;

   // Define Variables
      arma::field<arma::rowvec> ret(2);
      arma::rowvec mu_vec(n_cands);
      arma::rowvec variance_vec(n_cands);
    
   // Loop over Lambda
      int counter = 0;
      for (unsigned int l = 0; l < lambda_grid.n_elem; l++) {

      // Set Lambda
         double lambda = lambda_grid(l);

      // Loop over Kappa
         for (unsigned int k = 0; k < kappa_grid.n_elem; k++) {

         // Set Kappa
            double kappa =  kappa_grid(k);

         // Loop over all candidates
            for (unsigned int j = 0; j < s_t.n_elem; j++) {
        
            // Set Signals
               const double    s_t_j = s_t(j);
               const double s_pred_j = s_pred(j);

            // Check if signal is NA or not 
               bool is_na = !arma::is_finite(s_t_j);
               if(!is_na) {
   
               // Apply TVC-Function
                  const arma::field<double> tvc_results = tvc_model_(y_t,
                                                                     s_t_j,
                                                                     s_pred_j, 
                                                                     lambda,
                                                                     kappa,
                                                                     theta_cube.slice(counter), 
                                                                     cov_mat_cube.slice(counter), 
                                                                     h_vec(counter)); 
               // Assign TVC-Model-Results 
                  mu_vec(counter)       = tvc_results(0); 
                  variance_vec(counter) = tvc_results(1);

               } else {

               // Assign TVC-Model-Results 
                  mu_vec(counter)       = arma::datum::nan; 
                  variance_vec(counter) = arma::datum::nan;
               }

               // Update Counter
                  counter++;
            }
         }
      }
      
      // Fill Field 
         ret(0) = mu_vec;
         ret(1) = variance_vec;
     
      // Return
         return ret;
}

// ###################################################################################
// ###################################################################################
// ###################################################################################

// 2) Dynamic Subset Combination
// Function I - Initialize DSC-Parameter
   arma::field<arma::field<arma::rowvec>> dsc_init_(int n_cands,
                                                    int n_combs,
                                                    int n_gamma,
                                                    arma::uvec na_idx) {

   // Define Variables
      arma::field<arma::field<arma::rowvec>> ret(2); 

   // Initialize Vector for Performance-Score (Subset Combinations) -> Ranking
      arma::rowvec score_combs(n_combs, arma::fill::zeros);
   
   // Initialize Vector for Performance-Score (Candidate Models) -> Ranking 
      arma::rowvec vec(n_cands, arma::fill::zeros);
                   vec.elem(na_idx).fill(arma::datum::nan);

   // Fill Field for Candidate Models
      arma::field<arma::rowvec> score_cands(n_gamma);
      for (int i = 0; i < n_gamma; ++i) {
          score_cands(i) = vec;
      }
  
   // Fill Return-Field 
      ret(0) = score_cands;
      ret(1) = arma::field<arma::rowvec>(1);
      ret(1)(0) = score_combs;

   // Return
      return ret;
}


// 2) Dynamic Subset Combination
// Function II - Rank and Set Active Model Set (Active Models)
   arma::field<arma::uvec> dsc_active_models_(const arma::rowvec& score_cands_gamma,  
                                              int psi) {

   // Define Variables
      arma::field<arma::uvec> ret(2);

   // Get indices of Non-NaN values
      const arma::uvec non_na = arma::find_finite(score_cands_gamma);
      const int non_na_ctr = non_na.n_elem;
      const int psi_ = std::min(non_na_ctr, psi);

   // Check if all elements are equal (without NAs)
      const arma::vec vec = score_cands_gamma.elem(non_na);
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

   // Fill Return-Field
      ret(0) = idx;
      ret(1) = arma::uvec(1).fill(non_na_ctr);

   // Return
      return ret;
}


// 2) Dynamic Subset Combination
// Function III - Compute Aggregated Predictive Distribution
   arma::field<double> dsc_agg_density_(const arma::rowvec& active_weights, 
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
      const double mu_comb = accu(active_weights % oos_forecast_tvp_sub / oos_variance_tvp_sub) *
                             variance_comb;
   
   // Fill Field  
      ret(0) = mu_comb;
      ret(1) = variance_comb;

   // Return
      return ret;
}


// 2) Dynamic Subset Combination
// Function IV - Calculate (exponentially down-weighted) Performance Scores (-> Ranking) for Candidate Forecasts
   void dsc_score_cands_(arma::rowvec& score_cands_gamma, 
                         double y_t, 
                         const arma::rowvec& forecast_tvc_t,             
                         const arma::rowvec& variance_tvc_t,
                         double gamma,
                         int metric,
                         double risk_aversion,
                         double min_weight,
                         double max_weight) {
                                     
   // Define Variables
      int n_cands = score_cands_gamma.n_elem;
      arma::rowvec performance_score(n_cands); performance_score.fill(arma::datum::nan);

   // Calculate Performance
      for (int i=0; i<n_cands; i++) {

      // Check for NA value
         if (arma::is_finite(forecast_tvc_t(i))) {

            // Optimization-metric
               switch (metric) {
                  case 1: {
                     // Predictive-Log-Likelihoods
                     performance_score(i) = arma::log_normpdf(y_t,   
                                                              forecast_tvc_t(i),
                                                              pow(variance_tvc_t(i), 0.5));
                     break;
                  }
                  case 2: {
                     // Squared-Errors
                     performance_score(i) = -pow(y_t - forecast_tvc_t(i), 2.0);
                     break;
                  }
                  case 3: {
                     // Absolute-Errors
                     performance_score(i) = -std::abs(y_t - forecast_tvc_t(i));
                     break;
                  }
                  case 4: {
                     // Compounded Returns
                     // Calculate Market Weight
                        double w = (1.0 / risk_aversion) * (forecast_tvc_t(i) / variance_tvc_t(i));
   
                     // Restrict Market Weight
                        double weight = std::min(std::max(w, min_weight), max_weight);
   
                     // Returns
                        if (weight * y_t <= -1.0) {
                           performance_score(i) = -10000;
                        } else {
                           performance_score(i) = log(1.0 + weight * y_t);
                        }
                        break;
                  }
                  case 5: {
                     // Continuous-Ranked-Probability-Scores
                     // Convert
                        double obs = y_t;
                        double mu = forecast_tvc_t(i);
                        double sig = pow(variance_tvc_t(i), 0.5);
   
                     // Standardize Prediction Error
                        double z = (obs - mu) / sig;
   
                     // PDF evaluated at normalized Prediction Error
                        double pdf = arma::normpdf(z);
   
                     // CDF evaluated at normalized Prediction Error
                        double cdf = arma::normcdf(z);
   
                     // Inverse of pi
                        double pi_inv = 1.0 / pow(arma::datum::pi, 0.5);
   
                     // Compute Continuous Ranked Probability Score
                        double crps = sig * (z * (2.0 * cdf - 1.0) + 2.0 * pdf - pi_inv);
   
                     // Assign
                        performance_score(i) = -crps;
                        break;
                  }
                  default:
                     throw std::invalid_argument("Error: Metric not available");
               }
         }
      }

   // Exponentially down-weighted past performance scores
      score_cands_gamma = score_cands_gamma * gamma + performance_score * gamma;
}


// 2) Dynamic Subset Combination
// Function V - Rank and Select Aggregated Forecast
   arma::field<double> rank_comb_(const arma::rowvec& score_combs, 
                                  const arma::rowvec& mu_comb_vec,
                                  const arma::rowvec& variance_comb_vec) {
                              
   // Define Variables
      arma::field<double> ret(3);
   
   // Get index of highest performance score
      arma::uword best_idx = index_max(score_combs);
   
   // Select STSC-Forecast
      double forecast_stsc = mu_comb_vec(best_idx);
      double variance_stsc = variance_comb_vec(best_idx);
   
   // Fill Field  
      ret(0) = forecast_stsc;
      ret(1) = variance_stsc;
      ret(2) = best_idx;

   // Return
      return ret;
   }


// 2) Dynamic Subset Combination
// Function VI - Calculate (exponentially down-weighted) performance scores (-> ranking) for combinations (Aggregated Predictive Distributions)
   void dsc_score_comb_(arma::rowvec& score_combs, 
                       double y_t, 
                       const arma::rowvec& forecasts_comb,             
                       const arma::rowvec& variances_comb,
                       double delta,
                       int metric,
                       double risk_aversion,
                       double min_weight,
                       double max_weight) {
                                     
   // Define Variables
      int n_combs = score_combs.n_cols;
      arma::rowvec performance_score(n_combs);

   // Calculate Performance
      for (int i=0; i<n_combs; i++) {

      // Optimization-metric
         switch(metric) {
            case 1: { 
               // Predictive-Log-Likelihoods
               performance_score(i) = arma::log_normpdf(y_t,   
                                                        forecasts_comb(i),
                                                        pow(variances_comb(i), 0.5));
               break;
            }
            case 2: { 
               // Squared-Errors
               performance_score(i) = -pow(y_t - forecasts_comb(i), 2.0);
               break;
            }
            case 3: {
               // Absolute-Errors
               performance_score(i) = -std::abs(y_t - forecasts_comb(i));
               break;
            }
            case 4: {
               // Compounded Returns
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
            case 5: {
               // Continuous-Ranked-Probability-Scores
               // Convert
                  double obs = y_t;
                  double mu = forecasts_comb(i);
                  double sig = pow(variances_comb(i), 0.5);

               // Standardize Prediction Error
                  double z = (obs - mu) / sig;

               // PDF evaluated at normalized Prediction Error
                  double pdf = arma::normpdf(z);

               // CDF evaluated at normalized Prediction Error
                  double cdf = arma::normcdf(z);

               // Inverse of pi
                  double pi_inv = 1.0 / pow(arma::datum::pi, 0.5);

               // Compute Continuous Ranked Probability Score
                  double crps = sig * (z * (2.0 * cdf - 1.0) + 2.0 * pdf - pi_inv);

               // Assign
                  performance_score(i) = -crps;
                  break;
               }
            default:
               throw std::invalid_argument("Error: Metric not available");
            }
      }
    
   // Exponentially down-weighted past performance scores
      score_combs = delta * score_combs + performance_score * delta;
   }


// 2) Dynamic Subset Combination
// Helper Function to remove duplicates between arma::uvec and include_idx
   arma::uvec remove_duplicates_(const arma::uvec& vec,
                                 const arma::uvec& values) {

      // Define variables
         std::unordered_set<arma::uword> seen;
         std::unordered_set<arma::uword> target_values(values.begin(), values.end());
         std::vector<arma::uword> result;
 
      // Remove duplicates for specified values
         for (arma::uword i = 0; i < vec.n_elem; ++i) {
           arma::uword val = vec[i];
           if (target_values.find(val) != target_values.end()) {
             // Only check for duplicates for specified values
             if (seen.find(val) == seen.end()) {
               seen.insert(val);
               result.push_back(val);
             }
           } else {
             // Directly add other values
             result.push_back(val);
           }
         }
 
      // Convert std::vector to arma::uvec and return
         return arma::uvec(result);
   }


// 2) Dynamic Subset Combination
// Function VII - Loop over Gamma and Psi
   arma::field<arma::rowvec> dsc_loop_(arma::field<arma::rowvec>& score_cands,
                                       arma::rowvec& score_combs, 
                                       arma::rowvec gamma_grid, 
                                       arma::irowvec psi_grid, 
                                       double y_t,
                                       const arma::rowvec& forecast_tvc_t,             
                                       const arma::rowvec& variance_tvc_t,
                                       double delta,
                                       int metric,
                                       bool equal_weight,
                                       arma::uvec incl_idx,
                                       double risk_aversion,
                                       double min_weight,
                                       double max_weight) {

   // Define Variables
      int n_combs = score_combs.n_cols;
      arma::rowvec forecasts_comb(n_combs); forecasts_comb.fill(arma::datum::nan);
      arma::rowvec variances_comb(n_combs); variances_comb.fill(arma::datum::nan);
      arma::field<arma::rowvec> chosen_cands(n_combs); 
      arma::field<arma::uvec> active_models(2);
      arma::field<double> agg_density(2);
      arma::field<double> stsc_results(3);
      arma::field<arma::rowvec> ret(4); 

   // Set highest value for psi
      int psi_max = max(psi_grid);

   // Loop over Gamma and Psi
      int ctr = 0;
      for (unsigned int g=0; g<gamma_grid.n_elem; g++) {

         // Set Gamma
            double gamma = gamma_grid(g);

         // Compute Maximum Set of Active Candidate Models
            active_models = dsc_active_models_(score_cands(g), psi_max);
            arma::uvec active_idx = active_models(0);
            const int non_na_ctr = active_models(1)(0);

         // Add Keep-Index & Remove Duplicates
            if (incl_idx.n_elem > 0) {
               active_idx = arma::join_cols(incl_idx, active_idx); // -> Keep must be empty arma::uvec if no keep
               active_idx = remove_duplicates_(active_idx, incl_idx);
            }

         // Loop over Psi
            for (unsigned int p=0; p<psi_grid.n_elem; p++) {

            // Set Psi
               const int psi = std::min(non_na_ctr, psi_grid(p));
               const arma::uvec seq_psi = arma::regspace<arma::uvec>(0, psi - 1);

            // Select Active Set of Candidate Models
               const arma::uvec active_idx_uvec = arma::conv_to<arma::uvec>::from(active_idx.elem(seq_psi));

            // Save Active Set of Candidate Models
               chosen_cands(ctr) = arma::conv_to<arma::rowvec>::from(active_idx_uvec);

            // Create Active Weight Vector
               arma::rowvec active_weights(psi);
               if (equal_weight) {
                  active_weights.fill(1.0 / psi);
               } else {
                  const arma::rowvec raw_weights = arma::conv_to<arma::rowvec>::from(score_cands(g).elem(active_idx.elem(seq_psi))); 
                  const arma::rowvec exp_raw_weights = exp(raw_weights);
                  active_weights = exp_raw_weights / accu(exp_raw_weights); 
               }

            // Calculate Aggregated Predictive Density
               agg_density = dsc_agg_density_(active_weights,
                                              forecast_tvc_t,
                                              variance_tvc_t,
                                              active_idx_uvec); 
   
            // Assign Results
               forecasts_comb(ctr) = agg_density(0);
               variances_comb(ctr) = agg_density(1);
               ctr++;
            }

         // Update score for Candidate Models
            dsc_score_cands_(score_cands(g), 
                             y_t, 
                             forecast_tvc_t,             
                             variance_tvc_t,
                             gamma,
                             metric,
                             risk_aversion,
                             min_weight,
                             max_weight);
      }

   // Select Aggregated Forecast
      stsc_results = rank_comb_(score_combs, 
                                forecasts_comb,
                                variances_comb);

   // Assign Results
      double stsc_forecast = stsc_results(0);
      double stsc_variance = stsc_results(1);
      int stsc_idx = stsc_results(2);           

   // Update score for Combinations (Aggregated Predictive Distributions)
      dsc_score_comb_(score_combs, 
                      y_t, 
                      forecasts_comb,             
                      variances_comb,
                      delta,
                      metric,
                      risk_aversion,
                      min_weight,
                      max_weight);

   // Fill field   
      ret(0) = stsc_forecast;
      ret(1) = stsc_variance;
      ret(2) = stsc_idx;
      ret(3) = chosen_cands(stsc_idx); 

   // Return
      return ret;
}

// ###################################################################################
// ###################################################################################
// ###################################################################################

// Helper Function -- Replication of Rcpp::setdiff()
   arma::uvec my_setdiff_(arma::uvec y, arma::uvec x) {
       std::sort(x.begin(), x.end());
       std::sort(y.begin(), y.end());
       arma::uvec diff(x.n_elem);
       arma::uvec::iterator it = std::set_difference(x.begin(), x.end(), y.begin(), y.end(), diff.begin());
       diff.resize(it - diff.begin());
       return diff;
   }

// Helper Function -- Median with NAs
   double compute_median_(arma::rowvec vec) {
      
      // Filter out NA values
         arma::vec finiteVec = vec(arma::find_finite(vec));

      // Check if empty
         if (finiteVec.empty()) {
           return arma::datum::nan; // Return NA if there are no finite values
         }

      // Calculate and return median
         return arma::median(finiteVec);
}

// ###################################################################################
// ###################################################################################
// ###################################################################################


// 3.) STSC 
// Function I - Loop over t
// [[Rcpp::export]]
   List stsc_loop_(const arma::vec& y, 
                   Nullable<const NumericMatrix&> X_, 
                   Nullable<const NumericMatrix&> Ext_F_, 
                   int init,
                   arma::vec lambda_grid,
                   arma::vec kappa_grid,
                   bool bias,
                   arma::rowvec gamma_grid, 
                   arma::irowvec psi_grid, 
                   double delta,
                   int burn_in,
                   int burn_in_dsc,
                   int metric,
                   bool equal_weight,
                   Nullable<IntegerVector> incl_,
                   Nullable<NumericVector> portfolio_params_) {
   
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

   // Check Nullable Objects for metric 4
      if (metric == 4 && portfolio_params_.isNull()) {
         throw std::invalid_argument("Error: Relevant parameter not provided!");
      }
   
   // Cast Nullable Objects for metric 4
      double risk_aversion = arma::datum::nan;
      double min_weight = arma::datum::nan;
      double max_weight = arma::datum::nan;
      if (metric == 4) {
         // Cast to NumericVector and extract values
         NumericVector combined_params = as<NumericVector>(portfolio_params_.get());
         if (combined_params.size() != 3) {
            throw std::invalid_argument("Error: portfolio_params_ must contain exactly 3 elements!");
         }
         risk_aversion = combined_params[0];
         min_weight = combined_params[1];
         max_weight = combined_params[2];
      }
      
   // Number of Candiate Models and Signals
      int n_signal = n_raw_sig + n_point_f;
      int n_cands  = n_signal * lambda_grid.n_elem * kappa_grid.n_elem;
   
   // Define Variables for TVC-Models
      int tlength = y.n_elem;
      arma::rowvec forecast_tvc_pred(n_cands); forecast_tvc_pred.fill(arma::datum::nan);
      arma::rowvec variance_tvc_pred(n_cands); variance_tvc_pred.fill(arma::datum::nan);
      arma::cube   theta_cube(2, 1, n_cands), cov_mat_cube(2, 2, n_cands);
      arma::rowvec h_vec(n_cands);
      arma::uvec   current_na_cm, new_na_cm;

   // Define Variables for Dynamic Subset Combinations
      int n_combs = gamma_grid.n_elem * psi_grid.n_elem;
      arma::vec  stsc_forecast(tlength); stsc_forecast.fill(arma::datum::nan);
      arma::vec  stsc_variance(tlength); stsc_variance.fill(arma::datum::nan);
      arma::vec stsc_idx(tlength); stsc_idx.fill(arma::datum::nan);  
      arma::field<arma::rowvec> score_cands(gamma_grid.n_elem);
      arma::rowvec score_combs(n_combs);
      List chosen_cands(tlength);

   // Include: TVC-Models that must be included in the Subsets
      arma::uvec incl_idx;
      if (incl_.isNotNull()) {

      // Cast Values to uvec incl
         arma::uvec incl = as<arma::uvec>(incl_.get());

      // Number of TVC-Models per Signal
         int grid_size = lambda_grid.n_elem * kappa_grid.n_elem;

      // Resize incl_idx (-> Number of TVC-Models that must be included)
         incl_idx.set_size(grid_size * incl.n_elem);

      // Calculate the Indices of the TVC-Models
         int ctr = 0;
         for (arma::uword k = 0; k < incl.n_elem; ++k) {
           for (int i = 0; i < grid_size; ++i) {
             arma::uword index = (incl[k] - 1) + n_signal * i;
             incl_idx(ctr) = index;
             ctr++;
           }
         }
      }

   // --- 
   // Apply TVC-Init-Function
      List init_tvc_results(1);
      init_tvc_results = init_tvc_(y, 
                                   S,
                                   n_raw_sig,
                                   init,
                                   lambda_grid,
                                   kappa_grid,
                                   bias);

   // Assign Results
      theta_cube    = as<arma::cube>(init_tvc_results(0));
      cov_mat_cube  = as<arma::cube>(init_tvc_results(1));
      h_vec         = as<arma::rowvec>(init_tvc_results(2));
      current_na_cm = as<arma::uvec>(init_tvc_results(3));

   // ---
   // Apply DSC-Init-Function
      arma::field<arma::field<arma::rowvec>> init_dsc_results(2);
      init_dsc_results = dsc_init_(n_cands,
                                   n_combs,
                                   gamma_grid.n_elem,
                                   current_na_cm); 

   // Assign Results
      score_cands = init_dsc_results(0); 
      score_combs = init_dsc_results(1)(0); 

   // --- 
   // Loop over t
      for (int t=0; t<(tlength-1); t++ ) {

      // Subset Data
         const double y_t = y[t];
         const double y_pred = y[t+1];
         const arma::rowvec s_t = S.row(t);
         const arma::rowvec s_pred = S.row(t+1);

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
            arma::uvec vec_diff = my_setdiff_(new_na_cm, current_na_cm);
                  current_na_cm = new_na_cm;    
            for (unsigned int g=0; g<gamma_grid.n_elem; g++ ) {
               for (auto i : vec_diff) {
                  score_cands(g)(i) = compute_median_(arma::conv_to<arma::rowvec>::from(score_cands(g))); // 0.0;  // -> Insert Value !!!
               }  
            }
         }
   
      // Check for Burn-In-Period
         if (t == (burn_in-1)) {

            arma::field<arma::field<arma::rowvec>> init_dsc_results_after_burn_in = dsc_init_(n_cands, n_combs, gamma_grid.n_elem, current_na_cm);
            score_cands = init_dsc_results_after_burn_in(0); 
            score_combs = init_dsc_results_after_burn_in(1)(0); 
            stsc_forecast.fill(arma::datum::nan); 
            stsc_variance.fill(arma::datum::nan); 
            stsc_idx.fill(arma::datum::nan);   
            IntegerVector idx = Range(0, t);
            chosen_cands[idx] = NA_INTEGER;
         }
         
         if (t == (burn_in_dsc-1)) {

            arma::field<arma::field<arma::rowvec>> init_dsc_results_after_burn_in = dsc_init_(n_cands, n_combs, gamma_grid.n_elem, current_na_cm);
            score_combs = init_dsc_results_after_burn_in(1)(0); 
            stsc_forecast.fill(arma::datum::nan); 
            stsc_variance.fill(arma::datum::nan); 
            stsc_idx.fill(arma::datum::nan);   
            IntegerVector idx = Range(0, t);
            chosen_cands[idx] = NA_INTEGER;
         }

      // Apply TVC-Function
         arma::field<arma::rowvec> tvc_model_cand_results(2);
         tvc_model_cand_results = tvc_model_cand_(y_t, 
                                                  s_t,
                                                  s_pred, 
                                                  lambda_grid, 
                                                  kappa_grid, 
                                                  theta_cube, 
                                                  cov_mat_cube, 
                                                  h_vec); 

      // Assign Results
         forecast_tvc_pred = tvc_model_cand_results(0); 
         variance_tvc_pred = tvc_model_cand_results(1); 

      // Apply DSC-Function
         arma::field<arma::rowvec> dsc_results(4);
         dsc_results = dsc_loop_(score_cands,
                                 score_combs, 
                                 gamma_grid, 
                                 psi_grid, 
                                 y_pred,
                                 forecast_tvc_pred,
                                 variance_tvc_pred,
                                 delta,
                                 metric,
                                 equal_weight,
                                 incl_idx,
                                 risk_aversion,
                                 min_weight,
                                 max_weight); 

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
