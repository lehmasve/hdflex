#include <RcppArmadillo.h>
using namespace Rcpp;

// 1) TV-C Regressions
// Function I - Initialize TV-C-Parameter
// [[Rcpp::export]]
   List init_tvc_(const arma::vec& y, 
                  const arma::mat& S,
                  const int& n_sim_sig,
                  const int& sample_length,
                  const arma::vec& lambda_grid,
                  const arma::vec& kappa_grid) {

   // Define Variables      
      int n_signal = S.n_cols;
      int n_cands  = n_signal * lambda_grid.size() * kappa_grid.size();
      List theta_all(n_cands);
      List cov_mat_all(n_cands);
      List h_all(n_cands);
      IntegerVector na_idx;
   
   // Loop over all candidates
      int counter = 0;
         for (int l = 0; l < lambda_grid.size(); l++) {
           for (int k = 0; k < kappa_grid.size(); k++) {
               for (int j = 0; j < n_signal; j++) {
   
               // Define Variables
                  arma::mat x_sample_one, theta, cov_mat;
                  double intercept, var_y, var_x, h;
   
               // Select Signal
                  arma::vec x = S.col(j);

               // Check and Count for NA-Values
                  arma::uvec non_finite = arma::find_nonfinite(x);
                  int na_ctr = non_finite.size();
                  if (na_ctr > 0) {
                     na_idx.push_back(counter);
                  }
   
               // Index for subsetting
                  arma::uvec index = arma::regspace<arma::uvec>(0 + na_ctr, na_ctr + sample_length - 1);
   
               // Define and prepare matrices for regression  
                  arma::vec y_sample = y.elem(index);  
                  arma::mat x_sample(x.elem(index));   
                  x_sample_one = arma::ones<arma::mat>(x_sample.n_rows, 1); // muss nicht im Loop passieren
                  x_sample.insert_cols(0, x_sample_one);

               // Check if all elements are equal
                  if(arma::all(x_sample.col(1) == x_sample(0, 1))) {
                     Rcout << "Warning: Consider increasing sample_length - Column " << j + 1 << " only contains equal values. \n Initialization of corresponding TV-C-Model might be affected. \n";
                  }

               // While loop to differentiate between simple signals and point forecasts
                  if(j < n_sim_sig) {
   
                  // Initialize - Theta
                     theta = arma::zeros<arma::mat>(2,1); // muss nicht im Loop passieren

                  // Initialize - System Covariance
                     arma::colvec coef = solve(x_sample, y_sample);
                     intercept = coef(0);
                     var_y     = arma::var(y_sample);         
                     var_x     = arma::var(x_sample.col(1));  
                     cov_mat   = arma::zeros<arma::mat>(2, 2);  // muss nicht im Loop passieren
                     cov_mat(0, 0) = pow(intercept, 2) + var_y;
                     if(arma::is_finite(var_y / var_x)) {
                          cov_mat(1, 1) =  var_y / var_x; 
                     } else {
                          cov_mat(1, 1) = 0;
                     }
               
                  } else {

                  // Initialize - Theta
                     theta = arma::zeros<arma::mat>(2,1);  // muss nicht im Loop passieren
                     theta(1, 0) = 1;

                  // Initialize - System Covariance
                     arma::colvec coef = solve(x_sample, y_sample);
                     intercept = coef(0);
                     var_y     = arma::var(y_sample);         
                     var_x     = arma::var(x_sample.col(1));  
                     cov_mat   = arma::zeros<arma::mat>(2, 2);  // muss nicht im Loop passieren
                     cov_mat(0, 0) =  pow(intercept, 2) + var_y;  // Set to Zero for Constant Intercept
                     cov_mat(1, 1) =  0;
                  }

               // Initialize - Observational Variance
                  h = var_y;

               // Fill Return-List  
                  theta_all[counter]   = theta;
                  cov_mat_all[counter] = cov_mat; 
                  h_all[counter]       = h;
                  counter++;
          }
        }
      }
   
      // Fill Return-List  
         List ret_all(4);
         ret_all[0] = theta_all;
         ret_all[1] = cov_mat_all; 
         ret_all[2] = h_all;
         ret_all[3] = na_idx;
     
      // Return Results
         return ret_all;
   }
// ----------


// 1) TV-C Regressions
// Function II - Predictive and Update Step for Signal j and time t
// [[Rcpp::export]]
   List tvc_model_(const double& y_t, 
                   const double& s_t_j, 
                   const double& s_tplus1_j,
                   const double& lambda, 
                   const double& kappa, 
                   const arma::mat& theta, 
                   const arma::mat& cov_mat, 
                   const double& h) { 
  
   // Define Variables
      arma::mat z_t, z_pred, r_t, theta_new, cov_mat_new;
      double mu, variance, h_new;
  
   // Get Predictor for Time t
      z_t  = arma::zeros<arma::mat>(1, 2);  
      z_t(0, 0) = 1;
      z_t(0, 1) = s_t_j; 
   
   // Get Predictor for Predicting t + 1
      z_pred = arma::zeros<arma::mat>(1, 2);  
      z_pred(0, 0) = 1;
      z_pred(0, 1) = s_tplus1_j; 
   
   // Calculate R for Time t (Equation 5)
      //r_t = cov_mat / lambda;
      r_t = (1 - lambda) * cov_mat / lambda;

    
   // Update Theta for Time t (Equation 7)
      double inverse = arma::as_scalar(1 / (h + z_t * r_t * z_t.t()));
      theta_new = theta + r_t * z_t.t() * inverse * (y_t - z_t * theta);
   // theta_new = theta + r_t * z_t.t() * inv(h + z_t * r_t * z_t.t()) * (y_t - z_t * theta);
     
   // Update Var-Cov-Matrix for Time t (Equation 8)
      cov_mat_new = r_t - r_t * z_t.t() * inverse * (z_t * r_t);
   // cov_mat_new = r_t - r_t * z_t.t() * inv(h + z_t * r_t * z_t.t()) * (z_t * r_t);
    
   // Update H for Time t (Equation 10)
      double z_times_theta = arma::as_scalar(z_t * theta_new);
                     h_new = arma::as_scalar(kappa * h + (1 - kappa) * pow(y_t - z_times_theta, 2));
       
   // Get Predictive Density for Predicting t + 1 (Equation 9)
            mu = arma::as_scalar(z_pred * theta_new);
      variance = arma::as_scalar(h_new + z_pred * ((1 / lambda) * cov_mat_new) * z_pred.t());
      
   // Fill Return-List  
      List ret(5); 
      ret[0] = theta_new;
      ret[1] = cov_mat_new;
      ret[2] = h_new;
      ret[3] = mu;
      ret[4] = variance;
   
   // Return List  
      return ret;
}


// 1) TV-C Regressions
// Function III - Loop over Signals, Lambda and Kappa
// [[Rcpp::export]]
   List tvc_model_cand_(const double& y_t, 
                        const arma::rowvec& s_t,
                        const arma::rowvec& s_tplus1, 
                        const arma::vec& lambda_grid, 
                        const arma::vec& kappa_grid, 
                        List& theta_all, 
                        List& cov_mat_all, 
                        List& h_all) {      
  
   // Define Variables
      List tvc_results(1);
      int n_cands = s_t.size() * lambda_grid.size() * kappa_grid.size();
      NumericVector mu_vec(n_cands);
      NumericVector variance_vec(n_cands);
    
   // Loop over Lambda
      int counter = 0;
      for (int l = 0; l < lambda_grid.size(); l++) {

      // Set Lambda
         double lambda = lambda_grid(l);

      // Loop over Kappa
         for (int k = 0; k < kappa_grid.size(); k++) {

         // Set Kappa
            double kappa =  kappa_grid(k);

         // Loop over all candidates
            for (int j = 0; j < s_t.size(); j++) {
        
            // Set Signals
               double      s_t_j = s_t(j);
               double s_tplus1_j = s_tplus1(j);

            // Set Theta, Coefficient Variance, Observational Variance
               arma::mat theta;
               arma::mat cov_mat;
                  theta = as<arma::mat>(theta_all(counter));
                cov_mat = as<arma::mat>(cov_mat_all(counter));
                  double h = h_all(counter);

            // Check if signal is NA or not 
               if( !NumericVector::is_na(s_t_j) ) {
   
               // Apply TV-C-Function
                  tvc_results = tvc_model_(y_t,
                                           s_t_j,
                                           s_tplus1_j, 
                                           lambda,
                                           kappa,
                                           theta ,
                                           cov_mat,
                                           h);

               // Assign TV-C-Model-Results 
                  theta_all(counter)    = as<arma::mat>(tvc_results(0));
                  cov_mat_all(counter)  = as<arma::mat>(tvc_results(1));
                  h_all(counter)        = tvc_results(2);
                  mu_vec(counter)       = tvc_results(3); 
                  variance_vec(counter) = tvc_results(4);
                  counter++;

               } else if ( NumericVector::is_na(s_t_j) ) {

               // Assign TV-C-Model-Results 
                  theta_all(counter)    = theta;
                  cov_mat_all(counter)  = cov_mat;
                  h_all(counter)        = h;
                  mu_vec(counter)       = NA_REAL; 
                  variance_vec(counter) = NA_REAL;
                  counter++;

               } else {
                  throw std::invalid_argument("Something is wrong with the NA-Values");
               }
            }
         }
      }
      
      // Fill list 
         List ret(2);
         ret[0] = mu_vec;
         ret[1] = variance_vec;
     
      // Return list
         return ret;
}


// ################################################################


// 2) Dynamic Subset Combination
// Function I - Initialize DSC-Parameter
//[[Rcpp::export]]
   List dsc_init_(const int& n_cands,
                  const int& n_combs,
                  const int& n_gamma,
                  IntegerVector na_idx) { 

   // Initialize Vector for Discounted-Log-Likelihood-Score (Subset Combinations)  
      NumericVector dpll_combs(n_combs, 0.0);
   
   // Initialize List for Discounted-Log-Likelihood-Score (Candidate Models)  
      List dpll_cands(n_gamma);
      for (int g=0; g<n_gamma; g++) {
          NumericVector vec(n_cands, 0.0);
          vec[na_idx] = NA_REAL;
          dpll_cands(g) = vec; 
      }
  
   // Fill Return-List 
      List ret(2);
      ret[0] = dpll_cands;
      ret[1] = dpll_combs;

   // Return Vector
      return ret;
}


// Function II - Rank and Set Active Model Set (Active Models)
//[[Rcpp::export]]
   IntegerVector dsc_active_models_(const NumericVector& dpll_cands_gamma,  
                                    const int& psi){
 
   // Make sure Psi is equal or smaller to the number of non-na-values in dpll_cands_gamma
      int non_na_ctr = sum(!is_na(dpll_cands_gamma));  
      int psi_ = std::min(non_na_ctr, psi);
     
   // Set up Vector with Col-Indices
      IntegerVector idx = seq(0, dpll_cands_gamma.size()-1);
   
   // Check if all elements are equal (without NAs)
      NumericVector vec = dpll_cands_gamma[!is_na(dpll_cands_gamma)];
      bool isequal = is_true(all(vec == vec[0]));
 
   // If all elements are equal, return the indices from left to right
      int i = 0;
      if (isequal) {
         IntegerVector indices;
         while(indices.length() < psi_) {
            if (!std::isnan(dpll_cands_gamma[i])) {
               indices.push_back(i); 
            }
            i++;
         }

      // Fill idx with indices
         idx = indices;

    // If the elementd are not equal, use partial sort 
      } else {

      // Get psi-highest values (-> indices)
         std::nth_element(idx.begin(), idx.begin()+psi_, idx.end(), 
                          [&](int i, int j) {
                           if( std::isnan(dpll_cands_gamma[i]) ) return false;
                           if( std::isnan(dpll_cands_gamma[j]) ) return true;
                           return dpll_cands_gamma[i] > dpll_cands_gamma[j];
                           });
 
      // Sort the highest Values (-> indices)
         std::sort(idx.begin(), idx.begin()+psi_,
                   [&](int i, int j) {
                     if( std::isnan(dpll_cands_gamma[i]) ) return false;
                     if( std::isnan(dpll_cands_gamma[j]) ) return true;
                     return dpll_cands_gamma[i] > dpll_cands_gamma[j];
                     });
      }

    // Only return the first psi_ indices
       idx = idx[seq(0, psi_ - 1)];
   
   // Return list
      return idx;
}


// Function III - Compute Aggregated Predictive Distribution
//[[Rcpp::export]]
   List dsc_agg_density_(const NumericVector& active_weights, 
                         const NumericVector& forecast_tvc_t,             
                         const NumericVector& variance_tvc_t, 
                         const IntegerVector& idx_sub) {
                                     
   // Define Variables
      List ret(2);
      double mu_comb, variance_comb; 

   // Subset Matrices (select only active models)
      NumericVector oos_forecast_tvp_sub = forecast_tvc_t[idx_sub];
      NumericVector oos_variance_tvp_sub = variance_tvc_t[idx_sub];

   // Calculate Combined Predictive Density (Logarithmic Combination Rule)
      variance_comb = 1 / sum(active_weights / oos_variance_tvp_sub);
            mu_comb = sum(active_weights * oos_forecast_tvp_sub / oos_variance_tvp_sub) *
                          variance_comb;
   
   // Fill list  
      ret[0] = mu_comb;
      ret[1] = variance_comb;

   // Return list
       return ret;
}


// Function IV - Calculate (exponentially down-weighted) Log-Predictive-Likelihoods for Candidate Forecasts
//[[Rcpp::export]]
   NumericVector dsc_dpll_tvc_(NumericVector dpll_cands_gamma, 
                               const double& y_t, 
                               const NumericVector& forecast_tvc_t,             
                               const NumericVector& variance_tvc_t,
                               const double& gamma,
                               const int& method = 1,
                               Nullable<const NumericVector&> risk_aversion_ = R_NilValue,
                               Nullable<const NumericVector&> min_weight_ = R_NilValue,
                               Nullable<const NumericVector&> max_weight_ = R_NilValue) {
                                     
   // Define Variables
      int n_cands = dpll_cands_gamma.length();
      NumericVector performance_score(n_cands, NA_REAL);

   // Calculate Performance
      for (int i=0; i<n_cands; i++) {

      // Check for NA value
         if (!NumericVector::is_na(forecast_tvc_t(i))) {  

         // Optimization-Method
         // 1) Log-Likelihood
            if (method == 1) {
               performance_score(i) = R::dnorm(y_t,   
                                               forecast_tvc_t(i),
                                               pow(variance_tvc_t(i), 0.5),
                                               true);

         // 2) MSE
            } else if (method == 2) {
               performance_score(i) = -pow(y_t - forecast_tvc_t(i), 2.0);

         // 3) AE
            } else if (method == 3) {
               performance_score(i) = -abs(y_t - forecast_tvc_t(i));
         
         // 4) Compounded Returns
            } else if (method == 4) {

               // Check if relevant parameter are provided
                  if(risk_aversion_.isNull() && min_weight_.isNull() && max_weight_.isNull()) {
                     throw std::invalid_argument("Error: Relevant parameter not provided!");
                  }

               // Cast to type double
                  double risk_aversion = as<double>(risk_aversion_);
                  double min_weight = as<double>(min_weight_);
                  double max_weight = as<double>(max_weight_);

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
         
         // 5) Catch Error   
            } else {
               throw std::invalid_argument("Error: Method not available");
            }
         }
      }

   // Exponentially down-weighted past predictive log-likelihoods
      dpll_cands_gamma = dpll_cands_gamma * gamma + performance_score * gamma;
     
   // Return Updated Weights  
      return dpll_cands_gamma;
}


// Function V - Rank and Select Aggregated Forecast
//[[Rcpp::export]]
   List rank_comb_(const NumericVector& dpll_combs, 
                   const NumericVector& mu_comb_vec,
                   const NumericVector& variance_comb_vec) {
                              
   // Define Variables
      List ret(3);
      int best_idx;
      double forecast_stsc, variance_stsc;
    
   // Get Index of highest Rank
      best_idx = which_max(dpll_combs);
   
   // Select STSC-Forecast
      forecast_stsc = mu_comb_vec(best_idx);
      variance_stsc = variance_comb_vec(best_idx);

   // Fill list  
      ret[0] = forecast_stsc;
      ret[1] = variance_stsc;
      ret[2] = best_idx;

   // Return list
      return ret;
}


// Function VI - Calculate (exponentially down-weighted) Log-Predictive-Likelihoods for Combinations (Aggregated Predictive Distributions)
//[[Rcpp::export]]
   NumericVector dsc_dpll_comb_(NumericVector& dpll_combs, 
                                const double& y_t, 
                                const NumericVector& forecasts_comb,             
                                const NumericVector& variances_comb,
                                const double& delta,
                                const int& method = 1,
                                Nullable<const NumericVector&> risk_aversion_ = R_NilValue,
                                Nullable<const NumericVector&> min_weight_ = R_NilValue,
                                Nullable<const NumericVector&> max_weight_ = R_NilValue) {
                                     
   // Define Variables
      int n_combs = dpll_combs.length();
      NumericVector dpll_combs_new(n_combs);
      NumericVector performance_score(n_combs);

   // Calculate Performance
      for( int i=0; i<n_combs; i++) {

      // Optimization-Method
      // 1) Log-Likelihood
         if (method == 1) {
         performance_score(i) = R::dnorm(y_t,
                                         forecasts_comb(i),
                                         pow( variances_comb(i), 0.5 ),
                                         true); 

      // 2) SE
         } else if (method == 2) {
            performance_score(i) = -pow(y_t - forecasts_comb(i), 2.0);

      // 3) AE
         } else if (method == 3) {
            performance_score(i) = -abs(y_t - forecasts_comb(i));

      // 4) Compounded Returns
         } else if (method == 4) {

         // Check if relevant parameter are provided
            if(risk_aversion_.isNull() && min_weight_.isNull() && max_weight_.isNull()) {
               throw std::invalid_argument("Error: Relevant parameter not provided!");
            }

         // Cast to type double
            double risk_aversion = as<double>(risk_aversion_);
            double min_weight = as<double>(min_weight_);
            double max_weight = as<double>(max_weight_);

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

      // 5) Catch Error   
         } else {
            throw std::invalid_argument("Error: Method not available");
         }
      }
    
   // Exponentially down-weighted Past Predictive Log-Likelihoods
      dpll_combs_new = delta * dpll_combs + performance_score;

    // Return list
       return dpll_combs_new;
}


// Function VII - Loop over Gamma and Psi
//[[Rcpp::export]]
   List dsc_loop_(List dpll_cands,
                  NumericVector& dpll_combs, 
                  const NumericVector& gamma_grid, 
                  const IntegerVector& psi_grid, 
                  const double& y_t,
                  const NumericVector& forecast_tvc_t,             
                  const NumericVector& variance_tvc_t,
                  const double& delta,
                  const int& method = 1,
                  Nullable<const NumericVector&> risk_aversion_ = R_NilValue,
                  Nullable<const NumericVector&> min_weight_ = R_NilValue,
                  Nullable<const NumericVector&> max_weight_ = R_NilValue) { 

   // Define Variables
      int n_combs = dpll_combs.length();
      List chosen_cands(n_combs);
      NumericVector forecasts_comb(n_combs), variances_comb(n_combs);
      List stsc_results(1);
    
   // Loop over Gamma and Psi
      int ctr = 0;
      for (int g=0; g<gamma_grid.length(); g++) {

         // Set Gamma
            double gamma = gamma_grid(g);

         // Compute Active Set of Candidate Models
                  int psi_max = max(psi_grid);
            IntegerVector idx = dsc_active_models_(as<NumericVector>(dpll_cands(g)), psi_max);                                    

         // Loop over Psi
            for (int p=0; p<psi_grid.length(); p++) {

            // Set Psi
               int non_na_ctr =sum(!is_na(as<NumericVector>(dpll_cands(g))));
               int psi = std::min(non_na_ctr, psi_grid(p));

            // Define Variables
               List agg_density(1);
   
            // Get the Active Set of Candidate Models
               NumericVector active_weights(psi, 1.0 / psi);
               chosen_cands(ctr) = idx[seq(0, psi - 1)];
   
            // Calculate Aggregated Predictive Density
               agg_density = dsc_agg_density_(active_weights,
                                              forecast_tvc_t,
                                              variance_tvc_t,
                                              chosen_cands(ctr)); 
   
            // Assign Results
               forecasts_comb(ctr) = agg_density(0);
               variances_comb(ctr) = agg_density(1);
               ctr++;
            }

         // Update DPLL for Candidate Models
            dpll_cands(g) = dsc_dpll_tvc_(dpll_cands(g), 
                                          y_t, 
                                          forecast_tvc_t,             
                                          variance_tvc_t,
                                          gamma,
                                          method,
                                          risk_aversion_,
                                          min_weight_,
                                          max_weight_);
      }

      // Select Aggregated Forecast
         stsc_results = rank_comb_(dpll_combs, 
                                   forecasts_comb,
                                   variances_comb);

      // Assign Results
         double stsc_forecast = stsc_results(0);
         double stsc_variance = stsc_results(1);
         int stsc_idx = stsc_results(2);

      // Update DPLL for Combinations (Aggregated Predictive Distributions)
         dpll_combs = dsc_dpll_comb_(dpll_combs, 
                                     y_t, 
                                     forecasts_comb,             
                                     variances_comb,
                                     delta,
                                     method,
                                     risk_aversion_,
                                     min_weight_,
                                     max_weight_);
   
      // Fill list   
         List ret(4); 
         ret[0] = stsc_forecast;
         ret[1] = stsc_variance;
         ret[2] = stsc_idx;
         ret[3] = chosen_cands(stsc_idx);

      // Return list
         return ret;
}


// 3.) STSC 
// Function I - Loop over t
// [[Rcpp::export]]
   List stsc_loop(const arma::vec& y, 
                  Nullable<const NumericMatrix&> X_, 
                  Nullable<const NumericMatrix&> F_, 
                  const int& sample_length,
                  const arma::vec& lambda_grid,
                  const arma::vec& kappa_grid,
                  const int& burn_in_tvc,
                  const NumericVector& gamma_grid,
                  const IntegerVector& psi_grid,
                  const double& delta,
                  const int& burn_in_dsc,
                  const int& method = 1,
                  Nullable<const NumericVector&> risk_aversion_ = R_NilValue,
                  Nullable<const NumericVector&> min_weight_ = R_NilValue,
                  Nullable<const NumericVector&> max_weight_ = R_NilValue) { 

   
   // Check whether Simple Signals and / or Point Forecasts are provided ...
   // ... and create combined Signal-Matrix
      arma::mat S;
      int n_sim_sig, n_point_f, n_signal, n_cands;
     
      if (X_.isNull() && F_.isNull()) {
         throw std::invalid_argument("Error: No signals provided");
      } else if (X_.isNotNull() && F_.isNotNull()) {
         n_sim_sig = as<arma::mat>(X_).n_cols;
         n_point_f = as<arma::mat>(F_).n_cols;
         n_signal  = n_sim_sig + n_point_f;
         n_cands   = n_signal * lambda_grid.size() * kappa_grid.size();
         S = arma::join_rows(as<arma::mat>(X_), as<arma::mat>(F_));
      } else if (F_.isNull()) {
         n_sim_sig = as<arma::mat>(X_).n_cols;
         n_point_f = 0;
         n_signal  = n_sim_sig + n_point_f;
         n_cands   = n_signal * lambda_grid.size() * kappa_grid.size();
         S = as<arma::mat>(X_);
      } else {
         n_sim_sig = 0;
         n_point_f = as<arma::mat>(F_).n_cols;
         n_signal  = n_sim_sig + n_point_f;
         n_cands   = n_signal * lambda_grid.size() * kappa_grid.size();
         S = as<arma::mat>(F_);
      }
   
   // Define Variables for TV-C-Models
      List theta_all(n_cands), cov_mat_all(n_cands), h_all(n_cands);
      IntegerVector na_idx;
      NumericVector forecast_tvc_tplus1(n_cands, 0.0), variance_tvc_tplus1(n_cands, 0.0);

   // Define Variables for Dynamic Subset Combinations
      int tlength = y.size();
      int n_combs = gamma_grid.length() * psi_grid.length();
      NumericVector stsc_forecast(tlength, NumericVector::get_na());
      NumericVector stsc_variance(tlength, NumericVector::get_na());
      IntegerVector stsc_idx(tlength, IntegerVector::get_na());
      List chosen_cands(tlength);

   // --- 
   // Apply TV-C-Init-Function
      List init_tvc_results(1);
      init_tvc_results = init_tvc_(y, 
                                   S,
                                   n_sim_sig,
                                   sample_length,
                                   lambda_grid,
                                   kappa_grid);

   // Assign Results
      theta_all   = init_tvc_results(0);
      cov_mat_all = init_tvc_results(1);
      h_all       = init_tvc_results(2);
      na_idx      = init_tvc_results(3);

   // ---
   // Apply DSC-Init-Function
      List init_dsc_results(1);
      init_dsc_results = dsc_init_(n_cands,
                                   n_combs,
                                   gamma_grid.length(),
                                   na_idx); 

   // Assign Results
      List          dpll_cands = init_dsc_results(0);          
      NumericVector dpll_combs = clone(as<NumericVector>(init_dsc_results(1))); 
    
   // --- 
   // Loop over t
      for ( int t=0; t<(tlength-1); t++ ) {

      // Subset Data
         double             y_t = y[t];
         double        y_tplus1 = y[t+1];
         arma::rowvec       s_t = S.row(t);
         arma::rowvec  s_tplus1 = S.row(t+1);

      // Check for NA-Values in Candidate Models in t  !! Muss nur gemacht werden, wenn NA-Werte in Candidate Models !!
         IntegerVector na_idx_t;
         int ctr = 0;
            for (int l = 0; l < lambda_grid.size(); l++) {
               for (int k = 0; k < kappa_grid.size(); k++) {
                  for (int j = 0; j < S.n_cols; j++) {

                     // Check and Count for NA-Values
                        if (NumericVector::is_na(s_t(j))) {
                           na_idx_t.push_back(ctr);
                        }
                        ctr++;
                  }
               } 
            }

      // Identify Candidate Models that went to Non-Na
         if (na_idx_t.length() < na_idx.length()) {

            // Get the Index for the Signals that are not NA anymore
               IntegerVector vec_diff = setdiff(na_idx, na_idx_t);
                               na_idx = clone(na_idx_t);              
               for ( int g=0; g<gamma_grid.length(); g++ ) {
                   as<NumericVector>(dpll_cands(g))[vec_diff] = 0.0; }
         }

      // Check for Burn-In-Period
         if (t == (burn_in_tvc-1)) {
            dpll_cands = dsc_init_(n_cands, n_combs, gamma_grid.length(), na_idx)(0);     
            dpll_combs = clone(as<NumericVector>(init_dsc_results(1)));
            stsc_forecast.fill(NA_REAL); 
            stsc_variance.fill(NA_REAL); 
            stsc_idx.fill(NA_INTEGER);   
            IntegerVector idx = Range(0, t);
            chosen_cands[idx] = NA_INTEGER;
         }
         
         if (t == (burn_in_dsc-1)) {
            dpll_combs = clone(as<NumericVector>(init_dsc_results(1)));
            stsc_forecast.fill(NA_REAL); 
            stsc_variance.fill(NA_REAL); 
            stsc_idx.fill(NA_INTEGER);   
            IntegerVector idx = Range(0, t);
            chosen_cands[idx] = NA_INTEGER;
         }

      // Apply TV-C-Function
         List tvc_model_cand_results(1);
         tvc_model_cand_results = tvc_model_cand_(y_t, 
                                                  s_t,
                                                  s_tplus1, 
                                                  lambda_grid, 
                                                  kappa_grid, 
                                                  theta_all, 
                                                  cov_mat_all, 
                                                  h_all); 

      // Assign Results
         forecast_tvc_tplus1 = tvc_model_cand_results(0);
         variance_tvc_tplus1 = tvc_model_cand_results(1);

      // Apply DSC-Function
         List dsc_results(1);
         dsc_results = dsc_loop_(dpll_cands,
                                 dpll_combs, 
                                 gamma_grid, 
                                 psi_grid, 
                                 y_tplus1,
                                 forecast_tvc_tplus1,             
                                 variance_tvc_tplus1,
                                 delta,
                                 method,
                                 risk_aversion_,
                                 min_weight_,
                                 max_weight_); 

      // Assign Results
         stsc_forecast(t+1) = dsc_results(0);
         stsc_variance(t+1) = dsc_results(1);
         stsc_idx(t+1)      = dsc_results(2);
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
