#include <RcppArmadillo.h>
using namespace Rcpp;

// 2) Dynamic Subset Combination
// Function I - Initialize DSC-Parameter
   arma::field<arma::field<arma::rowvec>> dsc_init(int n_cands,
                                                   int n_combs,
                                                   int n_gamma) {

   // Define Variables
      arma::field<arma::field<arma::rowvec>> ret(2); 

   // Initialize Vector for Performance-Score (Subset Combinations) -> Ranking
      arma::rowvec score_combs(n_combs, arma::fill::zeros);
   
   // Initialize Vector for Performance-Score (Candidate Models) -> Ranking 
      arma::rowvec vec(n_cands, arma::fill::zeros);

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


// Function II - Rank and Set Active Model Set (Active Models)
   arma::field<arma::uvec> dsc_active_models(const arma::rowvec& score_cands_gamma,  
                                             int psi) {

   // Define Variables
      arma::field<arma::uvec> ret(2);
      const int n_elem = score_cands_gamma.n_elem;
      const int psi_ = std::min(n_elem, psi);

   // Check if all elements are equal
      const bool isequal = arma::all(score_cands_gamma == score_cands_gamma(0));

   // If all elements are equal, return the indices from left to right
      arma::uvec idx;
      if (isequal) {
          idx = arma::regspace<arma::uvec>(0, psi_ - 1);
      } else {
          // Get psi-highest values (-> indices)
          const arma::uvec sorted_idx = arma::sort_index(score_cands_gamma, "descend");
          idx = sorted_idx.head(psi_);
      }

   // Fill Return-Field
      ret(0) = idx;
      ret(1) = arma::uvec(1).fill(n_elem);

   // Return
      return ret;
}

// Function III - Compute Aggregated Predictive Distribution
   arma::field<double> dsc_agg_density(const arma::rowvec& active_weights, 
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


// Function IV - Calculate (exponentially down-weighted) Performance Scores (-> Ranking) for Candidate Forecasts
   void dsc_score_cands(arma::rowvec& score_cands_gamma, 
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

   // Exponentially down-weighted past performance scores
      score_cands_gamma = score_cands_gamma * gamma + performance_score * gamma;
}


// Function V - Rank and Select Aggregated Forecast
   arma::field<double> rank_comb(const arma::rowvec& score_combs, 
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


// Function VI - Calculate (exponentially down-weighted) performance scores (-> ranking) for combinations (Aggregated Predictive Distributions)
   void dsc_score_comb(arma::rowvec& score_combs, 
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

// Helper Function to remove duplicates between arma::uvec and include_idx
   arma::uvec remove_duplicates(const arma::uvec& vec,
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

// Function VII - Loop over Gamma and Psi
   arma::field<arma::rowvec> dsc_loop(arma::field<arma::rowvec>& score_cands,
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
      const int psi_max = max(psi_grid);

   // Loop over Gamma and Psi
      int ctr = 0;
      for (unsigned int g=0; g<gamma_grid.n_elem; g++) {

         // Set Gamma
            double gamma = gamma_grid(g);

         // Compute Maximum Set of Active Candidate Models
            active_models = dsc_active_models(score_cands(g), psi_max);
            arma::uvec active_idx = active_models(0);
            const int non_na_ctr = active_models(1)(0);

         // Add Keep-Index & Remove Duplicates
            if (incl_idx.n_elem > 0) {
               active_idx = arma::join_cols(incl_idx, active_idx);
               active_idx = remove_duplicates(active_idx, incl_idx);
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
               agg_density = dsc_agg_density(active_weights,
                                             forecast_tvc_t,
                                             variance_tvc_t,
                                             active_idx_uvec); 
   
            // Assign Results
               forecasts_comb(ctr) = agg_density(0);
               variances_comb(ctr) = agg_density(1);
               ctr++;
            }

         // Update score for Candidate Models
            dsc_score_cands(score_cands(g), 
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
      stsc_results = rank_comb(score_combs, 
                               forecasts_comb,
                               variances_comb);

   // Assign Results
      double stsc_forecast = stsc_results(0);
      double stsc_variance = stsc_results(1);
      int stsc_idx = stsc_results(2);           

   // Update score for Combinations (Aggregated Predictive Distributions)
      dsc_score_comb(score_combs, 
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

// 3.) Wrapper Dynamic Subset Combination 
// Function I - Loop over t
// [[Rcpp::export]]
   List dsc_(const arma::vec& y, 
             const arma::mat& point_forecasts, 
             const arma::mat& variance_forecasts, 
             arma::rowvec gamma_grid, 
             arma::irowvec psi_grid, 
             double delta,
             int burn_in,
             int burn_in_dsc,
             int metric,
             bool equal_weight,
             Nullable<IntegerVector> incl_,
             Nullable<NumericVector> portfolio_params_) {

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
      
   // Define Variables for Dynamic Subset Combinations
      int tlength = y.n_elem;
      int n_cands  = point_forecasts.n_cols;
      int n_combs = gamma_grid.n_elem * psi_grid.n_elem;
      arma::vec  stsc_forecast(tlength); stsc_forecast.fill(arma::datum::nan);
      arma::vec  stsc_variance(tlength); stsc_variance.fill(arma::datum::nan);
      arma::vec stsc_idx(tlength); stsc_idx.fill(arma::datum::nan);  
      arma::field<arma::rowvec> score_cands(gamma_grid.n_elem);
      arma::rowvec score_combs(n_combs);
      List chosen_cands(tlength);

   // Include: CFM that must be included in the Subsets
      arma::uvec incl_idx;
      if (incl_.isNotNull()) {

      // Cast Values to uvec incl
         arma::uvec incl = as<arma::uvec>(incl_.get());

      // Number of CFM per Signal
         int grid_size = 1;

      // Resize incl_idx (-> Number of CFM that must be included)
         incl_idx.set_size(grid_size * incl.n_elem);

      // Calculate the Indices of the CFM
         int ctr = 0;
         for (arma::uword k = 0; k < incl.n_elem; ++k) {
           for (int i = 0; i < grid_size; ++i) {
             arma::uword index = (incl[k] - 1);
             incl_idx(ctr) = index;
             ctr++;
           }
         }
      }

   // ---
   // Apply DSC-Init-Function
      arma::field<arma::field<arma::rowvec>> init_dsc_results(2);
      init_dsc_results = dsc_init(n_cands,
                                  n_combs,
                                  gamma_grid.n_elem); 

   // Assign Results
      score_cands = init_dsc_results(0); 
      score_combs = init_dsc_results(1)(0); 

   // --- 
   // Loop over t
      for (int t=0; t<tlength; t++ ) {

      // Subset Data
         const double y_t = y[t];
         const arma::rowvec point_forecasts_t = point_forecasts.row(t);
         const arma::rowvec variance_forecasts_t = variance_forecasts.row(t);
         
      // Check for Burn-In-Period
         if (t == (burn_in-1)) {

            arma::field<arma::field<arma::rowvec>> init_dsc_results_after_burn_in = dsc_init(n_cands, n_combs, gamma_grid.n_elem);
            score_cands = init_dsc_results_after_burn_in(0); 
            score_combs = init_dsc_results_after_burn_in(1)(0); 
            stsc_forecast.fill(arma::datum::nan); 
            stsc_variance.fill(arma::datum::nan); 
            stsc_idx.fill(arma::datum::nan);   
            IntegerVector idx = Range(0, t);
            chosen_cands[idx] = NA_INTEGER;
         }
         
         if (t == (burn_in_dsc-1)) {

            arma::field<arma::field<arma::rowvec>> init_dsc_results_after_burn_in = dsc_init(n_cands, n_combs, gamma_grid.n_elem);
            score_combs = init_dsc_results_after_burn_in(1)(0); 
            stsc_forecast.fill(arma::datum::nan); 
            stsc_variance.fill(arma::datum::nan); 
            stsc_idx.fill(arma::datum::nan);   
            IntegerVector idx = Range(0, t);
            chosen_cands[idx] = NA_INTEGER;
         }

      // Apply DSC-Function
         arma::field<arma::rowvec> dsc_results(4);
         dsc_results = dsc_loop(score_cands,
                                score_combs, 
                                gamma_grid, 
                                psi_grid, 
                                y_t,
                                point_forecasts_t,
                                variance_forecasts_t,
                                delta,
                                metric,
                                equal_weight,
                                incl_idx,
                                risk_aversion,
                                min_weight,
                                max_weight); 

      // Assign Results
         stsc_forecast(t) = dsc_results(0)(0);
         stsc_variance(t) = dsc_results(1)(0);
         stsc_idx(t)      = dsc_results(2)(0);
         chosen_cands(t)  = dsc_results(3);
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
