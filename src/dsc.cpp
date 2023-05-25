#include <Rcpp.h>
using namespace Rcpp;

// 2) Dynamic Subset Combination
// Function I - Initialize all Predictive Densities
//[[Rcpp::export]]
   NumericVector init_dsc(int number_forecasts){

  // Define Variables  
     NumericVector weights(number_forecasts, 1.0 / number_forecasts);

  // Return Vector
     return weights;
}

// Function II - Re-Weight Probabilities by Forgetting Factor Gamma
//[[Rcpp::export]]
   NumericVector forget_dsc(NumericVector weights, 
                            double gamma){

    // Define Variables
       weights = pow(weights, gamma) / sum(pow(weights, gamma));
  
    // Return Updated Probabilities
       return weights;
}

// Function III - Rank and Subset Forecasting Models (Active Models)
//[[Rcpp::export]]
   List active_models_dsc(NumericVector weights, 
                          int psi){
                              
    // Define Variables
       List ret(2);
       NumericVector active_weights(psi);
    
    // (Partial) Sort Weights
       IntegerVector idx = seq(0, weights.size()-1);
       std::partial_sort(idx.begin(), idx.begin()+psi, idx.end(), 
                         [&](int i, int j) {return weights[i] > weights[j]; });
                      
    // Get index of the 'psi' highest weights
       IntegerVector idx_sub = idx[seq(0, psi-1)];
   
    // Calculate Active Weights
       active_weights.fill(1.0 / psi); 

        //active_weights = weights[idx_sub];
        //active_weights = active_weights / sum(active_weights); 
    
    // Fill list  
       ret[0] = active_weights;
       ret[1] = idx_sub;

    // Return list
       return ret;
}

//-----------------
// Helper Function - Subset Matrix based on Column- and Row-Index
// [[Rcpp::export]]
   NumericVector matrix_subset_idx(NumericMatrix mat, 
                                   IntegerVector col_idx,
                                   int t) { 

    // Determine Number of Columns 
       int n_cols_out = col_idx.size();

    // Create Output Vector
       NumericVector out(n_cols_out);

    // Loop through each column and copy data
       for(int i = 0; i < n_cols_out; ++i) {
               out(i) = mat(t, col_idx[i]);
       }
       return out;
}
// -----------------

// Function IV - Compute Combined Predictive Distribution
//[[Rcpp::export]]
   List agg_density_dsc(NumericVector active_weights, 
                        NumericVector oos_target_var, 
                        NumericMatrix oos_forecast_tvp,             
                        NumericMatrix oos_variance_tvp, 
                        IntegerVector idx_sub,
                        int t) {
                                     
    // Define Variables
       List ret(3);
       double mu_agg, variance_agg, ln_score;

    // Subset Matrices
       NumericVector oos_forecast_tvp_sub = matrix_subset_idx(oos_forecast_tvp, idx_sub, t);
       NumericVector oos_variance_tvp_sub = matrix_subset_idx(oos_variance_tvp, idx_sub, t);

    // Calculate Combined Predictive Density
       variance_agg = 1 / sum(active_weights / oos_variance_tvp_sub);
             mu_agg = sum(active_weights * oos_forecast_tvp_sub / oos_variance_tvp_sub) *
                          variance_agg;

    // Calculate Predictive Log Score
       ln_score = R::dnorm(oos_target_var(t), mu_agg, pow( variance_agg, 0.5 ), true ); 
    
    // Fill list  
       ret[0] = mu_agg;
       ret[1] = variance_agg;
       ret[2] = ln_score;

    // Return list
       return ret;
}

// Function V - Update Probabilities using Bayes' Rule
//[[Rcpp::export]]
   NumericVector update_dsc(NumericVector weights, 
                            NumericVector oos_target_var, 
                            NumericMatrix oos_forecast_tvp,             
                            NumericMatrix oos_variance_tvp,
                            int n_models, 
                            int t) {
                                     
    // Define Variables
       NumericVector pred_lik(n_models);
       double oos_target_var_t = oos_target_var(t);

    // Calculate Predictive Likelihood
       for( int i=0; i<n_models; i++) {
           pred_lik(i) = R::dnorm( oos_target_var_t,   
                                   oos_forecast_tvp(t, i), 
                                   pow(oos_variance_tvp(t, i), 0.5),
                                   false);
       }

    // Update Using Bayes' rule
       weights = (weights * pred_lik) / sum(weights * pred_lik);
    
    // Return Updated Weights  
       return weights;
}

// Function VI - Loop over Predictive and Update Step
// [[Rcpp::export]]
    List dsc_loop(NumericVector weights, 
                  double gamma, 
                  int psi, 
                  NumericVector oos_target_var,
                  NumericMatrix oos_forecast_tvp,             
                  NumericMatrix oos_variance_tvp,
                  int len_para_grid,
                  int oos_length, 
                  int n_models) { 

  // Define Variables
     List active_results(1), agg_density(1), ret(4);
     NumericVector active_weights(psi), forecasts_agg(oos_length), variances_agg(oos_length), ln_scores(oos_length);
     List selected_models(oos_length);
    
  // Start loop
      for (int t = 0; t < oos_length; t++) {

        // Forget Weights
           weights = forget_dsc(weights, gamma);

        // Active Models
           active_results     = active_models_dsc(weights, psi);                                    
           active_weights     = active_results(0);
           selected_models(t) = active_results(1);


        // Aggregated Predictive Density and Log Score
           agg_density = agg_density_dsc(active_weights,
                                         oos_target_var,
                                         oos_forecast_tvp,
                                         oos_variance_tvp,
                                         active_results(1),
                                         t); 

        // Assign Results
           forecasts_agg(t) = agg_density(0);
           variances_agg(t) = agg_density(1);
           ln_scores(t)     = agg_density(2);

        // Update Weights 
           weights = update_dsc(weights,
                                oos_target_var,
                                oos_forecast_tvp,
                                oos_variance_tvp,
                                n_models,
                                t);
     }
      
  // Fill list    
     ret[0] = forecasts_agg;
     ret[1] = variances_agg;
     ret[2] = ln_scores;
     ret[3] = selected_models;

 // Return list
    return ret;
}
