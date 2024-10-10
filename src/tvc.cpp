#include <RcppArmadillo.h>
using namespace Rcpp;

// ###################################################################################
// 1) TV-C Regressions
// Function I - Initialize TV-C-Parameter
   List init_tvc(const arma::vec& y, 
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
      arma::cube   theta_cube(2, 1, n_cands);
      arma::cube   cov_mat_cube(2, 2, n_cands);
      arma::rowvec h_vec(n_cands);
      arma::vec    y_sample, x;
      arma::uvec   non_finite, init_idx;
      arma::mat    x_sample_one, theta, cov_mat;
      arma::colvec coef;
      double       intercept, var_y, var_x, h;
      List         ret_all(3);  

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
   
               // Index for subsetting the Initialisation Sample
                  init_idx = arma::regspace<arma::uvec>(0 + na_ctr, na_ctr + init - 1);
   
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
                     if (var_x != 0.0) {
                          cov_mat(1, 1) =  var_y / var_x; 
                     } else {
                          cov_mat(1, 1) = var_y;
                     }
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
                  theta_cube.slice(ctr)   = theta;
                  cov_mat_cube.slice(ctr) = cov_mat; 
                  h_vec(ctr)              = h;
                  ctr++;
          }
        }
      }
   
      // Fill Return-List  
         ret_all[0] = theta_cube;
         ret_all[1] = cov_mat_cube; 
         ret_all[2] = h_vec;
     
      // Return Results
         return ret_all;
   }
// ----------

// Function II - Predictive and Update Step for Signal j and time t
   arma::field<double> tvc_model(double y_t, 
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
// ----------

// Function III - Loop over Signals, Lambda and Kappa
   arma::field<arma::rowvec> tvc_model_cand(double y_t, 
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
   
               // Apply TV-C-Function
                  const arma::field<double> tvc_results = tvc_model(y_t,
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

               // Update Counter
                  counter++;
            }
         }
      }
      
      // Fill list 
         ret(0) = mu_vec;
         ret(1) = variance_vec;
     
      // Return list
         return ret;
}
// ----------

// Function IV - Loop over t
// [[Rcpp::export]]
   List tvc_(const arma::vec& y, 
             Nullable<const NumericMatrix&> X_, 
             Nullable<const NumericMatrix&> Ext_F_, 
             int init,
             const arma::vec& lambda_grid,
             const arma::vec& kappa_grid,
             bool bias) { 
   
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

   // Number of Candiate Models and Signals
      int n_signal = n_raw_sig + n_point_f;
      int n_cands  = n_signal * lambda_grid.n_elem * kappa_grid.n_elem;
   
   // Define Variables for TV-C-Models
      int tlength = y.n_elem;
      arma::mat tvc_forecast(tlength, n_cands); tvc_forecast.fill(arma::datum::nan);
      arma::mat tvc_variance(tlength, n_cands); tvc_variance.fill(arma::datum::nan);
      arma::cube theta_cube(2, 1, n_cands), cov_mat_cube(2, 2, n_cands);
      arma::rowvec h_vec(n_cands);
      arma::uvec current_na_cm, new_na_cm;

   // --- 
   // Apply TV-C-Init-Function
      List init_tvc_results(1);
      init_tvc_results = init_tvc(y, 
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
    
   // --- 
   // Loop over t
      for (int t=0; t<(tlength-1); t++ ) {

      // Subset Data
         const double           y_t = y[t];
         const arma::rowvec     s_t = S.row(t);
         const arma::rowvec  s_pred = S.row(t+1);
         
      // Apply TV-C-Function
         arma::field<arma::rowvec> tvc_model_cand_results(2);
         tvc_model_cand_results = tvc_model_cand(y_t, 
                                                 s_t,
                                                 s_pred, 
                                                 lambda_grid, 
                                                 kappa_grid, 
                                                 theta_cube, 
                                                 cov_mat_cube, 
                                                 h_vec); 

      // Assign Results
         tvc_forecast.row(t+1) = tvc_model_cand_results(0);
         tvc_variance.row(t+1) = tvc_model_cand_results(1);
      }

   // Fill list 
      List ret(2);   
      ret[0] = tvc_forecast;
      ret[1] = tvc_variance;

   // Return list
      return ret;
}
