#include <RcppArmadillo.h>
using namespace Rcpp;

// 1) TV-C Regressions
// Function I a - Initialize TV-C Models ('Raw' Predictor Time Series)
// [[Rcpp::export]]
   List init_tvc(arma::vec& y_var, 
                 arma::vec& x_var, 
                 int sample_length) {
  
  // Define Variables
    List ret(3);
    arma::mat x_sample_one, theta, cov_mat;
    double intercept, var_y, var_x, h;

  // Index for subsetting
    arma::uvec index = arma::regspace<arma::uvec>(0, sample_length - 1);
    
  // Define and prepare matrices for regression  
    arma::vec y_sample = y_var.elem(index);  
    arma::mat x_sample(x_var.elem(index));   
    x_sample_one = arma::ones<arma::mat>(x_sample.n_rows, 1); 
    x_sample.insert_cols(0, x_sample_one);
  
  // Initialize - Theta
    theta = arma::zeros<arma::mat>(2,1); 
  
  // Initialize - System Covariance
     arma::colvec coef = solve(x_sample, y_sample);
     intercept = coef(0);
     var_y     = arma::var(y_sample);         
     var_x     = arma::var(x_sample.col(1));  
     cov_mat   = arma::zeros<arma::mat>(2, 2);
     cov_mat(0, 0) = pow(intercept, 2) + var_y;
     if(arma::is_finite(var_y / var_x)) {
          cov_mat(1, 1) =  var_y / var_x; 
     } else {
          cov_mat(1, 1) = 0;
     }
  
  // Initialize - Observational Variance
     h = var_y;
  
  // Fill Return-List  
    ret[0] = theta;
    ret[1] = cov_mat; 
    ret[2] = h;
  
  // Return Results
    return ret;
}

// Function I b - Initialize TV-C Models (Point Forecast Time Series)
// [[Rcpp::export]]
   List init_tvc_forecast(arma::vec& y_var, 
                          arma::vec& x_var, 
                          int sample_length) {
  
  // Define Variables
    List ret(3);
    arma::mat x_sample_one, theta, cov_mat;
    double intercept, var_y, var_x, h;

  // Index for subsetting
    arma::uvec index = arma::regspace<arma::uvec>(0, sample_length - 1);
    
  // Define and prepare matrices for regression  
    arma::vec y_sample = y_var.elem(index);  
    arma::mat x_sample(x_var.elem(index));   
    x_sample_one = arma::ones<arma::mat>(x_sample.n_rows, 1); 
    x_sample.insert_cols(0, x_sample_one);
  
  // Initialize - Theta
    theta = arma::zeros<arma::mat>(2,1); 
    theta(1, 0) = 1;
  
  // Initialize - System Covariance
     arma::colvec coef = solve(x_sample, y_sample);
     intercept = coef(0);
     var_y     = arma::var(y_sample);         
     var_x     = arma::var(x_sample.col(1));  
     cov_mat   = arma::zeros<arma::mat>(2, 2);
     cov_mat(0, 0) =  pow(intercept, 2) + var_y;  // Set to Zero for Constant Intercept
     cov_mat(1, 1) =  0;
  
  // Initialize - Observational Variance
     h = var_y;
  
  // Fill Return-List  
    ret[0] = theta;
    ret[1] = cov_mat; 
    ret[2] = h;
  
  // Return Results
    return ret;
}

// Function II - Predictive and Update Step for time t
// [[Rcpp::export]]
   List tvc_model(arma::vec& y_var, 
                  arma::vec& x_var, 
                  int i, 
                  double lambda, 
                  double kappa, 
                  arma::mat theta, 
                  arma::mat cov_mat, 
                  double h) { 
  
  // Define Variables
    List ret(5);
    arma::mat z_t, z_pred, r_t, theta_new, cov_mat_new;
    double mu, variance, y_t, h_new;
  
  // Get Equity Premium for Time t
    y_t = y_var(i - 1);
 
  // Get Predictor for Time t
    z_t  = arma::zeros<arma::mat>(1, 2);  
    z_t(0, 0) = 1;
    z_t(0, 1) = x_var(i - 1); 
  
  // Get Predictor for Predicting t + 1
    z_pred = arma::zeros<arma::mat>(1, 2);  
    z_pred(0, 0) = 1;
    z_pred(0, 1) = x_var(i); 
  
  // Calculate R for Time t (Equation 5)
     r_t = (1 / lambda) * cov_mat;
   
  // Update Theta for Time t (Equation 7)
     theta_new = theta + r_t * z_t.t() * inv(h + z_t * r_t * z_t.t()) * (y_t - z_t * theta);
    
  // Update Var-Cov-Matrix for Time t (Equation 8)
     cov_mat_new = r_t - r_t * z_t.t() * inv(h + z_t * r_t * z_t.t()) * (z_t * r_t);
   
  // Update H for Time t (Equation 10)
     double z_times_theta = arma::as_scalar(z_t * theta_new);
     h_new = arma::as_scalar(kappa * h + (1 - kappa) * pow(y_t - z_times_theta, 2));
      
  // Get Predictive Density for Predicting t + 1 (Equation 9)
     mu = arma::as_scalar(z_pred * theta_new);
     variance = arma::as_scalar(h_new + z_pred * ((1 / lambda) * cov_mat_new) * z_pred.t());
     
  // Fill List   
     ret[0] = theta_new;
     ret[1] = cov_mat_new;
     ret[2] = h_new;
     ret[3] = mu;
     ret[4] = variance;
  
 // Return List  
    return ret;
}

// Function III - Loop over t 
// [[Rcpp::export]]
   List tvc_model_loop(arma::vec& y_var, 
                       arma::vec& x_var, 
                       double lambda, 
                       double kappa, 
                       arma::mat theta, 
                       arma::mat cov_mat, 
                       double h,
                       int ts_length,
                       int drop_length,
                       int max_length) { 
  
  // Define Variables
     List tvc_results(1), ret(2);
     arma::vec mu_mat, variance_mat;
     mu_mat       = arma::mat(max_length, 1);  mu_mat.fill(arma::datum::nan);
     variance_mat = arma::mat(max_length, 1);  variance_mat.fill(arma::datum::nan);
    
  // Start loop
     for (int i = 0; i < ts_length; i++) {

       tvc_results = tvc_model(y_var,
                               x_var,
                               i + 1, 
                               lambda,
                               kappa,
                               theta ,
                               cov_mat,
                               h);
       
       // Assign TV-C-Model-Results 
          theta   = as<arma::mat>(tvc_results(0));
          cov_mat = as<arma::mat>(tvc_results(1));
           h      = tvc_results(2);
   
   // Save Predictive Density for time t + 1
       mu_mat(drop_length + i + 1)       = tvc_results(3);  
       variance_mat(drop_length + i + 1) = tvc_results(4);  
     }
      
  // Fill list    
     ret[0] = mu_mat;
     ret[1] = variance_mat;
     
 // Return list
    return ret;
}
