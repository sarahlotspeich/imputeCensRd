#' Multiple, semiparametric conditional mean imputation for a censored covariate
#'
#' Multiple, semiparametric conditional mean imputation for a censored covariate using a Cox proportional hazards model and the Breslow's estimator to estimate the conditional survival function.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param analysis_model Analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be a continuous outcome for normal linear regression.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param trapezoidal_rule A logical input for how to approximate the integral in the imputed values. Default is \code{trapezoidal_rule=FALSE} to use adaptive quadrature, but \code{"tr"}, and \code{trapezoidal_rule=TRUE} will use the trapezoidal rule.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"cf"} (carry forward, the default), \code{"wm"} (weighted mean), or \code{"m"} (mean).
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"d"} (immediate drop off), \code{"e"} (exponential extension, the default), or \code{"w"} (weibull extension).
#' @param B numeric, number of imputations. Default is \code{10}. 
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

cmi_sp_bootstrap = function(imputation_model, analysis_model, data, trapezoidal_rule = TRUE, Xmax = Inf, surv_between = "cf", surv_beyond = "e", maxiter = 100, B = 10) {
  # Size of resample -----------------------------------------------------------
  n = nrow(data)
  
  # Convert data to data frame -------------------------------------------------
  data = data.frame(data)
  
  # Extract variable names from imputation_model -------------------------------
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator
  Z = all.vars(imputation_model)[-c(1:2)] ## additional covariates
  
  # Take names and dimension from naive fit 
  data[, "imp"] = data[, W] ## start imputed value with observed 
  naive_fit = lm(formula = analysis_model, data = data) ## fit naive model 
  p = length(naive_fit$coefficients) ## dimension of beta vector 
  
  # Initialize empty dataframe to hold results 
  mult_fit = data.frame()
  
  # Loop through replicates 
  for (b in 1:B) {
    # Sample with replacement from the original data ---------------------------
    re_rows = ceiling(runif(n = n, min = 0, max = 1) * n)
    re_data = data[re_rows, ]
    
    # Calculate proportion of censored observations in resampled data ----------
    re_prop_cens = 1 - mean(re_data[, Delta])

    # Check for censoring in resampled data
    if (0 < re_prop_cens) {
      # Use imputeCensRd::cmi_sp() to impute censored x in re_data -------------
      re_data_imp = cmi_sp(imputation_model = imputation_model, 
                           data = re_data, 
                           trapezoidal_rule = trapezoidal_rule, 
                           Xmax = Xmax,
                           surv_between = surv_between, 
                           surv_beyond = surv_beyond)
    } else {
      # If no censored, just make imp = W and fit the usual model --------------
      re_data_imp = list(imputed_data = re_data,
                         code = TRUE)
      re_data_imp$imputed_data$imp = re_data[, W]
    }
    
    # Check for the Weibull extension not converging -------------------------
    while (!re_data_imp$code) {
      # Sample with replacement from the original data ---------------------------
      re_rows = ceiling(runif(n = n, min = 0, max = 1) * n)
      re_data = data[re_rows, ]
      
      # Calculate proportion of censored observations in resampled data ----------
      re_prop_cens = 1 - mean(re_data[, Delta])
      
      # Check for censoring in resampled data
      if (0 < re_prop_cens) {
        # Use imputeCensRd::cmi_sp() to impute censored x in re_data -------------
        re_data_imp = cmi_sp(imputation_model = imputation_model, 
                             data = re_data, 
                             trapezoidal_rule = trapezoidal_rule, 
                             Xmax = Xmax,
                             surv_between = surv_between, 
                             surv_beyond = surv_beyond)
      } else {
        # If no censored, just make imp = W and fit the usual model --------------
        re_data_imp = list(imputed_data = re_data,
                           code = TRUE)
        re_data_imp$imputed_data$imp = re_data[, W]
      }
    }
    
    # Once imputation was successful, fit the analysis model -----------------
    re_fit = lm(formula = analysis_model, 
                data = re_data_imp$imputed_data) 
    
    # Save coefficients to results matrix --------------------------------------
    mult_fit = rbind(mult_fit, 
                     summary(re_fit)$coefficients)
  } 
  
  ## Pool estimates
  pooled_est = rowsum(x = mult_fit[, 1],
                      group = rep(x = 1:p, times = B)) / B
  
  ## Calculate within-imputation variance 
  within_var = rowsum(x = mult_fit[, 2] ^ 2,
                      group = rep(x = 1:p, times = B)) / B
  
  ## Calculate between-imputation variance
  between_var = rowsum(x = (mult_fit[, 1] - rep(pooled_est, times = B)) ^ 2,
                       group = rep(x = 1:p, times = B)) / (B - 1)
  
  ## Pool variance
  pooled_var = within_var + between_var + (between_var / B)
  
  # Return table of pooled estimates
  tab = data.frame(Coefficient = names(naive_fit$coefficients),
                   Est = pooled_est,
                   SE = sqrt(pooled_var))
  return(tab)  
}