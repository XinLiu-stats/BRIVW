#' BRIVW Estimation and BMEI-C Test
#'
#' @description
#' Performs bivariate rerandomized inverse variance weighted (BRIVW)
#' estimation for the causal effect and the accompanying combined
#' bivariate modified Egger intercept (BMEI-C) test for detecting
#' directional pleiotropy.
#'
#' The input summary statistics are assumed to have already undergone
#' allele harmonization, randomized instrument selection, and LD pruning.
#'
#' Importantly, the exposure and outcome standard errors supplied to this
#' function must already have been calibrated for sample structure using
#' LDSC. Specifically, if the original GWAS standard errors are
#' \eqn{\sigma_{\hat\gamma_j}} and \eqn{\sigma_{\hat\Gamma_j}}, the
#' inputs should be
#'
#' \deqn{
#' \mathrm{gamma\_se}_j =
#' \sqrt{c_1}\sigma_{\hat\gamma_j}
#' }
#'
#' and
#'
#' \deqn{
#' \mathrm{Gamma\_se}_j =
#' \sqrt{c_2}\sigma_{\hat\Gamma_j}.
#' }
#'
#' Accordingly, the correlation parameter supplied through \code{rho}
#' must be
#'
#' \deqn{
#' \rho = \frac{c_{12}}{\sqrt{c_1 c_2}},
#' }
#'
#' where \eqn{c_1}, \eqn{c_2}, and \eqn{c_{12}} are the LDSC estimates
#' characterizing the exposure GWAS, outcome GWAS, and their cross-trait
#' sample structure, respectively.
#'
#' @param gamma_hat A numeric vector of estimated SNP-exposure
#'   associations, \eqn{\hat\gamma_j}.
#'
#' @param gamma_se A numeric vector of standard errors for
#'   \code{gamma_hat}. These must already be LDSC-calibrated standard
#'   errors, i.e., the original exposure-GWAS standard errors multiplied
#'   by \eqn{\sqrt{c_1}}.
#'
#' @param Gamma_hat A numeric vector of estimated SNP-outcome
#'   associations, \eqn{\hat\Gamma_j}.
#'
#' @param Gamma_se A numeric vector of standard errors for
#'   \code{Gamma_hat}. These must already be LDSC-calibrated standard
#'   errors, i.e., the original outcome-GWAS standard errors multiplied
#'   by \eqn{\sqrt{c_2}}.
#'
#' @param rho A scalar correlation parameter describing the dependence
#'   between the calibrated exposure- and outcome-GWAS estimation errors.
#'   It should be calculated as
#'   \eqn{\rho = c_{12}/\sqrt{c_1c_2}} and must lie strictly between
#'   -1 and 1.
#'
#' @param eta A positive scalar giving the standard deviation of the
#'   randomization variable used for randomized IV selection. This value
#'   must be identical to that used when selecting the input SNPs.
#'   Default is \code{0.5}.
#'
#' @param lambda A non-negative scalar giving the randomized IV-selection
#'   threshold. This value must be identical to that used when selecting
#'   the input SNPs. Default is \code{4.06}.
#'
#' @param sig.level Confidence level for the BRIVW confidence interval.
#'   Default is \code{0.95}.
#'
#' @details
#' BRIVW applies Rao--Blackwell correction to both the exposure and
#' outcome GWAS associations and accounts for the covariance induced by
#' sample structure.
#'
#' The BMEI-C test combines two BMEI statistics obtained under two
#' allele orientations. The first uses the input harmonized coding.
#' The second simultaneously flips the exposure and outcome associations
#' for SNPs with negative estimated exposure associations, so that all
#' exposure associations are oriented in the positive direction.
#' The combined P-value is based on the maximum absolute Z statistic
#' while accounting for the estimated correlation between the two
#' statistics.
#'
#' @return
#' A one-row data frame containing:
#'
#' \describe{
#'   \item{\code{estimate}}{BRIVW causal-effect estimate.}
#'   \item{\code{std_error}}{Standard error of the BRIVW estimate.}
#'   \item{\code{ci_lower}}{Lower confidence limit.}
#'   \item{\code{ci_upper}}{Upper confidence limit.}
#'   \item{\code{p_value}}{Two-sided P-value for the BRIVW causal effect.}
#'   \item{\code{n_snp}}{Number of SNPs included in the analysis.}
#'   \item{\code{bmei_p_value}}{P-value from the combined BMEI-C test.}
#'   \item{\code{rho}}{Correlation parameter used in the analysis.}
#'   \item{\code{eta}}{Randomization standard deviation.}
#'   \item{\code{lambda}}{Randomized IV-selection threshold.}
#'   \item{\code{sig_level}}{Confidence level used for the confidence interval.}
#' }
#'
#' @importFrom stats dnorm pnorm qnorm
#' @importFrom mvtnorm pmvnorm
#'
#' @examples
#' \dontrun{
#' fit <- brivw(
#'   gamma_hat = dat$beta.exposure,
#'   gamma_se  = dat$se.exposure,
#'   Gamma_hat = dat$beta.outcome,
#'   Gamma_se  = dat$se.outcome,
#'   rho = 0.10,
#'   eta = 0.5,
#'   lambda = 4.06,
#'   sig.level = 0.95
#' )
#'
#' fit
#' }
#'
#' @export
brivw <- function(gamma_hat,
                  gamma_se,
                  Gamma_hat,
                  Gamma_se,
                  rho,
                  eta = 0.5,
                  lambda = 4.06,
                  sig.level = 0.95) {
  
  # --------------------------------------------------------------------
  # 1. Input checks
  # --------------------------------------------------------------------
  
  x <- list(
    gamma_hat = gamma_hat,
    gamma_se = gamma_se,
    Gamma_hat = Gamma_hat,
    Gamma_se = Gamma_se
  )
  
  if (!all(vapply(x, is.numeric, logical(1)))) {
    stop("gamma_hat, gamma_se, Gamma_hat, and Gamma_se must be numeric.")
  }
  
  n_snp <- length(gamma_hat)
  
  if (n_snp < 2L) {
    stop("At least two SNPs are required.")
  }
  
  if (any(vapply(x, length, integer(1)) != n_snp)) {
    stop(
      "gamma_hat, gamma_se, Gamma_hat, and Gamma_se ",
      "must have the same length."
    )
  }
  
  if (any(!is.finite(gamma_hat)) ||
      any(!is.finite(gamma_se)) ||
      any(!is.finite(Gamma_hat)) ||
      any(!is.finite(Gamma_se))) {
    stop("All association estimates and standard errors must be finite.")
  }
  
  if (any(gamma_se <= 0) || any(Gamma_se <= 0)) {
    stop("All standard errors must be positive.")
  }
  
  if (length(rho) != 1L ||
      !is.finite(rho) ||
      rho <= -1 ||
      rho >= 1) {
    stop("rho must be a finite scalar strictly between -1 and 1.")
  }
  
  if (length(eta) != 1L ||
      !is.finite(eta) ||
      eta <= 0) {
    stop("eta must be a finite positive scalar.")
  }
  
  if (length(lambda) != 1L ||
      !is.finite(lambda) ||
      lambda < 0) {
    stop("lambda must be a finite non-negative scalar.")
  }
  
  if (length(sig.level) != 1L ||
      !is.finite(sig.level) ||
      sig.level <= 0 ||
      sig.level >= 1) {
    stop("sig.level must be strictly between 0 and 1.")
  }
  
  
  # --------------------------------------------------------------------
  # 2. Local function for Rao--Blackwell correction
  #
  # This is defined locally and is not exposed to package users.
  # --------------------------------------------------------------------
  
  rb_correct <- function(gamma,
                         Gamma,
                         gamma_se,
                         Gamma_se) {
    
    A_plus <-
      lambda / eta -
      gamma / (gamma_se * eta)
    
    A_minus <-
      -lambda / eta -
      gamma / (gamma_se * eta)
    
    # Selection probability term.
    # The upper-tail form is used for improved numerical stability.
    denominator <-
      stats::pnorm(A_plus, lower.tail = FALSE) +
      stats::pnorm(A_minus)
    
    if (any(!is.finite(denominator)) ||
        any(denominator <= 0)) {
      stop("Numerical failure in the Rao--Blackwell correction.")
    }
    
    phi_diff <-
      stats::dnorm(A_plus) -
      stats::dnorm(A_minus)
    
    A3 <-
      (
        A_plus * stats::dnorm(A_plus) -
          A_minus * stats::dnorm(A_minus)
      ) / denominator
    
    A4 <-
      (phi_diff / denominator)^2
    
    # Rao--Blackwellized exposure association
    gamma_rb <-
      gamma -
      gamma_se / eta *
      phi_diff / denominator
    
    # Conditional variance of the RB exposure association
    var_gamma_rb <-
      gamma_se^2 *
      (
        1 -
          A3 / eta^2 +
          A4 / eta^2
      )
    
    # Rao--Blackwellized outcome association
    Gamma_rb <-
      Gamma -
      rho * Gamma_se / eta *
      phi_diff / denominator
    
    # Conditional covariance between the RB associations
    cov_rb <-
      rho * gamma_se * Gamma_se *
      (
        1 -
          A3 / eta^2 +
          A4 / eta^2
      )
    
    list(
      gamma_rb = gamma_rb,
      Gamma_rb = Gamma_rb,
      var_gamma_rb = var_gamma_rb,
      cov_rb = cov_rb
    )
  }
  
  
  # --------------------------------------------------------------------
  # 3. BRIVW estimation
  # --------------------------------------------------------------------
  
  rb <- rb_correct(
    gamma = gamma_hat,
    Gamma = Gamma_hat,
    gamma_se = gamma_se,
    Gamma_se = Gamma_se
  )
  
  gamma_rb <- rb$gamma_rb
  Gamma_rb <- rb$Gamma_rb
  var_gamma_rb <- rb$var_gamma_rb
  cov_rb <- rb$cov_rb
  
  theta1 <-
    sum(
      (Gamma_rb * gamma_rb - cov_rb) /
        Gamma_se^2
    )
  
  theta2 <-
    sum(
      (gamma_rb^2 - var_gamma_rb) /
        Gamma_se^2
    )
  
  if (!is.finite(theta2) ||
      abs(theta2) < sqrt(.Machine$double.eps)) {
    stop("The BRIVW denominator is zero or numerically unstable.")
  }
  
  beta_hat <- theta1 / theta2
  
  residual <-
    Gamma_rb * gamma_rb -
    cov_rb -
    beta_hat *
    (gamma_rb^2 - var_gamma_rb)
  
  beta_var <-
    sum(
      residual^2 /
        Gamma_se^4
    ) /
    theta2^2
  
  if (!is.finite(beta_var) || beta_var < 0) {
    stop("The estimated BRIVW variance is invalid.")
  }
  
  beta_se <- sqrt(beta_var)
  
  alpha <- 1 - sig.level
  
  critical_value <-
    stats::qnorm(1 - alpha / 2)
  
  ci_lower <-
    beta_hat -
    critical_value * beta_se
  
  ci_upper <-
    beta_hat +
    critical_value * beta_se
  
  if (beta_se > 0) {
    
    beta_p <-
      2 *
      stats::pnorm(
        abs(beta_hat / beta_se),
        lower.tail = FALSE
      )
    
  } else {
    
    beta_p <- NA_real_
    
  }
  
  
  # --------------------------------------------------------------------
  # 4. Local function for a BMEI statistic under one allele coding
  #
  # This function is internal to brivw() and is not exposed to users.
  # --------------------------------------------------------------------
  
  bmei_statistic <- function(gamma,
                             Gamma) {
    
    rb_bmei <- rb_correct(
      gamma = gamma,
      Gamma = Gamma,
      gamma_se = gamma_se,
      Gamma_se = Gamma_se
    )
    
    g_rb <- rb_bmei$gamma_rb
    G_rb <- rb_bmei$Gamma_rb
    g_var_rb <- rb_bmei$var_gamma_rb
    gG_cov_rb <- rb_bmei$cov_rb
    
    theta1_b <-
      sum(
        (G_rb * g_rb - gG_cov_rb) /
          Gamma_se^2
      )
    
    theta2_b <-
      sum(
        (g_rb^2 - g_var_rb) /
          Gamma_se^2
      )
    
    if (!is.finite(theta2_b) ||
        abs(theta2_b) < sqrt(.Machine$double.eps)) {
      stop("The BMEI denominator is zero or numerically unstable.")
    }
    
    beta_b <- theta1_b / theta2_b
    
    T_hat <-
      sum(g_rb / Gamma_se^2)
    
    # Corrected BMEI numerator
    lambda_terms <-
      theta2_b *
      G_rb / Gamma_se^2 -
      theta1_b *
      g_rb / Gamma_se^2 +
      G_rb *
      g_var_rb / Gamma_se^4 -
      g_rb *
      gG_cov_rb / Gamma_se^4
    
    lambda_b <-
      sum(lambda_terms)
    
    # Leading-term variance estimator
    xi_hat <-
      g_rb * G_rb -
      gG_cov_rb -
      beta_b *
      (g_rb^2 - g_var_rb)
    
    omega_hat <-
      G_rb -
      beta_b * g_rb
    
    u_hat <-
      theta2_b *
      omega_hat / Gamma_se^2 -
      T_hat *
      xi_hat / Gamma_se^2
    
    var_lambda <-
      sum(u_hat^2)
    
    if (!is.finite(var_lambda) ||
        var_lambda <= 0) {
      stop("The estimated BMEI variance is non-positive or unstable.")
    }
    
    z_value <-
      lambda_b /
      sqrt(var_lambda)
    
    list(
      z = z_value,
      u = u_hat
    )
  }
  
  
  # --------------------------------------------------------------------
  # 5. BMEI under the input harmonized allele coding
  # --------------------------------------------------------------------
  
  bmei_1 <-
    bmei_statistic(
      gamma = gamma_hat,
      Gamma = Gamma_hat
    )
  
  
  # --------------------------------------------------------------------
  # 6. BMEI under exposure-oriented allele coding
  #
  # Exposure and outcome associations are simultaneously flipped
  # whenever gamma_hat < 0.
  # --------------------------------------------------------------------
  
  gamma_flip <- gamma_hat
  Gamma_flip <- Gamma_hat
  
  flip <- gamma_flip < 0
  
  gamma_flip[flip] <-
    -gamma_flip[flip]
  
  Gamma_flip[flip] <-
    -Gamma_flip[flip]
  
  bmei_2 <-
    bmei_statistic(
      gamma = gamma_flip,
      Gamma = Gamma_flip
    )
  
  
  # --------------------------------------------------------------------
  # 7. Combined BMEI-C test
  # --------------------------------------------------------------------
  
  max_abs_z <-
    max(
      abs(bmei_1$z),
      abs(bmei_2$z)
    )
  
  rho_bmei <-
    sum(
      bmei_1$u *
        bmei_2$u
    ) /
    sqrt(
      sum(bmei_1$u^2) *
        sum(bmei_2$u^2)
    )
  
  if (!is.finite(rho_bmei)) {
    stop(
      "Unable to estimate the correlation between ",
      "the two BMEI statistics."
    )
  }
  
  # Numerical safeguard for the correlation matrix
  rho_bmei <-
    max(
      min(rho_bmei, 1 - 1e-10),
      -1 + 1e-10
    )
  
  correlation_matrix <-
    matrix(
      c(
        1, rho_bmei,
        rho_bmei, 1
      ),
      nrow = 2
    )
  
  bmei_p <-
    1 -
    as.numeric(
      mvtnorm::pmvnorm(
        lower = rep(-max_abs_z, 2),
        upper = rep(max_abs_z, 2),
        mean = c(0, 0),
        corr = correlation_matrix
      )
    )
  
  # Protect against numerical values slightly outside [0, 1]
  bmei_p <-
    min(
      max(bmei_p, 0),
      1
    )
  
  
  # --------------------------------------------------------------------
  # 8. Output
  # --------------------------------------------------------------------
  
  result <-
    data.frame(
      estimate = beta_hat,
      std_error = beta_se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = beta_p,
      n_snp = n_snp,
      bmei_p_value = bmei_p,
      rho = rho,
      eta = eta,
      lambda = lambda,
      sig_level = sig.level,
      row.names = NULL
    )
  
  return(result)
}
