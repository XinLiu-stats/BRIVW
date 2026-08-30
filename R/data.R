#' BMI--SBP Mendelian Randomization Example Data
#'
#' Analysis-ready GWAS summary statistics for Mendelian randomization
#' analysis of body mass index (BMI) on systolic blood pressure (SBP).
#'
#' The SNPs in this dataset have undergone the preprocessing steps used
#' in the accompanying BRIVW analysis, including allele harmonization,
#' randomized instrument selection, and LD pruning.
#'
#' The standard errors have already been calibrated for sample structure
#' using LDSC. Specifically, \code{gamma_se} and \code{Gamma_se} equal
#' the original exposure- and outcome-GWAS standard errors multiplied by
#' \eqn{\sqrt{c_1}} and \eqn{\sqrt{c_2}}, respectively.
#'
#' @format A data frame with one row per SNP and four variables:
#' \describe{
#'   \item{\code{gamma_hat}}{
#'     Estimated SNP-exposure association.
#'   }
#'   \item{\code{gamma_se}}{
#'     LDSC-calibrated standard error of the SNP-exposure association.
#'   }
#'   \item{\code{Gamma_hat}}{
#'     Estimated SNP-outcome association.
#'   }
#'   \item{\code{Gamma_se}}{
#'     LDSC-calibrated standard error of the SNP-outcome association.
#'   }
#' }
#'
#' @source GWAS summary statistics used in the accompanying BRIVW analysis.
#'
#' @name BMI_UKB_SBP
#' @docType data
NULL


#' LDSC C Matrix for the BMI--SBP Example
#'
#' LDSC-estimated matrix characterizing the sample structure of the
#' exposure and outcome GWASs used in the \code{BMI_UKB_SBP} example.
#'
#' The diagonal elements correspond to \eqn{c_1} and \eqn{c_2}, and the
#' off-diagonal elements correspond to \eqn{c_{12}}. The correlation
#' parameter supplied to \code{brivw()} is therefore calculated as
#'
#' \deqn{
#' \rho = \frac{c_{12}}{\sqrt{c_1c_2}}.
#' }
#'
#' The standard errors in \code{BMI_UKB_SBP} have already been calibrated
#' using \eqn{c_1} and \eqn{c_2}.
#'
#' @format A 2 by 2 numeric matrix with row and column names
#' \code{"Exposure"} and \code{"Outcome"}.
#'
#' @name BMI_UKB_SBP_C
#' @docType data
NULL