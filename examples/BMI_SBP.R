# Run from the package root:
# source("examples/BMI_SBP.R")
# This uses the local source code, rather than an older installed package.

source("R/brivw.R")
load("data/BMI_UKB_SBP.rda")
load("data/BMI_UKB_SBP_C.rda")

# The example standard errors are already LDSC-calibrated.
# Use the exact two-sided IV-selection threshold.
lambda <- qnorm(5e-5 / 2, lower.tail = FALSE)
rho <- BMI_UKB_SBP_C[1, 2] /
  sqrt(BMI_UKB_SBP_C[1, 1] * BMI_UKB_SBP_C[2, 2])

fit <- brivw(
  gamma_hat = BMI_UKB_SBP$gamma_hat,
  gamma_se = BMI_UKB_SBP$gamma_se,
  Gamma_hat = BMI_UKB_SBP$Gamma_hat,
  Gamma_se = BMI_UKB_SBP$Gamma_se,
  rho = rho,
  eta = 0.5,
  lambda = lambda,
  sig.level = 0.95
)

print(fit)
