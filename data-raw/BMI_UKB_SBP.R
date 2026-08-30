# Prepare BMI--SBP example data ------------------------------------------

dat_raw <- read.delim(
  "path/to/BMI_UKB_SBP.txt"
)

BMI_UKB_SBP <- data.frame(
  gamma_hat = dat_raw$beta.exposure,
  gamma_se  = dat_raw$se.exposure,
  Gamma_hat = dat_raw$beta.outcome,
  Gamma_se  = dat_raw$se.outcome
)

# Check the resulting dataset
str(BMI_UKB_SBP)
summary(BMI_UKB_SBP)

# Save as package data
usethis::use_data(
  BMI_UKB_SBP,
  overwrite = TRUE
)


# Prepare LDSC C matrix --------------------------------------------------

C_raw <- read.delim(
  "path/to/BMI_UKB~SBP_C.txt",
  header = FALSE
)

# The original C matrix was saved using as.vector(C),
# which follows column-major ordering in R.
C_vec <- as.numeric(C_raw[1, ])

BMI_UKB_SBP_C <- matrix(
  C_vec,
  nrow = 2,
  ncol = 2,
  byrow = FALSE,
  dimnames = list(
    c("Exposure", "Outcome"),
    c("Exposure", "Outcome")
  )
)

BMI_UKB_SBP_C

usethis::use_data(
  BMI_UKB_SBP_C,
  overwrite = TRUE
)
