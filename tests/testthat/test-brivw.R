test_that("BRIVW BMI-SBP example returns valid results", {

  rho <- BMI_UKB_SBP_C[1, 2] /
    sqrt(
      BMI_UKB_SBP_C[1, 1] *
        BMI_UKB_SBP_C[2, 2]
    )

  fit <- brivw(
    gamma_hat = BMI_UKB_SBP$gamma_hat,
    gamma_se = BMI_UKB_SBP$gamma_se,
    Gamma_hat = BMI_UKB_SBP$Gamma_hat,
    Gamma_se = BMI_UKB_SBP$Gamma_se,
    rho = rho,
    eta = 0.5,
    lambda = 4.06,
    sig.level = 0.95
  )

  expect_true(is.data.frame(fit))
  expect_equal(nrow(fit), 1L)

  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$std_error))
  expect_true(is.finite(fit$bmei_p_value))

  expect_gt(fit$std_error, 0)

  expect_gte(fit$p_value, 0)
  expect_lte(fit$p_value, 1)

  expect_gte(fit$bmei_p_value, 0)
  expect_lte(fit$bmei_p_value, 1)

  expect_equal(
    fit$n_snp,
    nrow(BMI_UKB_SBP)
  )
})
