#test-svyIncProp.R
testthat::test_that("SvyIncProp runs correctly", {
  # Load necessary data
  load(system.file("data", "AlcanLink22.rda", package = "alcanR"))
  t1 <- svydesign(id=~1,strata=~SUD_NEST,fpc=~fpc,weights=~WTANAL,data=AlcanLink22)
  #source("C:\\Users\\jwparrish\\Documents\\alcanR\\R\\SvyIncProp.R")

  # Run the function
  frpt_IP <- SvyIncProp(
    SurveyDesign = t1,
    event = "FST_RPT_Indicator",
    eventime = "FST_RPT_CensorAge",
    nreps = 20
  )

  # Check if the output is a data frame
  testthat::expect_true(is.data.frame(frpt_IP), info = "Output should be a data frame")

  # Additional checks: Ensure specific columns exist
  testthat::expect_true(all(c("theta", "empirical.mean") %in% names(frpt_IP)),
                        info = "Output should contain expected columns")

  # Check the number of rows (if expected)
  testthat::expect_gt(nrow(frpt_IP), 0, info = "Data frame should not be empty")

})
