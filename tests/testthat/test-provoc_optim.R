test_that("Ensure estimate is within the bootstrap confidence interval.", {
  expect_results <- function(converg_info){
    i = 1
    while(i <= nrow(converg_info$res_df)){ #iterates through each row of res_df
      expect_lt(converg_info$res_df$rho[i], converg_info$res_df$ci_high[i]) #checks to see if rho is less than the upper bound of the CI
      expect_gt(converg_info$res_df$rho[i], converg_info$res_df$ci_low[i]) #checks to see if rho is less than the lower bound of the CI
      i = i + 1
    }
  }
  test_bootstrap <- function(rel_counts){
    varmat <- simulate_varmat() #simulates a new variant matrix
    coco <- simulate_coco(varmat, rel_counts)
    converg_info <- provoc_optim(coco, varmat, 1000, TRUE) #estimates the proportions of VOCs with 1000 bootstrap samples
    expect_results(converg_info)
  }
  j = 0
  while(j < 10){ #simulates data 10 different times with a relative counts of each VOC (60 unit tests total)
    test_bootstrap(sample.int(1000,3,replace=TRUE)) #takes a random array of integers up to 1000 of length 3
    j = j + 1
  }
})

test_that("output is of correct type",{
  varmat <- simulate_varmat()
  coco <- simulate_coco(varmat)
  res <- provoc_optim(coco, varmat, 0, FALSE)
  expect_type(res$res_df, "data.frame")
  expect_type(res$convergence, "logical")
})