test_that("Ensure estimate is within the bootstrap confidence interval.", {
  expect_results <- function(converg_info){
    i = 1
    while(i <= nrow(converg_info)){ #iterates through each row of res_df
      expect_lt(converg_info$rho[i], converg_info$ci_high[i]) #checks to see if rho is less than the upper bound of the CI
      expect_gt(converg_info$rho[i], converg_info$ci_low[i]) #checks to see if rho is less than the lower bound of the CI
      i = i + 1
    }
  }
  test_bootstrap <- function(rel_counts){
    varmat <- simulate_varmat() #simulates a new variant matrix
    coco <- simulate_coco(varmat, rel_counts)
    converg_info <- provoc(coco, varmat, NULL, 1, 1000, 20, TRUE) #estimates the proportions of VOCs with 1000 bootstrap samples
    expect_results(converg_info)
  }
  j = 0
  while(j < 10){ #simulates data 10 different times with a relative counts of each VOC (60 unit tests total)
    test_bootstrap(sample.int(1000,3,replace=TRUE)) #takes a random array of integers up to 1000 of length 3
    j = j + 1
  }
})

test_that("rho total equals 1", {
  rho_test <- function(converg_info){
    epsilon = 0.005
    i = 1
    rho_total = 0
    while(i <= nrow(converg_info)){
      expect_lt(converg_info$rho[i], 1)
      expect_gt(converg_info$rho[i], 0)
      rho_total = rho_total + converg_info$rho[i]
      i = i + 1
    }
    expect_lt(rho_total, 1+epsilon)
    expect_gt(rho_total, 1-epsilon)
  }
  create_provoc_table <- function(rel_counts){
    varmat <- simulate_varmat() #simulates a new variant matrix
    coco <- simulate_coco(varmat, rel_counts)
    converg_info <- provoc(coco, varmat, NULL, 1, 0, 20, TRUE) #don't need bootstrap samples for this unit test 
    rho_test(converg_info)
  }
  j = 0
  while (j < 10){
    create_provoc_table(sample.int(1000,3,replace=TRUE))
    j = j+1
  }
})