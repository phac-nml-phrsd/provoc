
test_that("rho total is less than 1", {
  epsilon = 0.005 #initializes epsilon, put outside of rho_test function since it doesn't change
  rho_test <- function(converg_info){
    i = 1
    rho_total = 0
    while(i <= nrow(converg_info)){
      expect_lt(converg_info$rho[i], 1) #checks if each rho is less than 1
      expect_gt(converg_info$rho[i], 0) #checks if each rho is greater than 1
      rho_total = rho_total + converg_info$rho[i]
      i = i + 1
    }
    expect_lt(rho_total, 1+epsilon) #checks to see if rho total is between 1-epsilon and 1+epsilon, since addition with significant digits can cause issues on preciseness.
  }
  create_provoc_table <- function(rel_counts){
    varmat <- simulate_varmat() #simulates a new variant matrix
    coco <- simulate_coco(varmat, rel_counts)
    count <- coco$count
    coverage <- coco$coverage
    data <- as.data.frame(Baaijens)
    data$mutation <- parse_mutations(data$label)
    converg_info <- provoc(cbind(count, coverage) ~ B.1.1.7 + B.1.617.2, data = data, by = "sample_name")
    rho_test(converg_info)
  }
  j = 0
  while (j < 10){
    create_provoc_table(sample.int(1000,3,replace=TRUE))
    j = j+1
  }
})