test_that("verifying data in matrix is reasonable", {
  epsilon = 0.0005 #created epsilon as a global constant
  #finds if each element of the matrix is 0 or 1
  #and finds if the matrix total is greater than or equal to 0 
  #and less than or equal to the number of elements in the matrix
  in_range <- function(matrix){ 
    for(row in 1:nrow(matrix)){
      for(col in 1:ncol(matrix)){
        if (matrix[row,col]>0.5){ #ensures that all elements above 0.5 are equal to 1
          expect_equal(matrix[row,col], 1)
        } else{ #ensures all elements less than or equal to 0.5 are equal to 0
          expect_equal(matrix[row,col], 0)
        }
      }
    }
    matrix_ratio <- sum(matrix)/(dim(matrix)[1]*dim(matrix)[2]) #finds the decimal of 1's amoung the matrix
    expect_lte(matrix_ratio, 1+epsilon) #sees if decimal is less than or equal to 1 plus episilon (for significant digits issues)
    expect_gte(matrix_ratio, 0) #sees if decimal is greater than or equal to 0
  }
  matrix <- simulate_varmat()
  i = 1
  while (i<20){ #iterations through in_range function 20 times
    in_range(matrix)
    i = i + 1
  }
})