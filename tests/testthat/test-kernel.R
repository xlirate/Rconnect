test_that("normalize_kernal normalizes to 1 if sum(k) != 0", {
  expect_equal(sum(normalize_kernal(matrix(c(2)))), 1)
  expect_equal(sum(normalize_kernal(matrix(c(1)))), 1)
  expect_equal(sum(normalize_kernal(matrix(c(-1)))), 1)
  expect_equal(sum(normalize_kernal(matrix(c(-2)))), 1)

  expect_equal(sum(normalize_kernal(matrix(c(0,1,2,3,4), 5, 5))), 1)
  expect_equal(sum(normalize_kernal(matrix(c(-2,-1,0,1), 4, 4))), 1)
})

test_that("normalize_kernal normalizes to 0 if sum(k) == 0", {
  expect_equal(sum(normalize_kernal(matrix(c(0)))), 0)
  expect_equal(sum(normalize_kernal(matrix(c(0,0,0), 3, 3))), 0)
  expect_equal(sum(normalize_kernal(matrix(c(-2, -1, 0, 1, 2), 5, 5))), 0)
})

test_that("hard_circle_kernal returns a matrix,
before normalizing,
with all cells who's center is within 'r'
of the center of the matrix (inclusive)
set to 1 and otherwise set to 0", {
  helper_f <- function(r){
    k <- hard_uniform_circle_kernal(r)
    x <- matrix(-floor(nrow(k)/2):floor(nrow(k)/2), nrow(k), ncol(k), byrow=TRUE)
    y <- matrix(-floor(ncol(k)/2):floor(ncol(k)/2), nrow(k), ncol(k))
    !(FALSE %in% mapply(function(r, k, x, y){if(r*r >= x*x + y*y){k==1}else{k==0}}, r, k, x, y))
  }
  expect_true(helper_f(0))
  expect_true(helper_f(1))
  expect_true(helper_f(2))
  expect_true(helper_f(3))
  expect_true(helper_f(10))
  expect_true(helper_f(100))
  expect_true(helper_f(0.8))
  expect_true(helper_f(1.5))
  expect_true(helper_f(2.3))
})

test_that("smooth_uniform_circle_kernal,
before normalizing,
creates a matrix who's sum ==
the area of a circle with the same r", {
  expect_equal(sum(smooth_uniform_circle_kernal(0)), 0)
  expect_equal(sum(smooth_uniform_circle_kernal(0.1)), pi*0.1^2)
  expect_equal(sum(smooth_uniform_circle_kernal(0.5)), pi*0.5^2)
  expect_equal(sum(smooth_uniform_circle_kernal(0.9)), pi*0.9^2)
  expect_equal(sum(smooth_uniform_circle_kernal(1.0)), pi*1^2)
  expect_equal(sum(smooth_uniform_circle_kernal(1.1)), pi*1.1^2)
  expect_equal(sum(smooth_uniform_circle_kernal(1.5)), pi*1.5^2)
  expect_equal(sum(smooth_uniform_circle_kernal(1.7)), pi*1.7^2)
  expect_equal(sum(smooth_uniform_circle_kernal(2)), pi*2^2)
  expect_equal(sum(smooth_uniform_circle_kernal(2.2)), pi*2.2^2)
  expect_equal(sum(smooth_uniform_circle_kernal(3.3)), pi*3.3^2)
  expect_equal(sum(smooth_uniform_circle_kernal(4.7)), pi*4.7^2)
  expect_equal(sum(smooth_uniform_circle_kernal(10.11)), pi*10.11^2)
  expect_equal(sum(smooth_uniform_circle_kernal(55.67)), pi*55.67^2)
  expect_equal(sum(smooth_uniform_circle_kernal(101)), pi*101^2)
})


test_that("discance_kernal creates a kernal where the center == 0
and each other cell contains it's euclidian distance from the center", {
  helper_f <- function(r){
    k <- distance_kernal(r)
    xoff <- ceiling(ncol(k)/2)
    yoff <- ceiling(nrow(k)/2)
    for(x in 1:ncol(k)){
      for(y in 1:nrow(k)){
        if(sqrt((x-xoff)^2 + (y-yoff)^2) != k[x,y]){
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }
  expect_true(helper_f(1))
  expect_true(helper_f(2))
  expect_true(helper_f(3))
  expect_true(helper_f(10))
  expect_true(helper_f(100))
  expect_true(helper_f(0.8))
  expect_true(helper_f(1.5))
  expect_true(helper_f(2.3))
})



test_that("exponential_kernel creates a kernal with all non negative values", {
  helper_f <- function(r){
    k <- exponential_kernel(r)
    xoff <- ceiling(ncol(k)/2)
    yoff <- ceiling(nrow(k)/2)
    for(x in 1:ncol(k)){
      for(y in 1:nrow(k)){
        if(sqrt((x-xoff)^2 + (y-yoff)^2) != k[x,y]){
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }
  expect_true(helper_f(1))
  expect_true(helper_f(2))
  expect_true(helper_f(3))
  expect_true(helper_f(10))
  expect_true(helper_f(100))
  expect_true(helper_f(0.8))
  expect_true(helper_f(1.5))
  expect_true(helper_f(2.3))
})

