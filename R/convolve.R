
convolve_fix_data <- function(data){
  if(!is.matrix(data)){
    errorCondition("Data must be a martix")
  }
  return(data)
}

convolve_fix_kernel <- function(kernel, normalize){
  if(!is.list(kernel)){
    kernel <- list(kernel)
  }
  for(k in kernel){
    if(!(ncol(k)%%2 && nrow(k)%%2)){
      errorCondition("Kernels must be an n by m martix where both n and m are odd")
    }
  }
  if(normalize){
    return(normalize_kernel(kernel))
  }else{
    return(kernel)
  }
}

convolve_fix_reps <- function(reps){
  if(reps < 0){
    errorCondition("You cannot have negitive repititions")
  }else if(length(reps) != 1){
    errorCondition("Repititons must have length 1")
  }else{
    return(reps)
  }
}

convolve_alg <- function(alg, data, kernel, reps = 1L, normalize = TRUE){
  data <- convolve_fix_data(data)
  kernel <- convolve_fix_kernel(kernel, normalize)
  reps <- convolve_fix_reps(reps)
  while(reps > 0){
    reps <- reps-1
    data <- alg(data, kernel)
  }
  return(data)
}


#
# a a|a b c d e|e e
# a a|a b c d e|e e
# ---+---------+---
# a a|A B C D E|e e
# f f|F G H I J|j j
# k k|K L M N O|o o
# p p|P Q R S T|t t
# u u|U V W X Y|y y
# ---+---------+---
# u u|u v w x y|y y
# u u|u v w x y|y y
#
#' @export
convolve_stretch <- function(data, kernel, reps = 1L, normalize = TRUE){
  return(convolve_alg(Rconnect:::.convolve_stretch, data, kernel, reps, normalize))
}

#
# s t|p q r s t|p q
# x y|u v w x y|u v
# ---+---------+---
# d e|A B C D E|a b
# i j|F G H I J|f g
# n o|K L M N O|k l
# s t|P Q R S T|p q
# x y|U V W X Y|u v
# ---+---------+---
# d e|a b c d e|a b
# i j|f g h i j|f g
#
#' @export
convolve_wrap <- function(data, kernel, reps = 1L, normalize = TRUE){
  return(convolve_alg(Rconnect:::.convolve_wrap, data, kernel, reps, normalize))
}
#
# h f|f g h i j|j i
# b a|a b c d e|e d
# ---+---------+---
# b a|A B C D E|e d
# g f|F G H I J|j i
# l k|K L M N O|o n
# q p|P Q R S T|t s
# v u|U V W X Y|y x
# ---+---------+---
# v u|u v w x y|y x
# q u|p q r s t|t s
#
#' @export
convolve_reflect <- function(data, kernel, reps = 1L, normalize = TRUE){
  return(convolve_alg(Rconnect:::.convolve_refect, data, kernel, reps, normalize))
}
#
# 0 0|0 0 0 0 0|0 0
# 0 0|0 0 0 0 0|0 0
# ---+---------+---
# 0 0|A B C D E|0 0
# 0 0|F G H I J|0 0
# 0 0|K L M N O|0 0
# 0 0|P Q R S T|0 0
# 0 0|U V W X Y|0 0
# ---+---------+---
# 0 0|0 0 0 0 0|0 0
# 0 0|0 0 0 0 0|0 0
#
#' @export
convolve_zero <- function(data, kernel, reps = 1L, normalize = TRUE){
  return(convolve_alg(Rconnect:::.convolve_zero, data, kernel, reps, normalize))
}

# the output is shrunk down by enough that it never reaches outside the data matrix in the first place
#' @export
convolve_shrink <- function(data, kernel, reps = 1L, normalize = TRUE){
  return(convolve_alg(Rconnect:::.convolve_shrink, data, kernel, reps, normalize))
}
