
.beyer_fix_data <- function(data){
  if(!is.matrix(data)){
    errorCondition("Data must be a martix")
  }
  return(data)
}

.beyer_alg <- function(alg, data, z = 0.5, beta = 0.2, r0=0.05){
  return(alg(.beyer_fix_data(data)^z, flatten_kernel(normalize_kernel(beyer_kernel(beta, r0)))))
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
beyer_stretch <- function(data, kernel, beta = 0.2, z = 0.5, threshold=NULL){
  return(.beyer_alg(convolve_stretch, data, z, beta, r0))
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
beyer_wrap <- function(data, beta = 0.2, z = 0.5, r0=0.05){
  return(.beyer_alg(convolve_wrap, data, z, beta, r0))
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
beyer_reflect <- function(data, beta = 0.2, z = 0.5, r0=0.05){
  return(.beyer_alg(convolve_refect, data, z, beta, r0))
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
beyer_zero <- function(data, beta = 0.2, z = 0.5, r0=0.05){
  return(.beyer_alg(convolve_zero, data, z, beta, r0))
}

#
# NaN NaN|NaN NaN NaN NaN NaN|NaN NaN
# NaN NaN|NaN NaN NaN NaN NaN|NaN NaN
# -------+-------------------+-------
# NaN NaN| A   B   C   D   E |NaN NaN
# NaN NaN| F   G   H   I   J |NaN NaN
# NaN NaN| K   L   M   N   O |NaN NaN
# NaN NaN| P   Q   R   S   T |NaN NaN
# NaN NaN| U   V   W   X   Y |NaN NaN
# -------+-------------------+-------
# NaN NaN|NaN NaN NaN NaN NaN|NaN NaN
# NaN NaN|NaN NaN NaN NaN NaN|NaN NaN
#
#' @export
convolve_nan <- function(data, beta = 0.2, z = 0.5, r0=0.05){
  return(.convolve_alg(convolve_nan, data, z, beta, r0))
}

# the output is shrunk down by enough that it never reaches outside the data matrix in the first place
#' @export
beyer_shrink <- function(data, beta = 0.2, z = 0.5, r0=0.05){
  return(.beyer_alg(convolve_shrink, data, z, beta, r0))
}
