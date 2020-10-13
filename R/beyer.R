
beyer_fix_data <- function(data){
  if(!is.matrix(data)){
    errorCondition("Data must be a martix")
  }
  return(data)
}

beyer_fix_kernel <- function(kernel, threshold=NULL){
  if(is.list(kernel)){
    kernel <- Reduce(function(k1, k2){
      k <- matrix(0, nrow=(nrow(k1)+nrow(k2)-1), ncol=(nrow(k1)+nrow(k2)-1))
      for(x1 in 1:ncol(k1)){
        for(y1 in 1:nrow(k1)){
          for(x2 in 1:ncol(k2)){
            for(y2 in 1:nrow(k2)){
              k[y1+y2-1, x1+x2-1] = k[y1+y2-1, x1+x2-1]+(k1[y1,x1]*k2[y2,x2])
            }
          }
        }
      }
      return(k)
    }, kernel)
  }
  if(!(ncol(kernel)%%2 && nrow(kernel)%%2)){
    errorCondition("Kernels must be an n by m martix where both n and m are odd")
  }
  if(!is.null(threshold)){
    for(x in 1:ncol(kernel)){
      for(y in 1:nrow(kernel)){
        kernel[y,x] = +(threshold<=kernel[y,x])
      }
    }
  }
  return(kernel)
}


beyer_alg <- function(alg, data, kernel, beta = 0.2, z = 0.5, threshold=NULL){
  return(alg(beyer_fix_data(data), beyer_fix_kernel(kernel, threshold), beta, z))
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
  return(beyer_alg(Rconnect:::.beyer_stretch, data, kernel, beta, z, threshold))
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
beyer_wrap <- function(data, kernel, beta = 0.2, z = 0.5, threshold=NULL){
  return(beyer_alg(Rconnect:::.beyer_wrap, data, kernel, beta, z, threshold))
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
beyer_reflect <- function(data, kernel, beta = 0.2, z = 0.5, threshold=NULL){
  return(beyer_alg(Rconnect:::.beyer_refect, data, kernel, beta, z, threshold))
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
beyer_zero <- function(data, kernel, beta = 0.2, z = 0.5, threshold=NULL){
  return(beyer_alg(Rconnect:::.beyer_zero, data, kernel, beta, z, threshold))
}

# the output is shrunk down by enough that it never reaches outside the data matrix in the first place
#' @export
beyer_shrink <- function(data, kernel, beta = 0.2, z = 0.5, threshold=NULL){
  return(beyer_alg(Rconnect:::.beyer_shrink, data, kernel, beta, z, threshold))
}
