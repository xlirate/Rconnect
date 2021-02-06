#' @export
normalize_kernel <- function(k){
  if(is.list(k)){
    lapply(k, normalize_kernel)
  }else{
    if(s <- sum(k)){
      k/s
    }else{
      k
    }
  }
}

#' @export
flatten_kernel <- function(kernel){
  if(is.list(kernel)){
    kernel <- Reduce(function(k1, k2){return(k1 %*% k2)}, kernel)
  }
  return(kernel)
}

.qkernel_to_kernel <- function(q){
  if(typeof(q) == "list"){
    lapply(q, .qkernel_to_kernel)
  }else{
    rows <- nrow(q)
    cols <- ncol(q)
    if(rows == 1 && cols == 1){
      # 1 x 1
      q
    }else if(rows == 1){
      # 1 x cols
      matrix(c(q[1, cols:2], q), 1)
    }else if(cols == 1){
      # rows x 1
      matrix(c(q[rows:2], q))
    }else if(rows == 2 && cols == 2){
      # 2 x 2
      matrix(
        sapply(
          c(4,3,4,2,1,2,4,3,4),
          function(n){q[n]}
        ),
        3)
    }else if(rows == 2){
      # 2 x cols
      rbind(
        c(rev(q[-1,-1]), q[-1,]),
        c(rev(q[ 1,-1]), q[ 1,]),
        c(rev(q[-1,-1]), q[-1,])
      )
    }else if(cols == 2){
      # rows x 2
      cbind(
        c(rev(q[-1,-1]),q[,-1]),
        c(rev(q[-1, 1]),q[, 1]),
        c(rev(q[-1,-1]),q[,-1])
      )
    }else{
      # 3x3 input or larger
      rbind(cbind(q[rows:2, cols:2], q[rows:2,]),
            cbind(q[      , cols:2], q))
    }
  }
}

.hard_uniform_circle_qkernel <- function(radius) {
  qm <- matrix(0:radius, radius+1, radius+1)
  +(qm^2 + t(qm)^2 <= radius^2)
}

#' @export
hard_uniform_circle_kernel <- function(radius) {
  .qkernel_to_kernel(.hard_uniform_circle_qkernel(radius))
}

.half_circle_integral <- function(r, a, b) {
  # This is $$\int_{a}^b \sqrt{r^2 - x^2} dx$$ simplified slightly
  # integrate from a to b for sqrt(r^2 - x^2) dx
  (b*sqrt(r^2-b^2)-a*sqrt(r^2-a^2)+r^2*(atan2(b, sqrt(r^2-b^2))-atan2(a, sqrt(r^2-a^2))))/2
}

.square_covered_portion <- function(r, x, y){
  x <- abs(x)-1
  y <- abs(y)-1
  if (x < y){
    t <- x
    x <- y
    y <- t
  }
  #now 0,0 is the origin, and 0 <= y <= x
  #only needs to be correct if r > 1.5
  if(x == 0){
    # At the origin
    if (r <= (1/2)){
      # The entire circle fits in the cell
      #      |
      #   O--O--O
      #   |  |  |
      # --O-(+)-$--
      #   |  |  |
      #   O--O--O
      #      |

      pi*r*r
    }else if(r <= sqrt(2)/2){
      # The circle reaches out by the right term on each of the 4 sides
      #     /-\
      #   O/-X-\$
      #   /  |  \
      #--(X--+--X)--
      #   \  |  /
      #   O\-X-/O
      #     \-/
      pi*r*r-4*(r*r*acos(0.5/r)-(0.5)*sqrt(r*r-0.5*0.5))
    }else{
      # The circle covers the entire origin cell
      #  /   |   \
      # / X--X--X \
      #   |  |  |
      # --X--+--X--
      #   |  |  |
      # \ X--X--X /
      #  \   |   /
      1
    }
  }else if (y == 0){
    # On the axis, but not at the origin
    if(r < (x-0.5)){
      # not touched by the circle
      #
      #\  O-----O
      # \ |     |
      #--)$-----O--
      # / |     |
      #/  O-----O
      #
      0
    }else if(r^2 <= (((x-0.5)^2)+(0.5^2))){
      # Not covering the first pair of corners
      #
      #  \$-----O
      #   \     |
      #---X)----O--
      #   /     |
      #  /O-----O
      #
      crossover <- sqrt(r^2-(x-0.5)^2)
      2*(.half_circle_integral(r, 0, crossover)-crossover*(x-0.5))
    }else if(r <= (x+0.5)){
      # Covering the first pair of corners, but not reaching the next cell over
      #    \
      #   X-\---O
      #   |  \  |
      #---X---)-$--
      #   |  /  |
      #   X-/---O
      #    /

      .half_circle_integral(r, -0.5, 0.5)-(x-0.5)
    }else if(r^2 < (((x+0.5)^2)+(0.5^2))){
      # Reaching in to the next cell over, but not covering this cell completely
      #       \
      #   X----\$
      #   |     \
      #---X-----X)-
      #   |     /
      #   X----/O
      #       /
      crossover <- sqrt(r^2-(x+0.5)^2)
      2*(.half_circle_integral(r, crossover, 0.5)-(0.5-crossover)*(x-0.5)+(crossover))
    }else{
      # The cell is covered
      #          \
      #   X-----X \
      #   |     |  \
      #---X-----X---)
      #   |     |  /
      #   X-----X /
      #          /
      1
    }

  }else if(r^2 < ((x-0.5)^2+(y-0.5)^2)){
    # not covered at all, or are on the axis
    #
    #   O    O
    #\
    # \
    #  \$    O
    #   \
    0
  }else if(r^2 <= (x-0.5)^2+(y+0.5)^2){
    # only first corner is covered
    # \
    #  \$    O
    #   \
    #    \
    #   X \  O
    #      \
    crossover <- sqrt(r^2-(x-0.5)^2)
    .half_circle_integral(r, y-0.5, crossover)-(crossover-(y-0.5))*(x-0.5)
  }else if(r^2 <= (x+0.5)^2+(y-0.5)^2){
    # First 2 corners are covered
    #   \
    #   X\   O
    #     \
    #      \
    #   x   \$
    #        \
    .half_circle_integral(r, y-0.5, y+0.5)-(x-0.5)
  }else if(r^2 <= (x+0.5)^2+(y+0.5)^2){
    # Three corners are covered
    #      \
    #   X   \$
    #        \
    #         \
    #   X    X \
    #           \
    crossover <- sqrt(r^2-(x+0.5)^2)
    .half_circle_integral(r, crossover, y+0.5)-((y+0.5)-crossover)*(x-0.5)+(crossover-(y-0.5))
  }else{
    #Completely covered
    #         \
    #   X    X \
    #           \
    #            \
    #   X    X    \
    #              \
    1
  }
}

.smooth_uniform_circle_qkernel <- function(r) {
  qmx = matrix(1:(r+1.5), r+1.5, r+1.5)
  qmy = matrix(1:(r+1.5), r+1.5, r+1.5, byrow=TRUE)
  matrix(mapply(.square_covered_portion, r, qmx, qmy), r+1.5, r+1.5)
}

#' @export
smooth_uniform_circle_kernel <- function(r){
  .qkernel_to_kernel(.smooth_uniform_circle_qkernel(r))
}

.distance_qkernel <- function(r){
  #r <- ceiling(r)
  x <- matrix(0:r, r+1, r+1)
  y <- matrix(0:r, r+1, r+1, byrow=TRUE)
  matrix(mapply(function(x,y){sqrt(x*x+y*y)}, x, y),r+1, r+1)
}

#' @export
distance_kernel <- function(r){
  .qkernel_to_kernel(.distance_qkernel(r))
}

.exponential_qkernel<-function(dbar,cellDim=1,negligible=10^-10,returnScale=F,dmax=NULL){
  #Exponential kernel from Hughes et al 2015 American Naturalist
  dbarCell <- dbar/cellDim
  if(is.null(dmax)){
    dmax <- -0.5*dbarCell*log(pi*dbarCell^2*negligible/2)
    if(dmax<0){
      stop("Set negligible so that pi*dbar^2*negligible/2 <=1")
    }
  }
  m <- (2/(pi*dbarCell^2))*exp(-2*.distance_qkernel(dmax)/dbarCell)
  m[m<negligible] <- 0
  if(returnScale){
    return((sum(m[-1,])*4)+m[1,1])#this would be sum(m) if we had the entire matrix, but we do not
  }else{
    return(m)
  }
}

#' @export
exponential_kernel <- function(dbar,cellDim=1,negligible=10^-10,returnScale=F,dmax=NULL){
  .qkernel_to_kernel(.exponential_qkernel(dbar, cellDim, negligible, returnScale, dmax))
}

.gaussian_qkernel <- function(sd, r0=0.05){
  k <- matrix(exp(-(0:sqrt(-2*sd*sd*log(r0)))/(2*sd*sd)))
  list(k, t(k))
}

#' @export
gaussian_kernel <- function(sd, r0=0.05){
  .qkernel_to_kernel(.gaussian_qkernel(sd, r0))
}

.binomial_qkernal <- function(index){
  if(!is.integer(index) && !(index == floor(index))){
    stop("index must be an integer")
  }else if(index < 0){
    stop("index must be >= 0")
  }else if(index%%2){
    stop("index must be even. If it were odd, the kernal would not have a center")
  }else if(index == 0){
    return(matrix(1, 1))
  }else{
    
    full_row <- (lapply(index, function(i) choose(i, 0:i))[[1]])
    k <- matrix(full_row[-(index/2):0])
    return(list(k, t(k)))
  }
}

#' @export
binomial_kernal <- function(index){
  .qkernel_to_kernel(.binomial_qkernal(index))
}

.beyer_qkernel <- function(beta=0.2, r0=0.05){
  #assert that 0 < beta 
  #assert that 0 < r0 < 1
  
  # total area under the exponential distribution is (2*pi)/(beta*beta)
  # area within radius r is (2*pi)*(exp(-beta*r)*(-beta*r-1)+1)/(beta*beta)
  #
  # (2*pi)*(exp(-beta*r)*(-beta*r-1)+1)/(beta^2) < (1-r0)*(2*pi)/(beta^2)
  # exp(-beta*r)*(-beta*r-1) < -r0
  # 
  
  #count up the radius linearly until it is large enough. 
  #We could do a normal root finding, but this is good enough
  r = 1
  while(exp(-beta*r)*(-beta*r-1) < -r0){
    r = r+1
  }
  
  return(exp(-beta*.distance_qkernel(r))*.hard_uniform_circle_qkernel(r))
}


#' @export
beyer_kernel <- function(beta=0.2, r0=0.05){
  .qkernel_to_kernel(.beyer_qkernel(beta, r0))
}



