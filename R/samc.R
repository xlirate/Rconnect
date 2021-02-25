


samc_cache <- function(kernel, permiability, lethality=NULL){
  #kernel should be n*m where n and m are odd
  if(is.null(lethality)){
    lethality <- matrix(0, nrow(permiability), ncol(permiability))
  }
  
  return(.cache_samc(kernel, permiability, lethality))
}


samc_one_step <- function(cache, population, dead=NULL){
  sizes <- Rconnect:::.samc_cache_sizes(cache)
  # {ca->nrow, ca->ncol, ca->left_extra_cols, ca->right_extra_cols};
  s_pop <- NA
  s_dead <- NA
  
  if(nrow(population) != sizes[1]){
    stop("Population has the wrong height")
  }else if(ncol(population) == sizes[2]){
    s_pop <- matrix(0, nrow=sizes[1], ncol=sizes[2]+sizes[3]+sizes[4])
    
    s_pop[1:sizes[1], (1+sizes[3]):(sizes[2]+sizes[3])] <- population
  }else if(ncol(population) == sizes[2]+sizes[3]+sizes[4]){
    s_pop <- population
  }else{
    stop("Population has the wrong width")
  }
  
  if(is.null(dead)){
    s_dead <- matrix(0, nrow=sizes[1], ncol=sizes[2]+sizes[3]+sizes[4])
  }else if(nrow(dead) != sizes[1]){
    stop("Dead has the wrong height")
  }else if(ncol(dead) == sizes[2]){
    s_dead <- matrix(0, nrow=sizes[1], ncol=sizes[2]+sizes[3]+sizes[4])
    s_dead[1:sizes[1], (1+sizes[3]):(sizes[2]+sizes[3])] <- dead
  }else if(ncol(dead) == sizes[2]+sizes[3]+sizes[4]){
    s_dead <- dead
  }else{
    stop("Dead has the wrong width")
  }
  
  out <- Rconnect:::.samc_one_step(cache, s_pop, s_dead)
  #return(out)
  if(ncol(population) == sizes[2]){
    return(out[1:sizes[1], (1+sizes[3]):(sizes[2]+sizes[3])])
  }else{
    return(out)
  }
}