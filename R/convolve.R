
convolve <- function(data, kernal, reps = 1L, normalize = TRUE){
  if(reps < 0){
    errorCondition("You cannot have negitive repititions")
  }
  if(length(reps) != 1){
    errorCondition("Repititons must have length 1")
  }
  if(!is.matrix(data)){
    errorCondition("Data must be a martix")
  }
  if(!is.list(kernal)){
    kernal <- list(kernal)
  }
  if(normalize){
    kernal <- normalize_kernal(kernal)
  }
  return(convolve_cpp(data, kernal, reps))
}