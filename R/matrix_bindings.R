#library(magick)
library(imagine)

convolve <- function(data, kernal, times = 1, normalize = FALSE){
  if(normalize){
    kernal <- normalize_kernal(kernal)
  }
  if(is.list(kernal)){
    for(rep in 1:times){
      for(k in kernal){
        data <- convolution2D(data, k)
      }
    }
  }else{
    data <- convolution2D(data, kernal, times=times)
  }
  return(data)
}


  
#to_layer <- function(d, ...){
#  args = list(...)
  # 
#}


#to_image <- function(data, nrow=1, ncol=1, byrow=FALSE){
#  d <- dim(data)
#  if(is.null(d)){
#    data <- matrix(data, nrow, ncol, byrow)
#    d <- dim(data)
#  }
#  
#  defs <- c(size=paste(as.character(d[2]),"x",as.character(d[1]),sep=""), Colorspace="Gray", type="GRAY", format="Gray", quantum="format=floating-point", depth="64")
#  
#  return(magick:::magick_image_readbin(writeBin(data, raw(), size=8, endian = "little") , defines = defs))
#}


#@export
#to_image <- function(m, ...){
#  img <- .to_image(m, ...)
#  print(img)
#}

#crab <- image_read("./images/crab.png")
#block <- image_read("./images/block.bmp")

#import secrets; print(f"c(\n{(','+chr(10)).join([','.join([f'{secrets.randbelow(256):4d}' for _ in range(10)]) for _ in range(10)])})");

random_vec <- c(
   30, 243, 195, 240, 208, 118, 172,  14, 230,  32,
  128, 162, 227, 164, 172, 230, 225, 114,  46, 231,
  210, 123,  59, 193, 108,  77, 226, 185, 119,  79,
   54, 116, 189, 210, 240, 184,  89, 185, 185, 243,
  170, 217, 134, 234,  85,  24, 122, 112, 187, 175,
  116, 146,  17, 184, 247,  93,  43,  49, 174,  42,
   64,   3, 250,  36, 169, 222, 139, 177, 150, 172,
  150,  17, 236,  54,  37, 138, 250, 161,  46, 174,
   21, 247, 109,  57, 225,   9,  32,  81,   9, 134,
   85, 251, 223, 161,   8,  41, 141, 159, 140, 119)/256

centered_vec <- c(
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0)

# random_float <- random_vec/256

random_2d <- matrix(random_vec, 10)

centered_2d <- matrix(centered_vec, 11)

#random_3d <- array(random_vec, c(10, 10, 3))




#plot(as.raster(random_block))


#block <- array(random_block, c(10, 10, 3))/255

# img <- magick::image_read(random_3d) %>% magick::image_scale("500x500")

#plot(img)

