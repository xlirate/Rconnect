#include <Rcpp.h> 
#include <execution> 
#include <ranges>
#include <iostream>
using namespace Rcpp;
using namespace std;
// https://www.researchgate.net/publication/321325064_2D_Image_Convolution_using_Three_Parallel_Programming_Models_on_the_Xeon_Phi

using namespace std;

// [[Rcpp::plugins(cpp2a)]]



void _convolve(
    const NumericMatrix& kernel, 
    const NumericMatrix& data, 
    NumericMatrix& output,
    const int_least32_t width,
    const int_least32_t height,
    const int_least32_t k_width,
    const int_least32_t k_height){
  
  const auto x_indexes = std::ranges::iota_view(0, width);
  const auto y_indexes = std::ranges::iota_view(0, height);
  
  const auto x_k_indexes = std::ranges::iota_view(0, k_width);
  const auto y_k_indexes = std::ranges::iota_view(0, k_height);
  
  const auto policy = std::execution::par_unseq;
  //const auto policy = std::execution::seq;
  
  for_each(
    policy,
    begin(x_indexes), 
    end(x_indexes), 
    [&] (auto x){
      for_each(
        policy,
        begin(y_indexes), 
        end(y_indexes),
        [&] (auto y){
          output[y+x*height] = transform_reduce(
            policy,
            begin(y_k_indexes), 
            end(y_k_indexes),
            0.0,
            std::plus<>{},
            [&](auto ky){
              return transform_reduce(
                policy,
                begin(x_k_indexes), 
                end(x_k_indexes),
                0.0,
                std::plus<>{},
                [&](auto kx){
                  return kernel[ky+kx*k_height]*data[clamp(y+ky-k_height/2, 0, height-1)+height*clamp(x+kx-k_width/2, 0, width-1)];
                });
            });
        });
    });
}

// [[Rcpp::export]]
NumericMatrix convolve(NumericMatrix data, NumericMatrix kernel) {
  auto d_width = data.ncol();
  auto d_height = data.nrow();
  auto k_width = kernel.ncol();
  auto k_height = kernel.nrow();
  
  NumericMatrix output(d_height, d_width);
  
  _convolve(kernel, data, output, d_width, d_height, k_width, k_height);
  
  return output;
  
}