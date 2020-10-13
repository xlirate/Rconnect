#include <Rcpp.h>
#include <execution>
#include <ranges>
#include <iostream>

#include "matrix_edge_cases.h"

using namespace Rcpp; 
using namespace std;

// [[Rcpp::plugins(cpp2a)]]

template < int (CLAMP)(int, int, int, int), bool NARROW=false, bool ZERO_OUT_OF_BOUNDS=false>
NumericMatrix _convolve(
    const NumericMatrix& data,
    const NumericMatrix& kernel){
  
  auto k_width  = kernel.ncol();
  auto k_height = kernel.nrow();
  
  auto d_width  = data.ncol();
  auto d_height = data.nrow();
  
  auto o_width  = d_width;
  auto o_height = d_height;
  
  if constexpr(NARROW){
    if(d_height-(k_height-1) < 1 || d_width-(k_height-1) < 1){
      NumericMatrix output(0, 0);
      return output;
    }else{
      o_width  = d_width-(k_width -1);
      o_height = d_height-(k_height-1);
    }
  }
  
  NumericMatrix output(o_height, o_width);
  
  const auto x_indexes = ranges::iota_view<int, int>(0, o_width);
  const auto y_indexes = ranges::iota_view<int, int>(0, o_height);
  
  vector<tuple<int, int, double>> k_points;
  
  for(int x = 0; x<k_width; x++){
    for(int y = 0; y<k_height; y++){
      if(kernel[y+x*k_height]){
        k_points.emplace_back(x,y, kernel[y+x*k_height]);
      }
    }
  }
  
  /*
   * Parallel stratagy
   *
   * Use threads per column, in any order
   *
   * within each column, do each data point sequentially and in any order
   *
   * within each data point, perform the multiply and add in any order and
   *  if possible all at the same time using wide vector instructions
   *
   */
  for_each(
    //execution::seq,
    //execution::par,
    execution::par_unseq,
    begin(x_indexes),
    end(x_indexes),
    [&] (auto o_x){
      for_each(
        execution::seq,
        begin(y_indexes),
        end(y_indexes),
        [&] (auto o_y){
          output[o_y+o_x*o_height] = transform_reduce(
            //execution::seq,
            execution::unseq,
            begin(k_points),
            end(k_points),
            0.0,
            plus<>{},
            [&](auto point){
              auto[k_x, k_y, k_v]= point;
              auto d_x = CLAMP(o_x, k_x, d_width,  k_width);
              auto d_y = CLAMP(o_y, k_y, d_height, k_height);
              if constexpr(ZERO_OUT_OF_BOUNDS){
                if(d_x < 0 || d_width <= d_x || d_y < 0 || d_height <= d_y){
                  return 0.0;
                }
              }
              return data[d_y+d_height*d_x]*k_v;
            });
        });
    });
  
  return output;
}




// [[Rcpp::export(.convolve_stretch)]]
NumericMatrix convolve_stretch(NumericMatrix data, std::vector<NumericMatrix> kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve<&_stretch>);
}



// [[Rcpp::export(.convolve_wrap)]]
NumericMatrix convolve_wrap(NumericMatrix data, const std::vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve<&_wrap>);
}


// [[Rcpp::export(.convolve_refect)]]
NumericMatrix convolve_refect(NumericMatrix data, const std::vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve<&_reflect>);
}


// [[Rcpp::export(.convolve_zero)]]
NumericMatrix convolve_zero(NumericMatrix data, const std::vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve<&_zero, false, true>);
}


// [[Rcpp::export(.convolve_shrink)]]
NumericMatrix convolve_shrink(NumericMatrix data, const std::vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve<&_shrink, true>);
}












