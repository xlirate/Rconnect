#include <Rcpp.h>
#include <execution>
#include <ranges>
#include <iostream>
#include <ratio>
#include <limits>

#include "matrix_edge_cases.h"

using namespace Rcpp; 
using namespace std;

// [[Rcpp::plugins(cpp2a)]]

template < int (CLAMP)(int, int, int, int), bool NARROW=false, bool DEFAULT_OUT_OF_BOUNDS=false, intmax_t DEFAULT_NUM=0, intmax_t DEFAULT_DEN=1>
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
      /* 
       * The mask is large enough that there is no output, so we just return an empty matrix 
       * 
       * This could be made into an error condition instead, or filtered out up in R before getting here, 
       * but it is a cheap check and missing it might lead to a segfault
       *
       */
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
              if constexpr(DEFAULT_OUT_OF_BOUNDS){
                if(d_x < 0 || d_width <= d_x || d_y < 0 || d_height <= d_y){
                  if constexpr(DEFAULT_DEN != 0){
                    return double(DEFAULT_NUM)/DEFAULT_DEN;
                  }else if constexpr(DEFAULT_NUM > 0){
                    return std::numeric_limits<double>::infinity();
                  }else if constexpr(DEFAULT_NUM < 0){
                    return -std::numeric_limits<double>::infinity();
                  }else{
                    return std::numeric_limits<double>::quiet_NaN();
                  }
                }
              }
              return data[d_y+d_height*d_x]*k_v;
            });
        });
    });
  
  return output;
}


// [[Rcpp::export(.convolve_stretch)]]
NumericMatrix convolve_stretch(const NumericMatrix& data, const NumericMatrix& kernel) {
  return _convolve<&_stretch>(data, kernel);
}

// [[Rcpp::export(.convolve_wrap)]]
NumericMatrix convolve_wrap(const NumericMatrix& data, const NumericMatrix& kernel) {
  return _convolve<&_wrap>(data, kernel);
}

// [[Rcpp::export(.convolve_reflect)]]
NumericMatrix convolve_reflect(const NumericMatrix& data, const NumericMatrix& kernel) {
  return _convolve<&_reflect>(data, kernel);
}

// [[Rcpp::export(.convolve_zero)]]
NumericMatrix convolve_zero(const NumericMatrix& data, const NumericMatrix& kernel) {
  return _convolve<&_default, false, true, 0, 1>(data, kernel);
}

// [[Rcpp::export(.convolve_nan)]]
NumericMatrix convolve_nan(const NumericMatrix& data, const NumericMatrix& kernel) {
  return _convolve<&_default, false, true, 0, 0>(data, kernel);
}

// [[Rcpp::export(.convolve_shrink)]]
NumericMatrix convolve_shrink(const NumericMatrix& data, const NumericMatrix& kernel) {
  return _convolve<&_shrink, true>(data, kernel);
}



template <NumericMatrix (CONV)(const NumericMatrix&, const NumericMatrix&)>
NumericMatrix _powered_convolve(
    const NumericMatrix& data,
    const NumericMatrix& kernel,
    const double power=1.0){
  if(power == 1.0){
    return CONV(data, kernel);
  }else{
    NumericMatrix powered_data(data.nrow(), data.ncol());
    transform(
      execution::par_unseq,
      begin(data),
      end(data),
      begin(powered_data),
      [power](auto d){
      return std::pow(d, power);
    });
    return CONV(powered_data, kernel);
  }
}

// [[Rcpp::export(.powered_convolve_stretch)]]
NumericMatrix powered_convolve_stretch(const NumericMatrix& data, const NumericMatrix& kernel, const double power=1.0) {
  return _powered_convolve<&convolve_stretch>(data, kernel, power);
}

// [[Rcpp::export(.powered_convolve_wrap)]]
NumericMatrix powered_convolve_wrap(const NumericMatrix& data, const NumericMatrix& kernel, const double power=1.0) {
  return _powered_convolve<&convolve_wrap>(data, kernel, power);
}

// [[Rcpp::export(.powered_convolve_refect)]]
NumericMatrix powered_convolve_refect(const NumericMatrix& data, const NumericMatrix& kernel, const double power=1.0) {
  return _powered_convolve<&convolve_reflect>(data, kernel, power);
}

// [[Rcpp::export(.powered_convolve_zero)]]
NumericMatrix powered_convolve_zero(const NumericMatrix& data, const NumericMatrix& kernel, const double power=1.0) {
  return _powered_convolve<&convolve_zero>(data, kernel, power);
}

// [[Rcpp::export(.powered_convolve_nan)]]
NumericMatrix powered_convolve_nan(const NumericMatrix& data, const NumericMatrix& kernel, const double power=1.0) {
  return _powered_convolve<&convolve_nan>(data, kernel, power);
}

// [[Rcpp::export(.powered_convolve_shrink)]]
NumericMatrix powered_convolve_shrink(const NumericMatrix& data, const NumericMatrix& kernel, const double power=1.0) {
  return _powered_convolve<&convolve_shrink>(data, kernel, power);
}










