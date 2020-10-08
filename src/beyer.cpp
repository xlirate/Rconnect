#include <Rcpp.h>
#include <execution>
#include <ranges>
#include <iostream>

#include "matrix_edge_cases.hpp"

using namespace Rcpp; 
using namespace std;

// [[Rcpp::plugins(cpp2a)]]

template < int (CLAMP)(int, int, int, int), bool NARROW=false, bool ZERO_OUT_OF_BOUNDS=false>
NumericMatrix _beyer(
    const NumericMatrix& data,
    const NumericMatrix& mask,
    const double beta = 0.2,
    const double z = 0.5){
  
  auto m_width  = mask.ncol();
  auto m_height = mask.nrow();
  
  auto d_width  = data.ncol();
  auto d_height = data.nrow();
  
  auto o_width  = d_width;
  auto o_height = d_height;
  
  if constexpr(NARROW){
    if(d_height-(m_height-1) < 1 || d_width-(m_height-1) < 1){
      NumericMatrix output(0, 0);
      return output;
    }else{
      o_width  = d_width-(m_width -1);
      o_height = d_height-(m_height-1);
    }
  }
  
  NumericMatrix output(o_height, o_width);
  
  const auto x_indexes = ranges::iota_view<int, int>(0, o_width);
  const auto y_indexes = ranges::iota_view<int, int>(0, o_height);
  
  vector<tuple<int, int, double, double>> m_points;
  //m_points.reserve(m_width*m_height);
  double max_e = 0;
  
  for(int x = 0; x<m_width; x++){
    for(int y = 0; y<m_height; y++){
      if(mask[y+x*m_height]){
        auto delta_m_e = std::exp(-beta*sqrt((x-m_width /2)*(x-m_width /2)+
                                  (y-m_height/2)*(y-m_height/2)));
        m_points.emplace_back(
          x,
          y,
          mask[y+x*m_height],
          delta_m_e);
        max_e+= delta_m_e;
      }
    }
  }
  
  const auto m0 = mask[m_height/2 + (m_width/2)*m_height];
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
          const auto d0 = data[CLAMP(o_y, m_width/2, d_height, m_height) + d_height*CLAMP(o_x, m_width/2, d_width, m_width)];
          output[o_y+o_x*o_height] = transform_reduce(
            //execution::seq,
            execution::unseq,
            begin(m_points),
            end(m_points),
            0.0,
            plus<>{},
            [&](auto point){
              auto[m_x, m_y, m_v, m_e]= point;
              auto d_x = CLAMP(o_x, m_x, d_width,  m_width);
              auto d_y = CLAMP(o_y, m_y, d_height, m_height);
              if constexpr(ZERO_OUT_OF_BOUNDS){
                if(d_x < 0 || d_width <= d_x || d_y < 0 || d_height <= d_y){
                  return 0.0;
                }
              }
              return std::pow(data[d_y+d_height*d_x]*m_v* d0, z)*m_e;
            })/max_e;
        });
    });
  
  return output;
}


// [[Rcpp::export(.beyer_stretch)]]
NumericMatrix beyer_stretch(
    const NumericMatrix& data, 
    const NumericMatrix& mask, 
    const double beta = 0.2,
    const double z = 0.5) {
  return _beyer<&_stretch>(data, mask, beta, z);
}


// [[Rcpp::export(.beyer_wrap)]]
NumericMatrix beyer_wrap(
    const NumericMatrix& data, 
    const NumericMatrix& mask, 
    const double beta = 0.2,
    const double z = 0.5) {
  return _beyer<&_wrap>(data, mask, beta, z);
}


// [[Rcpp::export(.beyer_refect)]]
NumericMatrix beyer_refect(
    const NumericMatrix& data, 
    const NumericMatrix& mask, 
    const double beta = 0.2,
    const double z = 0.5) {
  return _beyer<&_reflect>(data, mask, beta, z);
}


// [[Rcpp::export(.beyer_zero)]]
NumericMatrix beyer_zero(
    const NumericMatrix& data, 
    const NumericMatrix& mask, 
    const double beta = 0.2,
    const double z = 0.5) {
  return _beyer<&_zero, false, true>(data, mask, beta, z);
}


// [[Rcpp::export(.beyer_shrink)]]
NumericMatrix beyer_shrink(
    const NumericMatrix& data, 
    const NumericMatrix& mask, 
    const double beta = 0.2,
    const double z = 0.5) {
  return _beyer<&_shrink, true>(data, mask, beta, z);
}





