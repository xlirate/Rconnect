#include <Rcpp.h> 
#include <execution> 
#include <ranges>
#include <iostream>

using namespace Rcpp;
using namespace std;
// https://www.researchgate.net/publication/321325064_2D_Image_Convolution_using_Three_Parallel_Programming_Models_on_the_Xeon_Phi

using namespace std;

// [[Rcpp::plugins(cpp2a)]]


NumericMatrix _convolve_stretch(
    const NumericMatrix& data,
    const NumericMatrix& kernel){
  
  auto d_width  = data.ncol();
  auto d_height = data.nrow();
  
  auto k_width  = kernel.ncol();
  auto k_height = kernel.nrow();
  
  NumericMatrix output(d_height, d_width);
  
  const auto x_indexes = std::ranges::iota_view<int, int>(0, d_width);
  const auto y_indexes = std::ranges::iota_view<int, int>(0, d_height);
  
  vector<pair<int, int>> k_points;
  k_points.reserve(k_width*k_height);
  
  for(size_t x = 0; x<k_width; x++){
    for(size_t y = 0; y<k_height; y++){
      k_points.emplace_back(x,y);
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
    std::execution::par,
    begin(x_indexes), 
    end(x_indexes), 
    [&] (auto o_x){
      for_each(
        std::execution::seq,
        begin(y_indexes),
        end(y_indexes),
        [&] (auto o_y){
          output[o_y+o_x*d_height] = transform_reduce(
            std::execution::unseq,
            begin(k_points), 
            end(k_points),
            0.0,
            std::plus<>{},
            [&](auto point){
              auto[k_x, k_y]= point;
              auto d_x = clamp(o_x + k_x - k_width /2, 0, d_width -1);
              auto d_y = clamp(o_y + k_y - k_height/2, 0, k_height-1);
              return data[d_y+d_height*d_x]*kernel[k_y+k_height*k_x];
            });
        });
    });
  
  return output;
}


// [[Rcpp::export]]
NumericMatrix convolve_stretch_cpp(NumericMatrix data, const vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve_stretch);
}


NumericMatrix _convolve_wrap(
    const NumericMatrix& data,
    const NumericMatrix& kernel){
  
  auto d_width  = data.ncol();
  auto d_height = data.nrow();
  
  auto k_width  = kernel.ncol();
  auto k_height = kernel.nrow();
  
  NumericMatrix output(d_height, d_width);
  
  const auto x_indexes = std::ranges::iota_view<int, int>(0, d_width);
  const auto y_indexes = std::ranges::iota_view<int, int>(0, d_height);
  
  vector<pair<int, int>> k_points;
  k_points.reserve(k_width*k_height);
  
  for(size_t x = 0; x<k_width; x++){
    for(size_t y = 0; y<k_height; y++){
      k_points.emplace_back(x,y);
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
    std::execution::par,
    begin(x_indexes), 
    end(x_indexes), 
    [&] (auto o_x){
      for_each(
        std::execution::seq,
        begin(y_indexes),
        end(y_indexes),
        [&] (auto o_y){
          output[o_y+o_x*d_height] = transform_reduce(
            std::execution::unseq,
            begin(k_points), 
            end(k_points),
            0.0,
            std::plus<>{},
            [&](auto point){
              auto[k_x, k_y]= point;
              auto d_x = (o_x + k_x - k_width /2)%(d_width );
              auto d_y = (o_y + k_y - k_height/2)%(d_height);
              return data[d_y+d_height*d_x]*kernel[k_y+k_height*k_x];
            });
        });
    });
  
  return output;
}

// [[Rcpp::export]]
NumericMatrix convolve_wrap_cpp(NumericMatrix data, const vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve_wrap);
}

NumericMatrix _convolve_reflect(
    const NumericMatrix& data,
    const NumericMatrix& kernel){
  
  auto d_width  = data.ncol();
  auto d_height = data.nrow();
  
  auto k_width  = kernel.ncol();
  auto k_height = kernel.nrow();
  
  NumericMatrix output(d_height, d_width);
  
  const auto x_indexes = std::ranges::iota_view<int, int>(0, d_width);
  const auto y_indexes = std::ranges::iota_view<int, int>(0, d_height);
  
  vector<pair<int, int>> k_points;
  k_points.reserve(k_width*k_height);
  
  for(size_t x = 0; x<k_width; x++){
    for(size_t y = 0; y<k_height; y++){
      k_points.emplace_back(x,y);
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
    std::execution::par,
    begin(x_indexes), 
    end(x_indexes), 
    [&] (auto o_x){
      for_each(
        std::execution::seq,
        begin(y_indexes),
        end(y_indexes),
        [&] (auto o_y){
          output[o_y+o_x*d_height] = transform_reduce(
            std::execution::unseq,
            begin(k_points), 
            end(k_points),
            0.0,
            std::plus<>{},
            [&](auto point){
              auto[k_x, k_y]= point;
              
              auto d_x = abs(o_x + k_x - k_width /2);
              auto d_y = abs(o_y + k_y - k_height/2);
              
              d_x = ((d_x/d_width )%2)?((d_width -1)-((d_x)%d_width )):(d_x%d_width );
              d_x = ((d_y/d_height)%2)?((d_height-1)-((d_y)%d_height)):(d_y%d_height);
                
              return data[d_y+d_height*d_x]*kernel[k_y+k_height*k_x];
            });
        });
    });
  return output;
}


// [[Rcpp::export]]
NumericMatrix convolve_refect_cpp(NumericMatrix data, const vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve_reflect);
}

NumericMatrix _convolve_zero(
    const NumericMatrix& data,
    const NumericMatrix& kernel){
  
  auto d_width  = data.ncol();
  auto d_height = data.nrow();
  
  auto k_width  = kernel.ncol();
  auto k_height = kernel.nrow();
  
  NumericMatrix output(d_height, d_width);
  
  const auto x_indexes = std::ranges::iota_view<int, int>(0, d_width);
  const auto y_indexes = std::ranges::iota_view<int, int>(0, d_height);
  
  vector<pair<int, int>> k_points;
  k_points.reserve(k_width*k_height);
  
  for(size_t x = 0; x<k_width; x++){
    for(size_t y = 0; y<k_height; y++){
      k_points.emplace_back(x,y);
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
    std::execution::par,
    begin(x_indexes), 
    end(x_indexes), 
    [&] (auto o_x){
      for_each(
        std::execution::seq,
        begin(y_indexes),
        end(y_indexes),
        [&] (auto o_y){
          output[o_y+o_x*d_height] = transform_reduce(
            std::execution::unseq,
            begin(k_points), 
            end(k_points),
            0.0,
            std::plus<>{},
            [&](auto point){
              auto[k_x, k_y]= point;
              auto d_x = o_x + k_x - k_width /2;
              auto d_y = o_y + k_y - k_height/2;
              return (0 <= d_x && d_x < d_width && 0 <= d_y && d_y < d_height) ? (data[d_y+d_height*d_x]*kernel[k_y+k_height*k_x]) : (0);
            });
        });
    });
  
  return output;
}


// [[Rcpp::export]]
NumericMatrix convolve_zero_cpp(NumericMatrix data, const vector<NumericMatrix>& kernel, const unsigned long count = 1) {
  return accumulate(begin(kernel), end(kernel), data, _convolve_zero);
}



NumericMatrix _convolve_shrink(
    const NumericMatrix& data,
    const NumericMatrix& kernel){
  
  auto d_width  = data.ncol();
  auto d_height = data.nrow();
  
  auto k_width  = kernel.ncol();
  auto k_height = kernel.nrow();
  
  if(d_height-(k_height-1) < 1 || d_width-(k_height-1) < 1){
    NumericMatrix output(0, 0);
    return output;
  }
  
  NumericMatrix output(d_height-(k_height-1), d_width-(k_height-1));

  auto o_width  = output.ncol();
  auto o_height = output.nrow();
  
  const auto x_indexes = std::ranges::iota_view<int, int>(0, o_width);
  const auto y_indexes = std::ranges::iota_view<int, int>(0, o_height);
  
  vector<pair<int, int>> k_points;
  k_points.reserve(k_width*k_height);
  
  for(size_t x = 0; x<k_width; x++){
    for(size_t y = 0; y<k_height; y++){
      k_points.emplace_back(x,y);
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
    std::execution::par,
    begin(x_indexes), 
    end(x_indexes), 
    [&] (auto o_x){
      for_each(
        std::execution::seq,
        begin(y_indexes),
        end(y_indexes),
        [&] (auto o_y){
          output[o_y+o_x*o_height] = transform_reduce(
            std::execution::unseq,
            begin(k_points), 
            end(k_points),
            0.0,
            std::plus<>{},
            [&](auto point){
              auto[k_x, k_y]= point;
              auto d_x = o_x + k_x;
              auto d_y = o_y + k_y;
              return data[d_y+d_height*d_x]*kernel[k_y+k_height*k_x];
            });
        });
    });
  
  return output;
}


// [[Rcpp::export]]
NumericMatrix convolve_shrink_cpp(NumericMatrix data, const vector<NumericMatrix>& kernel) {
  return accumulate(begin(kernel), end(kernel), data, _convolve_shrink);
}



