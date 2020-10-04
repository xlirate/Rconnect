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
    const size_t width,
    const size_t height,
    const size_t k_width,
    const size_t k_height,
    const vector<pair<size_t, size_t>>& k_points){
  
  const auto x_indexes = std::ranges::iota_view<size_t, size_t>(0, width);
  const auto y_indexes = std::ranges::iota_view<size_t, size_t>(0, height);
  
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
    [&] (auto x){
      for_each(
        std::execution::seq,
        begin(y_indexes), 
        end(y_indexes),
        [&] (auto y){
          output[y+x*height] = transform_reduce(
            std::execution::unseq,
            begin(k_points), 
            end(k_points),
            0.0,
            std::plus<>{},
            [&](auto point){
                auto[k_x, k_y]= point;
                return kernel.end()[-1-(k_y+k_x*k_height)]*data[clamp(y+k_y-k_height/2, (size_t)0, height-1)+height*clamp(x+k_x-k_width/2, (size_t)0, width-1)];
              });
        });
    });
}

// [[Rcpp::export]]
NumericMatrix convolve_cpp(NumericMatrix data, vector<NumericMatrix> kernel, unsigned long count = 1) {
  auto d_width = data.ncol();
  auto d_height = data.nrow();
  
  NumericMatrix _output_1(d_height, d_width);
  NumericMatrix _output_2(d_height, d_width);
  
  NumericMatrix* output_1 = &_output_1;
  NumericMatrix* output_2 = &_output_2;
  
  vector<vector<pair<size_t, size_t>>> k_points;
  k_points.reserve(kernel.size());
  
  const auto k_indexes = std::ranges::iota_view<size_t, size_t>(0, kernel.size());
  
  bool first_rep = true;
  
  for(unsigned long c = 0; c<count; c++){
    for(auto k_index : k_indexes){  
      auto k_width = kernel[k_index].ncol();
      auto k_height = kernel[k_index].nrow();
      
      if(k_index >= k_points.size()){
        k_points.emplace_back();
        k_points.back().reserve(k_width*k_height);
        
        for(size_t x = 0; x<k_width; x++){
          for(size_t y = 0; y<k_height; y++){
            k_points.back().emplace_back(x,y);
          }
        }
      }
      
      if(first_rep){
        first_rep = false;
        _convolve(kernel[k_index],  data,     *output_2, d_width, d_height, k_width, k_height, k_points[k_index]);
      }else{
        _convolve(kernel[k_index], *output_1, *output_2, d_width, d_height, k_width, k_height, k_points[k_index]);
      }
      auto t = output_2;
      output_2 = output_1;
      output_1 = t;
    }
  }
  
  if(first_rep){
    return data;
  }else{
    return *output_1;
  }
}


