#ifndef __SAMC_H__
#define __SAMC_H__

#include <Rcpp.h>
#include <execution>
#include <ranges>

// [[Rcpp::plugins(cpp2a)]]

namespace samc{

struct kernel_point_t{
    std::ptrdiff_t x_off;
    std::ptrdiff_t y_off;
    double num = 1;
};

struct cache{
  std::size_t ncol;
  std::size_t nrow;
  std::size_t kernel_size;
  std::size_t left_extra_cols;
  std::size_t right_extra_cols;
  std::vector<double> movement_rate; /* size should be ncol*nrow*kernel_size */
  std::vector<double> death_rate; /* size should be ncol*nrow */
  std::vector<std::ptrdiff_t> kernel; /* a list of offsets from a given point to all of the places that it needs to look in the list of points */

};


inline std::ostream &operator << (std::ostream &os, const cache &ca) {
  os << "{\"size\":{\"ncol\":" << ca.ncol
    << ", \"nrow\": " << ca.nrow
    << ", \"kernel_size\": " << ca.kernel_size
    << ", \"movement_rate_map_size\": " << ca.movement_rate.size()
    << ", \"right_extra_cols\": " << ca.right_extra_cols
    << ", \"left_extra_cols\": " << ca.left_extra_cols
    << "}, \"offsets\":[";
  bool first = true;
  for(auto off : ca.kernel){
    if(first){first = false;}else{os << ", ";}
    os << off;
  }
  os << "], \"death_rate\":[";

  first = true;
  for(auto death : ca.death_rate){
    if(first){first = false;}else{os << ", ";}
    os << death;
  }

  os << "], \"movement_rate_map\":[";
  first = true;
  for(size_t i = 0; i< ca.movement_rate.size()/ca.kernel_size; i++){
    if(first){first = false;}else{os << ", ";}
    os << "[";
    bool second = true;
    for(size_t j = 0; j<ca.kernel_size; j++){
      if(second){second = false;}else{os << ", ";}
      os << ca.movement_rate[j+i*ca.kernel_size];
    }
    os << "]";
  }
  os << "]}\n";

  return os;
}

}/* namespace samc */

#endif /* __SAMC_H__ */
