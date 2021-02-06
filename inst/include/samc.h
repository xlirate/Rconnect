#ifndef __SAMC_H__
#define __SAMC_H__

#include <Rcpp.h>
#include <execution>
#include <ranges>

// [[Rcpp::plugins(cpp2a)]]

namespace samc{

struct kernal_point_t{
    std::ptrdiff_t x_off; 
    std::ptrdiff_t y_off;
    long num = 1;
};

template<std::size_t KERNAL_SIZE>
struct cache{
  std::size_t ncol;
  std::size_t nrow;
  std::size_t left_extra_cols;
  std::size_t right_extra_cols;
  std::vector<std::array<double, KERNAL_SIZE>> movement_rate; /* size should be ncol*nrow */
  std::vector<double> death_rate; /* size should be ncol*nrow */
  std::array<std::ptrdiff_t, KERNAL_SIZE> kernal; /* a list of offsets from a given point to all of the places that it needs to look in the list of points */
  
};


template<std::size_t KERNAL_SIZE>
std::ostream &operator << (std::ostream &os, const cache<KERNAL_SIZE> &ca) {
  os << "{ \"KERNAL_SIZE\":" << KERNAL_SIZE << ", \"mr_size\":" << ca.movement_rate.size()
    << ", \"size\":{\"ncol\":" << ca.ncol 
    << ", \"nrow\": " << ca.nrow 
    << ",\"left_extra_cols\": " << ca.left_extra_cols
    << ",\"right_extra_cols\": " << ca.right_extra_cols 
    << "}, \"offsets\":[";
  bool first = true;
  for(auto off : ca.kernal){
    if(first){first = false;}else{os << ", ";}
    os << off;
  }
  os << "], \"death_rate\":[";
  
  first = true;
  for(auto death : ca.death_rate){
    if(first){first = false;}else{os << ", ";}
    os << death;
  }
  
  os << "], \"movement_rate\":[";
  first = true;
  for(size_t i = 0; i< ca.movement_rate.size(); i++){
    if(first){first = false;}else{os << ", ";}
    os << "[";
    bool second = true;
    for(size_t j = 0; j<KERNAL_SIZE; j++){
      if(second){second = false;}else{os << ", ";}
      os << ca.movement_rate[i][j];
    }
    os << "]";
  }
  os << "]}\n";
  
  return os;
}

}/* namespace samc */

#endif /* __SAMC_H__ */
