#include <samc.h>
#include <Rcpp.h>
#include <execution>
#include <ranges>
#include <cstddef>
#include <vector>
#include <iostream>

// [[Rcpp::plugins(cpp2a)]]

namespace samc{

//constexpr std::size_t index_of(std::size_t x, std::ptrdiff_t x_off, std::size_t y, std::ptrdiff_t y_off, std::size_t height){
//  return y+y_off + (x+x_off)*height;
//}

//template<kernal_point_t ...KERNAL>
void construct_cache(cache& ca, const std::vector<kernal_point_t>& kernal, const Rcpp::NumericMatrix& permiability, const Rcpp::NumericMatrix& death_rate){

  //const size_t KERNAL_SIZE = sizeof...(KERNAL);
  //const std::array<kernal_point_t, KERNAL_SIZE> KERNAL_ARR{KERNAL...};
  ca.kernal_size = kernal.size();
  ca.nrow = permiability.nrow();
  ca.ncol = permiability.ncol();

  if(!(ca.ncol * ca.nrow)){
    throw "There are no cells, nrow*ncol == 0";
  }

  if(death_rate.ncol() != ca.ncol || death_rate.nrow() != ca.nrow){
    throw "death_rate's size does not match permiability's size";
  }


  ca.movement_rate.clear();
  ca.movement_rate.resize(ca.kernal_size*ca.nrow*ca.ncol, {0});
  //std::cout << ca.movement_rate.size();
  ca.death_rate.assign(death_rate.begin(),death_rate.end());


  std::ptrdiff_t max_offset = 0;
  std::ptrdiff_t min_offset = 0;
  for(const auto& k_point : kernal){
    const std::ptrdiff_t offset = k_point.y_off + k_point.x_off*ca.nrow;
    ca.kernal.push_back(-offset);
    max_offset = std::max(max_offset, offset);
    min_offset = std::min(min_offset, offset);
  }

  ca.left_extra_cols  = std::max((ca.nrow-1-min_offset)/ca.nrow, {0});
  ca.right_extra_cols = std::max((ca.nrow-1+max_offset)/ca.nrow, {0});

  const auto x_indexes = std::ranges::iota_view<size_t, size_t>(0, ca.ncol);
  const auto y_indexes = std::ranges::iota_view<size_t, size_t>(0, ca.nrow);
  const auto k_indexes = std::ranges::iota_view<size_t, size_t>(0, ca.kernal_size);


  std::for_each(
    //std::execution::par_unseq,
    std::execution::seq,
    std::begin(x_indexes),
    std::end(x_indexes),
    [&] (auto x){
      std::for_each(
        //std::execution::par_unseq,
        std::execution::seq,
        std::begin(y_indexes),
        std::end(y_indexes),
        [&] (auto y){
          //std::cout << x << ", " << y << ", " << ca.death_rate[y+x*ca.nrow] << "\n";
          //for each x,y

          double weighted_sum =
            std::transform_reduce(
            //std::execution::par_unseq,
            std::execution::seq,
            std::begin(kernal),
            std::end(kernal),
            0.0,
            std::plus<>{},
            [&] (auto& k_point){
              //const auto& k_point = kernal[k];
              const size_t k_x = x+k_point.x_off;
              const size_t k_y = y+k_point.y_off;
              //std::cout << "(" << x << ", " << y << ", " << k << ", ";
              if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
                //std::cout << permiability[k_y + k_x*ca.nrow] * k_point.num <<")\n";
                return permiability[k_y + k_x*ca.nrow] * k_point.num;
              }else{
                //std::cout << "0)\n";
                return 0.0;
              }
            });

          const double scalar = (weighted_sum)?((1.0-death_rate[y+x*ca.nrow])/weighted_sum):(0);
          //std::cout << "(" << x << "," << y << "," << scalar << ")\n";

          std::for_each(
            //std::execution::par_unseq,
            std::execution::seq,
            std::begin(k_indexes),
            std::end(k_indexes),
            [&] (auto k){
              const auto& k_point = kernal[k];
              const size_t k_x = x+k_point.x_off;
              const size_t k_y = y+k_point.y_off;


              if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
                ca.movement_rate[k+(k_y+k_x*ca.nrow)*ca.kernal_size] = scalar * permiability[k_y + k_x*ca.nrow];
              }
              //if no, leave it 0, 0 is the default value anyway
            });

        });
    });
}

} /* namespace samc  */

// [[Rcpp::export(.cache_samc_4)]]
Rcpp::XPtr<samc::cache> cache_samc_4(const Rcpp::NumericMatrix& permiability, const Rcpp::NumericMatrix& death_rate) {
  samc::cache* ca = new samc::cache;
  std::vector<samc::kernal_point_t> kernal{
                                   samc::kernal_point_t{-1, 0},
      samc::kernal_point_t{ 0,-1}, samc::kernal_point_t{ 0, 0}, samc::kernal_point_t{ 0, 1},
                                   samc::kernal_point_t{ 1, 0}};
  samc::construct_cache(*ca, kernal, permiability, death_rate);
  Rcpp::XPtr<samc::cache> xp(ca, true);
  return xp;
}

inline void samc_step(
    const samc::cache& ca,
    const Rcpp::NumericMatrix& pop_in,
    const Rcpp::NumericMatrix& dead_in,
    Rcpp::NumericMatrix& pop_out,
    Rcpp::NumericMatrix& dead_out){

  const auto indexes     = std::ranges::iota_view<size_t, size_t>(0, ca.death_rate.size());
  const auto connections = std::ranges::iota_view<size_t, size_t>(0, ca.kernal_size);

  const size_t offset = ca.nrow * ca.left_extra_cols;

  //std::cout << ca.death_rate.size() << " " << ca.kernal_size << " " << offset <<"\n";

  std::for_each(
    std::execution::par_unseq,
    //std::execution::seq,
    std::begin(indexes),
    std::end(indexes),
    [&] (auto i){
      dead_out[i+offset] = dead_in[i+offset]+ca.death_rate[i]*pop_in[i+offset];
      pop_out[i+offset] = transform_reduce(
        //std::execution::seq,
        std::execution::unseq,
        std::begin(connections),
        std::end(connections),
        0.0,
        std::plus<>{},
        [&](auto con){
          //std::cout << i << "\t" << con << "\t" << i*ca.kernal_size+con << "\t" << i+offset+ca.kernal[con]  <<"\n";
          //return ca.movement_rate[i][con]*pop_in[i+offset+ca.kernal[con]];
          return ca.movement_rate[i*ca.kernal_size+con]*pop_in[i+offset+ca.kernal[con]];
          //return 1.0;
        });
    });
}

// [[Rcpp::export(.samc_cache_sizes)]]
std::vector<size_t> samc_cache_sizes(const Rcpp::XPtr<samc::cache>& ca){
  return {ca->nrow, ca->ncol, ca->left_extra_cols, ca->right_extra_cols};
}

// [[Rcpp::export(.samc_print_cache)]]
void samc_print_cache(
    const Rcpp::XPtr<samc::cache>& ca){
  std::cout << *ca;
}


//Rcpp::SubMatrix<REALSXP> samc_one_step(
// [[Rcpp::export(.samc_one_step)]]
Rcpp::NumericMatrix samc_one_step(
    const Rcpp::XPtr<samc::cache>& ca,
    const Rcpp::NumericMatrix& pop_in,
    const Rcpp::NumericMatrix& dead_in){

  //return pop_in;

  Rcpp::NumericMatrix pop_out( ca->nrow, (ca->ncol + ca->left_extra_cols + ca->right_extra_cols));
  Rcpp::NumericMatrix dead_out(ca->nrow, (ca->ncol + ca->left_extra_cols + ca->right_extra_cols));

  samc_step(*ca, pop_in, dead_in, pop_out, dead_out);

  return pop_out;
}

/*

// [[Rcpp::export(.cache_samc_8)]]
Rcpp::XPtr<samc::cache<9>> cache_samc_8(const Rcpp::NumericMatrix& permiability, const Rcpp::NumericMatrix& death_rate) {
  samc::cache<9>* ca = new samc::cache<9>;
  samc::construct_cache<
    samc::kernal_point_t{-1,-1}, samc::kernal_point_t{-1, 0}, samc::kernal_point_t{-1, 1},
    samc::kernal_point_t{ 0,-1}, samc::kernal_point_t{ 0, 0}, samc::kernal_point_t{ 0, 1},
    samc::kernal_point_t{ 1,-1}, samc::kernal_point_t{ 1, 0}, samc::kernal_point_t{ 1, 1}
    >(*ca, permiability, death_rate);
  Rcpp::XPtr<samc::cache<9>> xp(ca, true);
  return xp;
}

// [[Rcpp::export(.cache_samc_binomial_2)]]
Rcpp::XPtr<samc::cache<9>> cache_samc_binomial_2(const Rcpp::NumericMatrix& permiability, const Rcpp::NumericMatrix& death_rate) {
  samc::cache<9>* ca = new samc::cache<9>;
  samc::construct_cache<
    samc::kernal_point_t{-1,-1, 1}, samc::kernal_point_t{-1, 0, 2}, samc::kernal_point_t{-1, 1, 1},
    samc::kernal_point_t{ 0,-1, 2}, samc::kernal_point_t{ 0, 0, 4}, samc::kernal_point_t{ 0, 1, 2},
    samc::kernal_point_t{ 1,-1, 1}, samc::kernal_point_t{ 1, 0, 2}, samc::kernal_point_t{ 1, 1, 1}
    >(*ca, permiability, death_rate);
  Rcpp::XPtr<samc::cache<9>> xp(ca, true);
  return xp;
}

// [[Rcpp::export(.cache_samc_binomial_4)]]
Rcpp::XPtr<samc::cache<25>> cache_samc_binomial_4(const Rcpp::NumericMatrix& permiability, const Rcpp::NumericMatrix& death_rate) {
  samc::cache<25>* ca = new samc::cache<25>;
  samc::construct_cache<
    samc::kernal_point_t{-2,-2, 1}, samc::kernal_point_t{-2,-1, 4}, samc::kernal_point_t{-2, 0, 6}, samc::kernal_point_t{-2, 1, 4}, samc::kernal_point_t{-2, 2, 1},
    samc::kernal_point_t{-1,-2, 4}, samc::kernal_point_t{-1,-1,16}, samc::kernal_point_t{-1, 0,24}, samc::kernal_point_t{-1, 1,16}, samc::kernal_point_t{-1, 2, 4},
    samc::kernal_point_t{ 0,-2, 6}, samc::kernal_point_t{ 0,-1,24}, samc::kernal_point_t{ 0, 0,36}, samc::kernal_point_t{ 0, 1,24}, samc::kernal_point_t{ 0, 2, 6},
    samc::kernal_point_t{ 1,-2, 4}, samc::kernal_point_t{ 1,-1,16}, samc::kernal_point_t{ 1, 0,24}, samc::kernal_point_t{ 1, 1,16}, samc::kernal_point_t{ 1, 2, 4},
    samc::kernal_point_t{ 2,-2, 1}, samc::kernal_point_t{ 2,-1, 4}, samc::kernal_point_t{ 2, 0, 6}, samc::kernal_point_t{ 2, 1, 4}, samc::kernal_point_t{ 2, 2, 1}
    >(*ca, permiability, death_rate);
  Rcpp::XPtr<samc::cache<25>> xp(ca, true);
  return xp;
}

// [[Rcpp::export(.cache_samc_binomial_6)]]
Rcpp::XPtr<samc::cache<49>> cache_samc_binomial_6(const Rcpp::NumericMatrix& permiability, const Rcpp::NumericMatrix& death_rate) {
  samc::cache<49>* ca = new samc::cache<49>;
  samc::construct_cache<
    samc::kernal_point_t{-3,-3, 1}, samc::kernal_point_t{-3,-2,  6}, samc::kernal_point_t{-3,-1, 15}, samc::kernal_point_t{-3, 0, 20}, samc::kernal_point_t{-3, 1, 15}, samc::kernal_point_t{-3, 2,  6},  samc::kernal_point_t{-3, 3, 1},
    samc::kernal_point_t{-2,-3, 6}, samc::kernal_point_t{-2,-2, 36}, samc::kernal_point_t{-2,-1, 90}, samc::kernal_point_t{-2, 0,120}, samc::kernal_point_t{-2, 1, 90}, samc::kernal_point_t{-2, 2, 36},  samc::kernal_point_t{-2, 3, 6},
    samc::kernal_point_t{-1,-3,15}, samc::kernal_point_t{-1,-2, 90}, samc::kernal_point_t{-1,-1,225}, samc::kernal_point_t{-1, 0,300}, samc::kernal_point_t{-1, 1,225}, samc::kernal_point_t{-1, 2, 90},  samc::kernal_point_t{-1, 3,15},
    samc::kernal_point_t{ 0,-3,20}, samc::kernal_point_t{ 0,-2,120}, samc::kernal_point_t{ 0,-1,300}, samc::kernal_point_t{ 0, 0,400}, samc::kernal_point_t{ 0, 1,300}, samc::kernal_point_t{ 0, 2,120},  samc::kernal_point_t{ 0, 3,20},
    samc::kernal_point_t{ 1,-3,15}, samc::kernal_point_t{ 1,-2, 90}, samc::kernal_point_t{ 1,-1,225}, samc::kernal_point_t{ 1, 0,300}, samc::kernal_point_t{ 1, 1,225}, samc::kernal_point_t{ 1, 2, 90},  samc::kernal_point_t{ 1, 3,15},
    samc::kernal_point_t{ 2,-3, 6}, samc::kernal_point_t{ 2,-2, 36}, samc::kernal_point_t{ 2,-1, 90}, samc::kernal_point_t{ 2, 0,120}, samc::kernal_point_t{ 2, 1, 90}, samc::kernal_point_t{ 2, 2, 36},  samc::kernal_point_t{ 2, 3, 6},
    samc::kernal_point_t{ 3,-3, 1}, samc::kernal_point_t{ 3,-2,  6}, samc::kernal_point_t{ 3,-1, 15}, samc::kernal_point_t{ 3, 0, 20}, samc::kernal_point_t{ 3, 1, 15}, samc::kernal_point_t{ 3, 2,  6},  samc::kernal_point_t{ 3, 3, 1}
    >(*ca, permiability, death_rate);
  Rcpp::XPtr<samc::cache<49>> xp(ca, true);
  return xp;
}
*/
