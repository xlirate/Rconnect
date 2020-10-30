
inline constexpr int _shrink(const int o_x, const int k_x, const int, const int){
  return o_x + k_x;
}

inline constexpr int _zero(const int o_x, const int k_x, const int d_width, const int k_width){
  return o_x + k_x - k_width /2;
}

inline int _reflect(const int o_x, const int k_x, const int d_width, const int k_width){
  auto d_x = abs(o_x + k_x - k_width /2);
  return ((d_x/d_width )%2)?((d_width -1)-((d_x)%d_width )):(d_x%d_width );
}

inline constexpr int _wrap(const int o_x, const int k_x, const int d_width, const int k_width){
  return (o_x + k_x - k_width /2)%(d_width);
}

inline int _stretch(const int o_x, const int k_x, const int d_width, const int k_width){
  return std::clamp(o_x+k_x-k_width/2, 0, d_width -1);
}