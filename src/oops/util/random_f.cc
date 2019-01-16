/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#include "oops/util/random_f.h"
#include "oops/util/Random.h"

namespace util {

// -----------------------------------------------------------------------------
void random_uniform_float_f(const std::size_t & N, const float & minv, 
                            const float & maxv, std::int32_t & seed, float* vec) {
  util::UniformDistribution<float> x(N, minv, maxv, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}
			      
void random_uniform_double_f(const std::size_t & N, const double & minv, 
                             const double & maxv, std::int32_t & seed, double* vec) {
  util::UniformDistribution<double> x(N, minv, maxv, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}
			      
void random_uniform_int_f(const std::size_t & N, const std::int32_t & minv, 
                          const std::int32_t & maxv, std::int32_t & seed,
                          std::int32_t* vec) {
  util::UniformIntDistribution<std::int32_t> x(N, minv, maxv, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}
			      
void random_uniform_long_f(const std::size_t & N, const std::int64_t & minv, 
                           const std::int64_t & maxv, std::int32_t & seed,
                           std::int64_t* vec) {
  util::UniformIntDistribution<std::int64_t> x(N, minv, maxv, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}
			      
void random_normal_float_f(const std::size_t & N, const float & mean, 
                           const float & sdev, std::int32_t & seed, float* vec) {
  util::NormalDistribution<float> x(N, mean, sdev, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}

void random_normal_double_f(const std::size_t & N, const double & mean, 
                            const double & sdev, std::int32_t & seed, double* vec) {
  util::NormalDistribution<float> x(N, mean, sdev, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}

// -----------------------------------------------------------------------------

}  // namespace util
  
