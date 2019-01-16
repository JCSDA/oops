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
                            const float & maxv, int32_t & seed, float* vec) {
  util::UniformDistribution<float> x(N, minv, maxv, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}
			      
void random_uniform_double_f(const std::size_t & N, const double & minv, 
                             const double & maxv, int32_t & seed, double* vec) {
  util::UniformDistribution<double> x(N, minv, maxv, static_cast<unsigned int>(seed));
  for (std::size_t jj = 0; jj < N; ++jj) vec[jj] = x[jj]; 
}
			      
// -----------------------------------------------------------------------------

}  // namespace util
  
