/*
 * (C) Copyright 2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#include <string>
#include <vector>
#include <boost/random.hpp>
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

namespace util {

// -----------------------------------------------------------------------------
// Uniform random distribution

template <typename T>
UniformDistribution<T>::UniformDistribution(size_t N, T minv, T maxv, int seed):
  Random<T>(N, seed), minv_(minv), maxv_(maxv) {
  boost::random::mt19937 generator(this->seed_);

  if (std::is_integral<T>::value) {
    boost::random::uniform_int_distribution<T> distribution(minv_, maxv_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  } else {
    boost::random::uniform_real_distribution<T> distribution(minv_, maxv_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }
}

// -----------------------------------------------------------------------------
// Normal random distribution

template <typename T>
NormalDistribution<T>::NormalDistribution(size_t N, T mean, T sdev, int seed):
  Random<T>(N, seed), mean_(mean), sdev_(sdev) {
  if (!std::is_floating_point<T>::value) {
    oops::Log::error() << "NormalDistribution only implemented for floating point data types"
                       << std::endl;
    ABORT("NormalDistribution only implemented for floating point data types");
  }

  boost::random::mt19937 generator(this->seed_);
  boost::random::normal_distribution<T> distribution(mean_, sdev_);
  for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
}

// -----------------------------------------------------------------------------

}  // namespace util
