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
/*! Class for generating uniformly-distributed random numbers
 *
 * \details *util::UniformDistribution* creates a vector of psedo-random numbers
 * that are uniformly distributed across a specified interval of values.
 *
 * \param[in] N The size of the desired array (default 1)
 * \param[in] minv The minimum value of the interval (default 0)
 * \param[in] minv The maximum value of the interval (default 1)
 * \param[in] optional seed to use for the random number generator.  If omitted,
 *            the code will define a seed based on the current (calendar) time.
 * 
 * \note If the data type is real, the interval is closed on the lower end and open
 * on the upper end, i.e. [minv,maxv).  However, if the type is integer the interval
 * is closed on both ends, i.e. [minv,maxv].  So, calling the integer constructor 
 * with vmin = 1 and vmax = 6 will return a random integer number between 1 and 6.
 *
 */

template <typename T>
UniformDistribution<T>::UniformDistribution(size_t N, T minv, T maxv, unsigned int seed):
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
/*! Class for generating Gaussian-distributed random numbers
 *
 * \details *util::NormalDistribution* creates a vector of psedo-random numbers
 * with a normal (Gaussian) distribution.
 *
 * \param[in] N The size of the desired array (default 1)
 * \param[in] mean The mean of the distribution (default 0)
 * \param[in] sdev The standard deviation of the destribution (default 1)
 * \param[in] seed seed to use for the random number generator.  If omitted,
 *            a seed will be generated based on the current (calendar) time.
 * 
 */

template <typename T>
NormalDistribution<T>::NormalDistribution(size_t N, T mean, T sdev, unsigned int seed):
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
