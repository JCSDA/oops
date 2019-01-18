/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#ifndef OOPS_UTIL_RANDOM_H_
#define OOPS_UTIL_RANDOM_H_

#include <iostream>
#include <string>
#include <vector>
#include <boost/random.hpp>
#include "oops/util/abor1_cpp.h"
#include "oops/util/formats.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------
/*! Base Class for generating compiler-indpependent random numbers
 *
 * \details The *Random* class and its derived classes are intended to provide 
 * a compiler independent pseudo random number generator.  So, if you use the
 * same random seed, you should get the same results, regardless of compiler
 * or platform.  Its sub-classes implement specific distributions, including
 * uniform and Gaussian (Normal) deviates.  
 * For usage tips, see each individual sub-class.
 *
 * \date Jan, 2019 (M. Miesch, created)
 *
 * \sa util::UniformDistribution, util::UniformIntDistribution, 
 * util::NormalDistribution
 */

template <typename T>
class Random : public util::Printable {
 public:
  const T & operator[](const std::size_t ii) const {return data_[ii];}
  std::vector<T> data() {return data_;}

 protected:
  Random(size_t N, unsigned int seed): N_(N), seed_(seed) {}
  virtual ~Random() {}

  std::size_t N_;
  std::vector<T> data_;
  unsigned int seed_;

 private:
  /*! This prints in a format that can be easily inserted into a yaml file for testing */
  void print(std::ostream & os) const {
    for (size_t jj=0; jj < N_; ++jj) {
      os << "   - " << util::full_precision(data_[jj]) << std::endl;
    }
  }
};

// -----------------------------------------------------------------------------
/*! Class for generating uniformly-distributed random numbers
 *
 * \details *util::UniformDistribution* creates a vector of psedo-random real
 * numbers that are uniformly distributed across a specified interval of values.
 *
 * \param[in] N The size of the desired array (default 1)
 * \param[in] minv The minimum value of the interval (default 0)
 * \param[in] minv The maximum value of the interval (default 1)
 * \param[in] optional seed to use for the random number generator.  If omitted,
 *            the code will define a seed based on the current (calendar) time.
 * 
 * \example Example usage:
 * util::UniformDistribution<double> x(N,1.0,100.0)  
 * std::cout << x[i] << std::endl;  // access one element
 * std::cout << x << std::endl;  // print full array
 *  
 * \note The interval is closed on the lower end and open on the upper end, i.e. [minv,maxv).  
 *
 * \warning Only implemented for floating point data types.  For a sequence of random 
 * integers use **util::UniformIntDistribution()**
 *
 */

template <typename T>
class UniformDistribution : public Random<T> {
 public:
  UniformDistribution(std::size_t N = 1, T minv = 0, T maxv = 1,
                      unsigned int seed = static_cast<std::uint32_t>(std::time(0))):
  Random<T>(N, seed), minv_(minv), maxv_(maxv) {
    static boost::random::mt19937 generator(this->seed_);
    boost::random::uniform_real_distribution<T> distribution(minv_, maxv_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }

  virtual ~UniformDistribution() {}

 private:
  T minv_;
  T maxv_;
};

// -----------------------------------------------------------------------------
/*! Class for generating uniformly-distributed random integers
 *
 * \details *util::UniformIntDistribution* creates a vector of psedo-random integers
 * that are uniformly distributed across a specified interval of values.
 *
 * \param[in] N The size of the desired array (default 1)
 * \param[in] minv The minimum value of the interval (default 0)
 * \param[in] minv The maximum value of the interval (default 100)
 * \param[in] optional seed to use for the random number generator.  If omitted,
 *            the code will define a seed based on the current (calendar) time.
 * 
 * \example Example usage:
 * util::UniformIntDistribution<int> x(N,1,100)  
 * std::cout << x[i] << std::endl;  // access one element
 * std::cout << x << std::endl;  // print full array
 *  
 * \note the interval is closed on both ends, i.e. [minv,maxv].  So, calling the integer 
 * constructor with vmin = 1 and vmax = 6 will return a random integer number between 1 and 6.
 *
 */

template <typename T>
class UniformIntDistribution : public Random<T> {
 public:
  UniformIntDistribution(std::size_t N = 1, T minv = 0, T maxv = 100,
                         unsigned int seed = static_cast<std::uint32_t>(std::time(0))):
  Random<T>(N, seed), minv_(minv), maxv_(maxv) {
    static boost::random::mt19937 generator(this->seed_);
    boost::random::uniform_int_distribution<T> distribution(minv_, maxv_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }

  virtual ~UniformIntDistribution() {}

 private:
  T minv_;
  T maxv_;
};

// ------------------------------------------------------------------------------
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
 * \example Example usage:
 * util::NormalDistribution<double> x(N,0.0,20.0)  
 * std::cout << x[i] << std::endl;  // access one element
 * std::cout << x << std::endl;  // print full array
 * 
 * \warning Only implemented for floating point data types
 */

template <typename T>
class NormalDistribution : public Random<T> {
 public:
  NormalDistribution(std::size_t N = 1, T mean = 0, T sdev = 1,
                     unsigned int seed = static_cast<std::uint32_t>(std::time(0))):
  Random<T>(N, seed), mean_(mean), sdev_(sdev) {
    if (!std::is_floating_point<T>::value) {
      oops::Log::error() << "NormalDistribution only implemented for floating point data types"
                         << std::endl;
      ABORT("NormalDistribution only implemented for floating point data types");
    }
    static boost::random::mt19937 generator(this->seed_);
    boost::random::normal_distribution<T> distribution(mean_, sdev_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }

  virtual ~NormalDistribution() {}

 private:
  T mean_;
  T sdev_;
};

// ------------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_RANDOM_H_
