/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#ifndef OOPS_UTIL_RANDOM_H_
#define OOPS_UTIL_RANDOM_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <boost/random.hpp>
#include "eckit/exception/Exceptions.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/formats.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------
/*! Base Class for generating compiler-independent random numbers
 *
 * \details The *Random* class and its derived classes are intended to provide
 * a compiler independent pseudo random number generator.  So, if you use the
 * same random seed, you should get the same results, regardless of compiler
 * or platform.  Its sub-classes implement specific distributions, including
 * uniform and Gaussian (Normal) deviates.
 * For usage tips, see each individual sub-class.
 *
 * \date Jan, 2019 (M. Miesch, created)
 * \date Jan, 2020 (added shuffle class and reset option)
 *
 * \sa util::UniformDistribution, util::UniformIntDistribution,
 * util::NormalDistribution
 */

template <typename datatype>
class Random : public util::Printable {
 public:
  const datatype & operator[](const std::size_t ii) const {return data_[ii];}
  std::vector<datatype> data() {return data_;}
  void sort() {std::sort(data_.begin(), data_.end());}

 protected:
  explicit Random(size_t N): N_(N) {data_.reserve(N);}
  virtual ~Random() {}

  std::size_t N_;
  std::vector<datatype> data_;

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
 * \details *util::UniformDistribution* creates a vector of pseudo-random real
 * numbers that are uniformly distributed across a specified interval of values.
 *
 * \param[in] N The size of the desired array (default 1)
 * \param[in] minv The minimum value of the interval (default 0)
 * \param[in] minv The maximum value of the interval (default 1)
 * \param[in] seed The seed to use for the random number generator upon first
 *            instantiation of the class (optional).  If omitted, a seed
 *            will be generated based on the calendar time.  This parameter is 
 *            only used on the first instantiation of this class for a
 *            particular data type and/or if the **reset** flag is set to true.
 * \param[in] reset If this is set to **true** then this forces the generator
 *            to re-initilize itself with the specified seed.  Otherwise, the
 *            seed is only used on first instantiation (this is the default
 *            behavior).
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

template <typename datatype>
class UniformDistribution : public Random<datatype> {
 public:
  UniformDistribution(std::size_t N = 1, datatype minv = 0, datatype maxv = 1,
                         unsigned int seed = static_cast<std::uint32_t>(std::time(nullptr)),
                         bool reset = false): Random<datatype>(N), minv_(minv), maxv_(maxv) {
    static boost::random::mt19937 generator(seed);
    if (reset) generator.seed(seed);
    boost::random::uniform_real_distribution<datatype> distribution(minv_, maxv_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }

  virtual ~UniformDistribution() {}

 private:
  datatype minv_;
  datatype maxv_;
};

// -----------------------------------------------------------------------------
/*! Class for generating uniformly-distributed random integers
 *
 * \details *util::UniformIntDistribution* creates a vector of pseudo-random integers
 * that are uniformly distributed across a specified interval of values.
 *
 * \param[in] N The size of the desired array (default 1)
 * \param[in] minv The minimum value of the interval (default 0)
 * \param[in] minv The maximum value of the interval (default 100)
 * \param[in] seed The seed to use for the random number generator upon first
 *            instantiation of the class (optional).  If omitted, a seed
 *            will be generated based on the calendar time.  This parameter is 
 *            only used on the first instantiation of this class for a
 *            particular data type and/or if the **reset** flag is set to true.
 * \param[in] reset If this is set to **true** then this forces the generator
 *            to re-initilize itself with the specified seed.  Otherwise, the
 *            seed is only used on first instantiation (this is the default
 *            behavior).
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

template <typename datatype>
class UniformIntDistribution : public Random<datatype> {
 public:
  UniformIntDistribution(std::size_t N = 1, datatype minv = 0, datatype maxv = 100,
                         unsigned int seed = static_cast<std::uint32_t>(std::time(nullptr)),
                         bool reset = false): Random<datatype>(N), minv_(minv), maxv_(maxv) {
    static boost::random::mt19937 generator(seed);
    if (reset) generator.seed(seed);
    boost::random::uniform_int_distribution<datatype> distribution(minv_, maxv_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }

  virtual ~UniformIntDistribution() {}

 private:
  datatype minv_;
  datatype maxv_;
};

// ------------------------------------------------------------------------------
/*! Class for generating Gaussian-distributed random numbers
 *
 * \details *util::NormalDistribution* creates a vector of pseudo-random numbers
 * with a normal (Gaussian) distribution.
 *
 * \param[in] N The size of the desired array (default 1)
 * \param[in] mean The mean of the distribution (default 0)
 * \param[in] sdev The standard deviation of the destribution (default 1)
 * \param[in] seed The seed to use for the random number generator upon first
 *            instantiation of the class (optional).  If omitted, a seed
 *            will be generated based on the calendar time.  This parameter is 
 *            only used on the first instantiation of this class for a
 *            particular data type and/or if the **reset** flag is set to true.
 * \param[in] reset If this is set to **true** then this forces the generator
 *            to re-initilize itself with the specified seed.  Otherwise, the
 *            seed is only used on first instantiation (this is the default
 *            behavior).
 *
 * \example Example usage:
 * util::NormalDistribution<double> x(N,0.0,20.0)
 * std::cout << x[i] << std::endl;  // access one element
 * std::cout << x << std::endl;  // print full array
 *
 * \warning Only implemented for floating point data types
 */

template <typename datatype>
class NormalDistribution : public Random<datatype> {
 public:
  NormalDistribution(std::size_t N = 1, datatype mean = 0, datatype sdev = 1,
                     unsigned int seed = static_cast<std::uint32_t>(std::time(nullptr)),
                     bool reset = false): Random<datatype>(N), mean_(mean), sdev_(sdev) {
    if (!std::is_floating_point<datatype>::value) {
      throw eckit::BadCast("NormalDistribution only implemented for floating point data types");
    }
    static boost::random::mt19937 generator(seed);
    if (reset) generator.seed(seed);
    boost::random::normal_distribution<datatype> distribution(mean_, sdev_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }

  virtual ~NormalDistribution() {}

 private:
  datatype mean_;
  datatype sdev_;
};

// -----------------------------------------------------------------------------
/*! \brief Shuffle (reorder randomly) a range of elements.
 *
 * The C++ standard doesn't require all STL implementations of std::shuffle to produce the same
 * results. This function provides a concrete implementation based on an example given on
 * https://en.cppreference.com/w/cpp/algorithm/random_shuffle, itself very similar to libc++'s
 * implementation.
 * \param[in] begin
 *   A random-access iterator addressing the position of the first element in the
 *   range to be shuffled.
 * \param[in] end
 *   A random-access iterator addressing the position one past the last element
 *  in the range to be shuffled.
 * \param[in] seed The seed to use for the random number generator upon first
 *            instantiation of the function (optional).  If omitted, a seed
 *            will be generated based on the calendar time.  This parameter is
 *            only used on the first instantiation of this class for a
 *            particular data type and/or if the **reset** flag is set to true.
 * \param[in] reset If this is set to **true** then this forces the generator
 *            to re-initialize itself with the specified seed.  Otherwise, the
 *            seed is only used on first instantiation (this is the default
 *            behavior).
 */
template<class RandomIt>
void shuffle(RandomIt begin, RandomIt end, unsigned int seed =
             static_cast<std::uint32_t>(std::time(nullptr)), bool reset = false) {
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
    typedef boost::random::uniform_int_distribution<diff_t> distr_t;
    typedef typename distr_t::param_type param_t;

    static boost::random::mt19937 generator(seed);
    if (reset) generator.seed(seed);

    distr_t distribution;
    diff_t n = end - begin;
    for (diff_t i = n - 1; i > 0; --i) {
        std::swap(begin[i], begin[distribution(generator, param_t(0, i))]);
    }
}

// ------------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_RANDOM_H_
