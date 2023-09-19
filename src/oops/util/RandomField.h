/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#ifndef OOPS_UTIL_RANDOMFIELD_H_
#define OOPS_UTIL_RANDOMFIELD_H_

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <boost/random.hpp>
#include "eckit/exception/Exceptions.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/formats.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------
/*! Base Class for generating compiler-independent random numbers for fields
 *
 * \details The *RandomField* class and its derived classes are intended to provide
 * a compiler independent pseudo random number generator.  So, if you use the
 * same random seed, you should get the same results, regardless of compiler
 * or platform.  Its sub-classes implement Gaussian (Normal) deviates.
 * For usage tips, see each individual sub-class.
 *
 * \date Aug, 2023 (copy from Random.h, remove reset option)
 *
 * \sa util::NormalDistribution
 */

class RandomField : public util::Printable {
 public:
  const double & operator[](const size_t ii) const {return data_[ii];}
  std::vector<double> data() {return data_;}

 protected:
  explicit RandomField(size_t N): N_(N) {data_.reserve(N);}
  virtual ~RandomField() {}

  size_t N_;
  std::vector<double> data_;

 private:
  /*! This prints in a format that can be easily inserted into a yaml file for testing */
  void print(std::ostream & os) const {
    for (size_t jj=0; jj < N_; ++jj) {
      os << "   - " << util::full_precision(data_[jj]) << std::endl;
    }
  }
};

// ------------------------------------------------------------------------------
/*! Class for generating Gaussian-distributed random numbers
 *
 * \details *util::NormalDistributionField* creates a vector of pseudo-random numbers
 * with a normal (Gaussian) distribution. Specific seed handling for fields (different
 * on each MPI task).
 *
 * \param[in] N The size of the desired array
 * \param[in] mean The mean of the distribution
 * \param[in] sdev The standard deviation of the destribution
 * \param[in] seed The seed to use for the random number generator (optional).
 *
 * \example Example usage:
 * util::NormalDistributionField x(N,0.0,20.0)
 * std::cout << x[i] << std::endl;  // access one element
 * std::cout << x << std::endl;  // print full array
 */

class NormalDistributionField : public RandomField {
 public:
  NormalDistributionField(size_t N, double mean, double sdev):
    RandomField(N), mean_(mean), sdev_(sdev) {
    static boost::random::mt19937 generator(oops::mpi::world().rank()+1);
    boost::random::normal_distribution<double> distribution(mean_, sdev_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
  }
  NormalDistributionField(size_t N, double mean, double sdev, unsigned int seed):
    RandomField(N), mean_(mean), sdev_(sdev) {
    static boost::random::mt19937 generator(oops::mpi::world().rank()+1);
    std::stringstream ss;
    ss << generator;
    generator.seed(seed);
    boost::random::normal_distribution<double> distribution(mean_, sdev_);
    for (size_t jj=0; jj < this->N_; ++jj) this->data_.push_back(distribution(generator));
    ss >> generator;
  }

  virtual ~NormalDistributionField() {}

 private:
  double mean_;
  double sdev_;
};

// ------------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_RANDOMFIELD_H_
