/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#ifndef OOPS_UTIL_RANDOM_H_
#define OOPS_UTIL_RANDOM_H_

#include <string>
#include <vector>
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
 * \sa util::UniformDistribution util::NormalDistribution
 */


template <typename T>
class Random : public util::Printable {
 public:

  void SetSeed(const unsigned int seed =
               static_cast<std::uint32_t>(std::time(0))) {seed_ = seed;}

  const T & operator[](const std::size_t ii) const {return data_[ii];}

 protected:
  Random(size_t N, unsigned int seed): N_(N), seed_(seed) {}
  virtual ~Random();

  std::size_t N_;
  std::vector<T> data_;
  unsigned int seed_;

 private:
  virtual void print(std::ostream &) const = 0;
  
};

// -----------------------------------------------------------------------------

template <typename T>
class UniformDistribution : Random<T> {

 public:
  /*! principal constructor */
  UniformDistribution(size_t, T, T, unsigned int);

  /*! constructors with default values */
  UniformDistribution(size_t N, T minv, T maxv):
  UniformDistribution(N, minv, maxv, static_cast<std::uint32_t>(std::time(0))) {}

  UniformDistribution(size_t N):
  UniformDistribution(N, 0.0, 1.0, static_cast<std::uint32_t>(std::time(0))) {}

  UniformDistribution():
  UniformDistribution(1, 0.0, 1.0, static_cast<std::uint32_t>(std::time(0))) {}

  ~UniformDistribution();

 private:
  T minv_;
  T maxv_;
};

// ------------------------------------------------------------------------------
// Normal distributions

template <typename T>
class NormalDistribution : Random<T> {

 public:
  /*! principal constructor */
  NormalDistribution(size_t, T, T, unsigned int);

  /*! constructors with default values */
  NormalDistribution(size_t N, T mean, T sdev):
  NormalDistribution(N, mean, sdev, static_cast<std::uint32_t>(std::time(0))) {}

  NormalDistribution(size_t N):
  NormalDistribution(N, 0.0, 1.0, static_cast<std::uint32_t>(std::time(0))) {}

  NormalDistribution():
  NormalDistribution(1, 0.0, 1.0, static_cast<std::uint32_t>(std::time(0))) {}

  ~NormalDistribution();

 private:
  T mean_;
  T sdev_;
};

// ------------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_RANDOM_H_
