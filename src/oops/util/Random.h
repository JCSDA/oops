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
/// Classes for generating compiler-indpependent random numbers

template <typename T>
class Random : public util::Printable {
 public:
  virtual const std::string classname();

  void newseed(const int seed = 77) {seed_ = seed;}

  const T & operator[](const std::size_t ii) const {return data_[ii];}

  virtual void print(std::ostream) const;

 protected:
  Random(size_t N, int seed): N_(N), seed_(seed) {}
  ~Random();

  std::size_t N_;
  std::vector<T> data_;
  int seed_;
};

// -----------------------------------------------------------------------------
// Uniform distributions

template <typename T>
class UniformDistribution : Random<T> {
  const std::string classname() {return "util::UniformDistribution";}

 public:
  UniformDistribution(size_t, T, T, int);
  ~UniformDistribution();

 private:
  T minv_;
  T maxv_;
};

// -----------------------------------------------------------------------------
// Normal distributions

template <typename T>
class NormalDistribution : Random<T> {
  const std::string classname() {return "util::NormalDistribution";}

 public:
  NormalDistribution(size_t, T, T, int);
  ~NormalDistribution();

 private:
  T mean_;
  T sdev_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_RANDOM_H_
