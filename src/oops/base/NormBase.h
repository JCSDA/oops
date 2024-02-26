/*
 * (C) Crown Copyright 2023, the Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_NORMBASE_H_
#define OOPS_BASE_NORMBASE_H_


#include <map>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/base/Increment.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL>
class NormBase {
 public:
  static const std::string classname() {return "oops::Norm";}

  NormBase() {}
  virtual ~NormBase() {}

  /* multiplyMatrix: the multiplication of the input increment dx by the matrix that gives rise to the inner
   * product that induces this norm. If norm(x) = sqrt(x^T A x), then this method computes
   * A dx. */
  virtual void multiplyMatrix(Increment<MODEL>&) = 0;

  /* multiplyMatrixInverse: the multiplication of the input increment dx by the inverse of the matrix that gives rise
   * to the inner product that induces this norm. If norm(x) = sqrt(x^T A x), then this method
   * computes A^(-1) dx. */
  virtual void multiplyMatrixInverse(Increment<MODEL>&) = 0;
};

/// \brief factory
template <typename MODEL>
class NormFactory {
 public:
  static NormBase<MODEL> * create(const eckit::Configuration &);
  virtual ~NormFactory() = default;
 protected:
  explicit NormFactory(const std::string &);
 private:
  virtual NormBase<MODEL> * make(const eckit::Configuration &) = 0;
  static std::map < std::string, NormFactory<MODEL> * > & getMakers() {
    static std::map < std::string, NormFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class NormMaker : public NormFactory<MODEL> {
  virtual NormBase<MODEL> * make(const eckit::Configuration & conf)
    { return new T(conf); }
 public:
  explicit NormMaker(const std::string & name) : NormFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
NormFactory<MODEL>::NormFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in Jc Matrix factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
NormBase<MODEL>* NormFactory<MODEL>::create(
                             const eckit::Configuration & config) {
  std::string id = config.getString("norm type");
  Log::info() << "Norm Type = " << id << std::endl;
  typename std::map<std::string, NormFactory<MODEL>*>::iterator j = getMakers().find(id);
  if (j == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in Norm factory.");
  }
  return (*j).second->make(config);
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_BASE_NORMBASE_H_
