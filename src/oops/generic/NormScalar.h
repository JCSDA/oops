/*
 * (C) Crown Copyright 2023, the Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#ifndef OOPS_GENERIC_NORMSCALAR_H_
#define OOPS_GENERIC_NORMSCALAR_H_

#include <string>

#include "oops/base/Increment.h"
#include "oops/base/NormBase.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL>
class NormScalar : public NormBase<MODEL> {
 public:
  static const std::string classname() {return "oops::NormScalar";}

  explicit NormScalar(const eckit::Configuration &);
  ~NormScalar() {}

  void multiplyMatrix(Increment<MODEL>&);
  void multiplyMatrixInverse(Increment<MODEL>&);

 private:
  double alpha_;
};

// =============================================================================

template<typename MODEL>
NormScalar<MODEL>::NormScalar(const eckit::Configuration & conf)
{
  alpha_ = conf.getDouble("alpha", 1.0);
}

template<typename MODEL>
void NormScalar<MODEL>::multiplyMatrix(Increment<MODEL> &dx) {
  dx *= alpha_;
}

template<typename MODEL>
void NormScalar<MODEL>::multiplyMatrixInverse(Increment<MODEL> &dx) {
  dx *= 1.0/alpha_;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_GENERIC_NORMSCALAR_H_
