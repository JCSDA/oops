/*
 * (C) Crown Copyright 2026, Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#ifndef OOPS_GENERIC_NORMSCALAR_H_
#define OOPS_GENERIC_NORMSCALAR_H_

#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/NormBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL>
class NormScalar : public NormBase<MODEL> {
  typedef Increment<MODEL>              Increment_;
  typedef State<MODEL>                  State_;

 public:
  static const std::string classname() {return "oops::NormScalar";}

  NormScalar(const Variables &,
             const eckit::Configuration &,
             const State_ &);
  ~NormScalar() = default;

  void computeGradient(Increment_ &) const override;
  void computeGradientInverse(Increment_ &) const override;
  const double computeNorm(const Increment_ &) const override;

 private:
  void print(std::ostream &) const override;
  double alpha_;
};

// =============================================================================

template<typename MODEL>
NormScalar<MODEL>::NormScalar(const Variables &,
                              const eckit::Configuration & conf,
                              const State_ &)
{
  alpha_ = conf.getDouble("alpha", 1.0);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void NormScalar<MODEL>::computeGradient(Increment_ & dx) const {
  dx *= alpha_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void NormScalar<MODEL>::computeGradientInverse(Increment_ & dx) const {
  dx *= 1.0/alpha_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const double NormScalar<MODEL>::computeNorm(const Increment_ & dx) const {
  return alpha_*dot_product(dx, dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void NormScalar<MODEL>::print(std::ostream & os) const {
  os << "NormScalar<MODEL> alpha_ = " << alpha_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_GENERIC_NORMSCALAR_H_
