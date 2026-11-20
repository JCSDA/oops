/*
 * (C) Crown Copyright 2026, Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_NORM_H_
#define OOPS_INTERFACE_NORM_H_

#include <memory>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Increment.h"
#include "oops/base/NormBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {
namespace interface {

// =============================================================================

template <typename MODEL, typename IMPL>
class Norm : public oops::NormBase<MODEL> {
  typedef oops::Increment<MODEL>        Increment_;
  typedef oops::State<MODEL>            State_;

 public:
  static const std::string classname() {return "oops::interface::Norm";}

// Constructor, destructor
  Norm(const Variables &, const eckit::Configuration &, const State_ &);
  virtual ~Norm();

// Apply norm to increment
  void computeGradient(Increment_ &) const override;
  void computeGradientInverse(Increment_ &) const override;
  const double computeNorm(const Increment_ &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<IMPL> norm_;
};

// =============================================================================

template <typename MODEL, typename IMPL>
Norm<MODEL, IMPL>::Norm(const Variables & vars,
                        const eckit::Configuration & conf,
                        const State_ & state) : norm_(nullptr)
{
  Log::trace() << "Norm<MODEL, IMPL>::Norm starting" << std::endl;
  util::Timer timer(classname(), "Norm");
  norm_.reset(new IMPL(vars, conf, state.state()));
  Log::trace() << "Norm<MODEL, IMPL>::Norm done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename IMPL>
Norm<MODEL, IMPL>::~Norm() {
  Log::trace() << "Norm<MODEL, IMPL>::~Norm starting" << std::endl;
  util::Timer timer(classname(), "~Norm");
  norm_.reset();
  Log::trace() << "Norm<MODEL, IMPL>::~Norm done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename IMPL>
void Norm<MODEL, IMPL>::computeGradient(Increment_ & dx) const {
  Log::trace() << "Norm<MODEL, IMPL>::computeGradient starting" << std::endl;
  util::Timer timer(classname(), "computeGradient");
  norm_->computeGradient(dx.increment());
  Log::trace() << "Norm<MODEL, IMPL>::computeGradient done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename IMPL>
void Norm<MODEL, IMPL>::computeGradientInverse(Increment_ & dx) const {
  Log::trace() << "Norm<MODEL, IMPL>::computeGradientInverse starting" << std::endl;
  util::Timer timer(classname(), "computeGradientInverse");
  norm_->computeGradientInverse(dx.increment());
  Log::trace() << "Norm<MODEL, IMPL>::computeGradientInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename IMPL>
const double Norm<MODEL, IMPL>::computeNorm(const Increment_ & dx) const {
  Log::trace() << "Norm<MODEL, IMPL>::computeNorm starting" << std::endl;
  util::Timer timer(classname(), "computeNorm");
  const double zz = norm_->computeNorm(dx.increment());
  Log::trace() << "Norm<MODEL, IMPL>::computeNorm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename IMPL>
void Norm<MODEL, IMPL>::print(std::ostream & os) const {
  os << *norm_;
}

// -----------------------------------------------------------------------------

}  // namespace interface
}  // namespace oops

#endif  // OOPS_INTERFACE_NORM_H_
