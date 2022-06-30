/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_NORMGRADIENT_H_
#define OOPS_INTERFACE_NORMGRADIENT_H_

#include <string>

#include <memory>
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// =============================================================================

template <typename MODEL>
class NormGradient : public util::Printable,
                     private util::ObjectCounter<NormGradient<MODEL> > {
  typedef typename MODEL::NormGradient  NormGradient_;
  typedef Geometry<MODEL>               Geometry_;
  typedef Increment<MODEL>              Increment_;
  typedef State<MODEL>                  State_;

 public:
  static const std::string classname() {return "oops::NormGradient";}

// Constructor, destructor
  NormGradient(const Geometry_ &, const State_ &, const eckit::Configuration &);
  virtual ~NormGradient();

// Compute gradient of some norm J(increment)
  void apply(Increment_ &) const;

// Interfacing
  const NormGradient_ & normgradient() const {return *normgradient_;}
  NormGradient_ & normgradient() {return *normgradient_;}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<NormGradient_> normgradient_;
};

// =============================================================================

template<typename MODEL>
NormGradient<MODEL>::NormGradient(const Geometry_ & resol, const State_ & xr,
                                  const eckit::Configuration & conf) : normgradient_()
{
  Log::trace() << "NormGradient<MODEL>::NormGradient starting" << std::endl;
  util::Timer timer(classname(), "NormGradient");
  normgradient_.reset(new NormGradient_(resol.geometry(), xr.state(), conf));
  Log::trace() << "NormGradient<MODEL>::NormGradient done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
NormGradient<MODEL>::~NormGradient() {
  Log::trace() << "NormGradient<MODEL>::~NormGradient starting" << std::endl;
  util::Timer timer(classname(), "~NormGradient");
  normgradient_.reset();
  Log::trace() << "NormGradient<MODEL>::~NormGradient done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void NormGradient<MODEL>::apply(Increment_ & dx) const {
  Log::trace() << "NormGradient<MODEL>::apply starting" << std::endl;
  util::Timer timer(classname(), "apply");
  normgradient_->apply(dx.increment());
  Log::trace() << "NormGradient<MODEL>::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void NormGradient<MODEL>::print(std::ostream & os) const {
  Log::trace() << "NormGradient<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *normgradient_;
  Log::trace() << "NormGradient<MODEL>::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_NORMGRADIENT_H_
