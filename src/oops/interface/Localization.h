/*
* Copyright 2011 ECMWF
* Copyright 2020-2020 UCAR
*
* This software was developed at ECMWF for evaluation
* and may be used for academic and research purposes only.
* The software is provided as is without any warranty.
*
* This software can be used, copied and modified but not
* redistributed or sold. This notice must be reproduced
* on each copy made.
*/

#ifndef OOPS_INTERFACE_LOCALIZATION_H_
#define OOPS_INTERFACE_LOCALIZATION_H_

#include <memory>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LocalizationBase.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Model-space localization class: intended for model-specific implementations
template<typename MODEL, typename LOC>
class Localization : public LocalizationBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
 public:
  static const std::string classname() {return "oops::Localization";}

  Localization(const Geometry_ &, const util::DateTime &, const eckit::Configuration &);
  ~Localization();

  void doRandomize(Increment_ &) const override;
  void doMultiply(Increment_ &) const override;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<LOC> loc_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename LOC>
Localization<MODEL, LOC>::Localization(const Geometry_ & geometry,
                                       const util::DateTime & time,
                                       const eckit::Configuration & conf) {
  Log::trace() << "Localization<MODEL>::Localization starting" << std::endl;
  util::Timer timer(classname(), "Localization");
  loc_.reset(new LOC(geometry.geometry(), conf));
  Log::trace() << "Localization<MODEL>::Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename LOC>
Localization<MODEL, LOC>::~Localization() {
  Log::trace() << "Localization<MODEL>::~Localization starting" << std::endl;
  util::Timer timer(classname(), "~Localization");
  loc_.reset();
  Log::trace() << "Localization<MODEL>::~Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename LOC>
void Localization<MODEL, LOC>::doRandomize(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  loc_->randomize(dx.increment());
  Log::trace() << "Localization<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename LOC>
void Localization<MODEL, LOC>::doMultiply(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  loc_->multiply(dx.increment());
  Log::trace() << "Localization<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename LOC>
void Localization<MODEL, LOC>::print(std::ostream & os) const {
  Log::trace() << "Localization<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *loc_;
  Log::trace() << "Localization<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCALIZATION_H_
