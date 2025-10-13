/*
* Copyright 2011 ECMWF
* Copyright 2020-2021 UCAR
*
* This software was developed at ECMWF for evaluation
* and may be used for academic and research purposes only.
* The software is provided as is without any warranty.
*
* This software can be used, copied and modified but not
* redistributed or sold. This notice must be reproduced
* on each copy made.
*/

#pragma once

#include <memory>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/generic/LocalizationBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {

namespace interface {

// -----------------------------------------------------------------------------

template<typename MODEL>
class Localization: public oops::LocalizationBase<MODEL> {
  typedef typename MODEL::Localization    Localization_;
  typedef oops::Geometry<MODEL>           Geometry_;
  typedef oops::Increment<MODEL>          Increment_;
 public:
  static const std::string classname() {return "oops::Localization";}

  Localization(const Geometry_ &, const oops::Variables &, const eckit::Configuration &);
  virtual ~Localization();

  /// Overrides for oops::LocalizationBase classes, passing MODEL-specific classes to the
  /// MODEL-specific implementations of Localization
  void randomize(Increment_ & dx) const override;
  void multiply(Increment_ & dx) const override;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<Localization_> loc_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & resol, const oops::Variables & vars,
                                  const eckit::Configuration & conf)
  : loc_()
{
  Log::trace() << "Localization<MODEL>::Localization starting" << std::endl;
  util::Timer timer(classname(), "Localization");
  loc_.reset(new Localization_(resol.geometry(), vars, conf));
  Log::trace() << "Localization<MODEL>::Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  Log::trace() << "Localization<MODEL>::~Localization starting" << std::endl;
  util::Timer timer(classname(), "~Localization");
  loc_.reset();
  Log::trace() << "Localization<MODEL>::~Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::randomize(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  loc_->randomize(dx.increment());
  Log::trace() << "Localization<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  loc_->multiply(dx.increment());
  Log::trace() << "Localization<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Localization<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *loc_;
  Log::trace() << "Localization<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace interface

}  // namespace oops

