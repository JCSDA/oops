/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LOCALIZATION_H_
#define OOPS_INTERFACE_LOCALIZATION_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LocalizationBase.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class Localization : public util::Printable,
                     private boost::noncopyable,
                     private util::ObjectCounter<Localization<MODEL> > {
  typedef LocalizationBase<MODEL>    LocalizationBase_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;

 public:
  static const std::string classname() {return "oops::Localization";}

  Localization(const Geometry_ &, const eckit::Configuration &);
  virtual ~Localization();

  void multiply(Increment_ &) const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<LocalizationBase_> local_;
};

// =============================================================================

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & resol,
                                  const eckit::Configuration & conf) : local_()
{
  Log::trace() << "Localization<MODEL>::Localization starting" << std::endl;
  util::Timer timer(classname(), "Localization");
  local_.reset(LocalizationFactory<MODEL>::create(resol, conf));
  Log::trace() << "Localization<MODEL>::Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  Log::trace() << "Localization<MODEL>::~Localization starting" << std::endl;
  util::Timer timer(classname(), "~Localization");
  local_.reset();
  Log::trace() << "Localization<MODEL>::~Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::mult starting" << std::endl;
  util::Timer timer(classname(), "mult");
  local_->multiply(dx);
  Log::trace() << "Localization<MODEL>::mult done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Localization<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *local_;
  Log::trace() << "Localization<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCALIZATION_H_
