/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_VARIABLES_H_
#define OOPS_INTERFACE_VARIABLES_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "util/Logger.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class Variables : public util::Printable,
                  private util::ObjectCounter<Variables<MODEL> > {
  typedef typename MODEL::Variables              Variables_;

 public:
  static const std::string classname() {return "oops::Variables";}

  explicit Variables(const eckit::Configuration &);
  explicit Variables(boost::shared_ptr<const Variables_>);
  Variables(const Variables &);
  ~Variables();

/// Interfacing
  const Variables_ & variables() const {return *variables_;}

 private:
  Variables & operator=(const Variables &);
  void print(std::ostream &) const;
  boost::shared_ptr<const Variables_> variables_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL>::Variables(const eckit::Configuration & conf): variables_() {
  Log::trace() << "Variables<MODEL>::Variables starting" << std::endl;
  util::Timer timer(classname(), "Variables");
  variables_.reset(new Variables_(conf));
  Log::trace() << "Variables<MODEL>::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL>::Variables(boost::shared_ptr<const Variables_> ptr): variables_(ptr) {
  Log::trace() << "Variables<MODEL>::Variables shared_ptr done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL>::Variables(const Variables & other): variables_(other.variables_) {
  Log::trace() << "Variables<MODEL>::Variables copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL>::~Variables() {
  Log::trace() << "Variables<MODEL>::~Variables starting" << std::endl;
  util::Timer timer(classname(), "~Variables");
  variables_.reset();
  Log::trace() << "Variables<MODEL>::~Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Variables<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Variables<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *variables_;
  Log::trace() << "Variables<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_VARIABLES_H_
