/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSAUXCONTROL_H_
#define OOPS_INTERFACE_OBSAUXCONTROL_H_

#include <iostream>
#include <string>

#include <boost/scoped_ptr.hpp>

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
class ObsAuxControl : public util::Printable,
                      private util::ObjectCounter<ObsAuxControl<MODEL> > {
  typedef typename MODEL::ObsAuxControl        ObsAuxControl_;

 public:
  static const std::string classname() {return "oops::ObsAuxControl";}

  explicit ObsAuxControl(const eckit::Configuration &);
  explicit ObsAuxControl(const ObsAuxControl &, const bool copy = true);
  ~ObsAuxControl();

/// Interfacing
  const ObsAuxControl_ & obsauxcontrol() const {return *aux_;}
  ObsAuxControl_ & obsauxcontrol() {return *aux_;}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

 private:
  ObsAuxControl & operator=(const ObsAuxControl &);
  void print(std::ostream &) const;
  boost::scoped_ptr<ObsAuxControl_> aux_;
};

// =============================================================================

template<typename MODEL>
ObsAuxControl<MODEL>::ObsAuxControl(const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ObsAuxControl<MODEL>::ObsAuxControl starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxControl");
  aux_.reset(new ObsAuxControl_(conf));
  Log::trace() << "ObsAuxControl<MODEL>::ObsAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsAuxControl<MODEL>::ObsAuxControl(const ObsAuxControl & other, const bool copy) : aux_()
{
  Log::trace() << "ObsAuxControl<MODEL>::ObsAuxControl copy starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxControl");
  aux_.reset(new ObsAuxControl_(*other.aux_, copy));
  Log::trace() << "ObsAuxControl<MODEL>::ObsAuxControl copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsAuxControl<MODEL>::~ObsAuxControl() {
  Log::trace() << "ObsAuxControl<MODEL>::~ObsAuxControl starting" << std::endl;
  util::Timer timer(classname(), "~ObsAuxControl");
  aux_.reset();
  Log::trace() << "ObsAuxControl<MODEL>::~ObsAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxControl<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ObsAuxControl<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(conf);
  Log::trace() << "ObsAuxControl<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxControl<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ObsAuxControl<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(conf);
  Log::trace() << "ObsAuxControl<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double ObsAuxControl<MODEL>::norm() const {
  Log::trace() << "ObsAuxControl<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ObsAuxControl<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxControl<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxControl<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ObsAuxControl<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSAUXCONTROL_H_
