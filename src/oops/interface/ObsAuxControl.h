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
#include <memory>
#include <string>

#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxControl : public util::Printable,
                      private util::ObjectCounter<ObsAuxControl<MODEL> > {
  typedef typename MODEL::ObsAuxControl        ObsAuxControl_;

 public:
  static const std::string classname() {return "oops::ObsAuxControl";}

  ObsAuxControl(const ObsSpace<MODEL> &, const eckit::Configuration &);
  explicit ObsAuxControl(const ObsAuxControl &, const bool copy = true);
  ~ObsAuxControl();

/// Interfacing
  const ObsAuxControl_ & obsauxcontrol() const {return *aux_;}
  ObsAuxControl_ & obsauxcontrol() {return *aux_;}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

/// Other
  const Variables & requiredGeoVaLs() const;
  const Variables & requiredHdiagnostics() const;

/// Operator
  ObsAuxControl & operator=(const ObsAuxControl &);

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsAuxControl_> aux_;
};

// =============================================================================

template<typename MODEL>
ObsAuxControl<MODEL>::ObsAuxControl(const ObsSpace<MODEL> & os,
                                    const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ObsAuxControl<MODEL>::ObsAuxControl starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxControl");
  aux_.reset(new ObsAuxControl_(os.obsspace(), conf));
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
const Variables & ObsAuxControl<MODEL>::requiredGeoVaLs() const {
  Log::trace() << "ObsAuxControl<MODEL>::requiredGeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "requiredGeoVaLs");
  Log::trace() << "ObsAuxControl<MODEL>::requiredGeoVaLs done" << std::endl;
  return aux_->requiredGeoVaLs();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const Variables & ObsAuxControl<MODEL>::requiredHdiagnostics() const {
  Log::trace() << "ObsAuxControl<MODEL>::requiredHdiagnostics starting" << std::endl;
  util::Timer timer(classname(), "requiredHdiagnostics");
  Log::trace() << "ObsAuxControl<MODEL>::requiredHdiagnostics done" << std::endl;
  return aux_->requiredHdiagnostics();
}

// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxControl<MODEL> & ObsAuxControl<MODEL>::operator=(const ObsAuxControl & rhs) {
  Log::trace() << "ObsAuxControl<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *aux_ = *rhs.aux_;
  Log::trace() << "ObsAuxControl<MODEL>::operator= done" << std::endl;
  return *this;
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
