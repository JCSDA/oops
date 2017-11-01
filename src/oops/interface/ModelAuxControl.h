/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_MODELAUXCONTROL_H_
#define OOPS_INTERFACE_MODELAUXCONTROL_H_

#include <iostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/Geometry.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ModelAuxControl : public util::Printable,
                        private util::ObjectCounter<ModelAuxControl<MODEL> > {
  typedef typename MODEL::ModelAuxControl      ModelAuxControl_;
  typedef Geometry<MODEL>            Geometry_;

 public:
  static const std::string classname() {return "oops::ModelAuxControl";}

  ModelAuxControl(const Geometry_ &, const eckit::Configuration &);
  ModelAuxControl(const Geometry_ &, const ModelAuxControl &);
  explicit ModelAuxControl(const ModelAuxControl &, const bool copy = true);
  ~ModelAuxControl();

/// Interfacing
  const ModelAuxControl_ & modelauxcontrol() const {return *aux_;}
  ModelAuxControl_ & modelauxcontrol() {return *aux_;}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

 private:
  ModelAuxControl & operator=(const ModelAuxControl &);
  void print(std::ostream &) const;
  boost::scoped_ptr<ModelAuxControl_> aux_;
};

// =============================================================================

template<typename MODEL>
ModelAuxControl<MODEL>::ModelAuxControl(const Geometry_ & resol,
                                        const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxControl");
  aux_.reset(new ModelAuxControl_(resol.geometry(), conf));
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxControl<MODEL>::ModelAuxControl(const Geometry_ & resol,
                                        const ModelAuxControl & other) : aux_()
{
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl interpolated starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxControl");
  aux_.reset(new ModelAuxControl_(resol.geometry(), *other.aux_));
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxControl<MODEL>::ModelAuxControl(const ModelAuxControl & other,
                                        const bool copy) : aux_()
{
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl copy starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxControl");
  aux_.reset(new ModelAuxControl_(*other.aux_, copy));
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxControl<MODEL>::~ModelAuxControl() {
  Log::trace() << "ModelAuxControl<MODEL>::~ModelAuxControl starting" << std::endl;
  util::Timer timer(classname(), "~ModelAuxControl");
  aux_.reset();
  Log::trace() << "ModelAuxControl<MODEL>::~ModelAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxControl<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ModelAuxControl<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(conf);
  Log::trace() << "ModelAuxControl<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxControl<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ModelAuxControl<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(conf);
  Log::trace() << "ModelAuxControl<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double ModelAuxControl<MODEL>::norm() const {
  Log::trace() << "ModelAuxControl<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ModelAuxControl<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxControl<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ModelAuxControl<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ModelAuxControl<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELAUXCONTROL_H_
