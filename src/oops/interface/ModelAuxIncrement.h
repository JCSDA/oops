/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_MODELAUXINCREMENT_H_
#define OOPS_INTERFACE_MODELAUXINCREMENT_H_

#include <iostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ModelAuxIncrement : public util::Printable,
                          private util::ObjectCounter<ModelAuxIncrement<MODEL> > {
  typedef typename MODEL::ModelAuxIncrement     ModelAuxIncrement_;
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxControl<MODEL>      ModelAuxControl_;

 public:
  static const std::string classname() {return "oops::ModelAuxIncrement";}

/// Constructor, destructor
  ModelAuxIncrement(const Geometry_ &, const eckit::Configuration &);
  explicit ModelAuxIncrement(const ModelAuxIncrement &, const bool copy = true);
  ModelAuxIncrement(const ModelAuxIncrement &, const eckit::Configuration &);
  ~ModelAuxIncrement();

/// Interfacing
  const ModelAuxIncrement_ & modelauxincrement() const {return *aux_;}
  ModelAuxIncrement_ & modelauxincrement() {return *aux_;}

/// Linear algebra operators
  void diff(const ModelAuxControl_ &, const ModelAuxControl_ &);
  void zero();
  ModelAuxIncrement & operator=(const ModelAuxIncrement &);
  ModelAuxIncrement & operator+=(const ModelAuxIncrement &);
  ModelAuxIncrement & operator-=(const ModelAuxIncrement &);
  ModelAuxIncrement & operator*=(const double &);
  void axpy(const double &, const ModelAuxIncrement &);
  double dot_product_with(const ModelAuxIncrement &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<ModelAuxIncrement_> aux_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelAuxControl<MODEL> & operator+=(ModelAuxControl<MODEL> & xx,
                                    const ModelAuxIncrement<MODEL> & dx) {
  Log::trace() << "operator+=(ModelAuxControl, ModelAuxIncrement) starting" << std::endl;
  util::Timer timer("oops::ModelAuxIncrement", "operator+=ModelAuxControl");
  xx.modelauxcontrol() += dx.modelauxincrement();
  Log::trace() << "operator+=(ModelAuxControl, ModelAuxIncrement) done" << std::endl;
  return xx;
}

// =============================================================================

template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const Geometry_ & resol,
                                            const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrement");
  aux_.reset(new ModelAuxIncrement_(resol.geometry(), conf));
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const ModelAuxIncrement & other,
                                            const bool copy) : aux_()
{
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement copy starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrement");
  aux_.reset(new ModelAuxIncrement_(*other.aux_, copy));
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement copy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const ModelAuxIncrement & other,
                                            const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement interpolated starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrement");
  aux_.reset(new ModelAuxIncrement_(*other.aux_, conf));
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement interpolated done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::~ModelAuxIncrement() {
  Log::trace() << "ModelAuxIncrement<MODEL>::~ModelAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "~ModelAuxIncrement");
  aux_.reset();
  Log::trace() << "ModelAuxIncrement<MODEL>::~ModelAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::diff(const ModelAuxControl_ & x1, const ModelAuxControl_ & x2) {
  Log::trace() << "ModelAuxIncrement<MODEL>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  aux_->diff(x1.modelauxcontrol(), x2.modelauxcontrol());
  Log::trace() << "ModelAuxIncrement<MODEL>::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::zero() {
  Log::trace() << "ModelAuxIncrement<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  aux_->zero();
  Log::trace() << "ModelAuxIncrement<MODEL>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator=(const ModelAuxIncrement & rhs) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *aux_ = *rhs.aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator+=(const ModelAuxIncrement & rhs) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *aux_ += *rhs.aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator-=(const ModelAuxIncrement & rhs) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *aux_ -= *rhs.aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator*=(const double & zz) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *aux_ *= zz;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::axpy(const double & zz, const ModelAuxIncrement & dx) {
  Log::trace() << "ModelAuxIncrement<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  aux_->axpy(zz, *dx.aux_);
  Log::trace() << "ModelAuxIncrement<MODEL>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ModelAuxIncrement<MODEL>::dot_product_with(const ModelAuxIncrement & dx) const {
  Log::trace() << "ModelAuxIncrement<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = aux_->dot_product_with(*dx.aux_);
  Log::trace() << "ModelAuxIncrement<MODEL>::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ModelAuxIncrement<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(conf);
  Log::trace() << "ModelAuxIncrement<MODEL>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ModelAuxIncrement<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(conf);
  Log::trace() << "ModelAuxIncrement<MODEL>::write done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ModelAuxIncrement<MODEL>::norm() const {
  Log::trace() << "ModelAuxIncrement<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ModelAuxIncrement<MODEL>::norm done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ModelAuxIncrement<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELAUXINCREMENT_H_
