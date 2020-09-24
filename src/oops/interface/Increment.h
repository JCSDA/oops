/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_INCREMENT_H_
#define OOPS_INTERFACE_INCREMENT_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Metadata.h"

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/State.h"
#if !ATLASIFIED
#include "oops/generic/UnstructuredGrid.h"
#endif
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"
#include "oops/util/Timer.h"

namespace oops {

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in an Increment.
 */

// -----------------------------------------------------------------------------

template <typename MODEL>
class Increment : public oops::GeneralizedDepartures,
                  public util::Printable,
                  public util::Serializable,
                  private util::ObjectCounter<Increment<MODEL> > {
  typedef typename MODEL::Increment  Increment_;
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::Increment";}

/// Constructor, destructor
  Increment(const Geometry_ &, const Variables &, const util::DateTime &);
  Increment(const Geometry_ &, const Increment &);
  Increment(const Increment &, const bool copy = true);
  virtual ~Increment();

/// Interfacing
  Increment_ & increment() {return *increment_;}
  const Increment_ & increment() const {return *increment_;}

/// Interactions with State
  void diff(const State_ &, const State_ &);

/// Time
  const util::DateTime validTime() const {return increment_->validTime();}
  void updateTime(const util::Duration & dt) {increment_->updateTime(dt);}

/// Linear algebra operators
  void zero();
  void zero(const util::DateTime &);
  void ones();
  void dirac(const eckit::Configuration &);
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &, const Increment &, const bool check = true);
  double dot_product_with(const Increment &) const;
  void schur_product_with(const Increment &);
  void random();
  void accumul(const double &, const State_ &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

  LocalIncrement getLocal(const GeometryIterator_ & iter) const;
  void setLocal(const LocalIncrement & gp, const GeometryIterator_ & iter);

/// Get geometry
  Geometry_ geometry() const;
  const Variables & variables() const {return variables_;}

#if ATLASIFIED
/// ATLAS FieldSet
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);
#else
/// Unstructured grid
  void ug_coord(UnstructuredGrid &) const;
  void field_to_ug(UnstructuredGrid &, const int & = 0) const;
  void field_from_ug(const UnstructuredGrid &, const int & = 0);
#endif

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<Increment_> increment_;
  const Variables variables_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
State<MODEL> & operator+=(State<MODEL> & xx, const Increment<MODEL> & dx) {
  Log::trace() << "operator+=(State, Increment) starting" << std::endl;
  util::Timer timer("oops::Increment", "operator+=(State, Increment)");
  xx.state() += dx.increment();
  Log::trace() << "operator+=(State, Increment) done" << std::endl;
  return xx;
}

// =============================================================================
/// Constructor, destructor
// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & resol, const Variables & vars,
                            const util::DateTime & time) : increment_(), variables_(vars)
{
  Log::trace() << "Increment<MODEL>::Increment starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(resol.geometry(), vars, time));
  Log::trace() << "Increment<MODEL>::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & resol, const Increment & other)
  : increment_(), variables_(other.variables_)
{
  Log::trace() << "Increment<MODEL>::Increment starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(resol.geometry(), *other.increment_));
  Log::trace() << "Increment<MODEL>::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Increment & other, const bool copy)
  : increment_(), variables_(other.variables_)
{
  Log::trace() << "Increment<MODEL>::Increment copy starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(*other.increment_, copy));
  Log::trace() << "Increment<MODEL>::Increment copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::~Increment() {
  Log::trace() << "Increment<MODEL>::~Increment starting" << std::endl;
  util::Timer timer(classname(), "~Increment");
  increment_.reset();
  Log::trace() << "Increment<MODEL>::~Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::diff(const State_ & x1, const State_ & x2) {
  Log::trace() << "Increment<MODEL>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  increment_->diff(x1.state(), x2.state());
  Log::trace() << "Increment<MODEL>::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero() {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  increment_->zero();
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero(const util::DateTime & tt) {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  increment_->zero(tt);
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::ones() {
  Log::trace() << "Increment<MODEL>::ones starting" << std::endl;
  util::Timer timer(classname(), "ones");
  increment_->ones();
  Log::trace() << "Increment<MODEL>::ones done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::dirac(const eckit::Configuration & config) {
  Log::trace() << "Increment<MODEL>::dirac starting" << std::endl;
  util::Timer timer(classname(), "dirac");
  increment_->dirac(config);
  Log::trace() << "Increment<MODEL>::dirac done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *increment_ = *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator+=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *increment_ += *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator-=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *increment_ -= *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator*=(const double & zz) {
  Log::trace() << "Increment<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *increment_ *= zz;
  Log::trace() << "Increment<MODEL>::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::axpy(const double & zz, const Increment & dx, const bool check) {
  Log::trace() << "Increment<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  increment_->axpy(zz, *dx.increment_, check);
  Log::trace() << "Increment<MODEL>::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::dot_product_with(const Increment & dx) const {
  Log::trace() << "Increment<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = increment_->dot_product_with(*dx.increment_);
  Log::trace() << "Increment<MODEL>::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::schur_product_with(const Increment & dx) {
  Log::trace() << "Increment<MODEL>::schur_product_with starting" << std::endl;
  util::Timer timer(classname(), "schur_product_with");
  increment_->schur_product_with(*dx.increment_);
  Log::trace() << "Increment<MODEL>::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::random() {
  Log::trace() << "Increment<MODEL>::random starting" << std::endl;
  util::Timer timer(classname(), "random");
  increment_->random();
  Log::trace() << "Increment<MODEL>::random done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::accumul(const double & zz, const State_ & xx) {
  Log::trace() << "Increment<MODEL>::accumul starting" << std::endl;
  util::Timer timer(classname(), "accumul");
  increment_->accumul(zz, xx.state());
  Log::trace() << "Increment<MODEL>::accumul done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalIncrement Increment<MODEL>::getLocal(const GeometryIterator_ & iter) const {
  Log::trace() << "Increment<MODEL>::getLocal starting" << std::endl;
  util::Timer timer(classname(), "getLocal");
  LocalIncrement gp = increment_->getLocal(iter.geometryiter());
  Log::trace() << "Increment<MODEL>::getLocal done" << std::endl;
  return gp;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment<MODEL>::setLocal(const LocalIncrement & gp,
                                const GeometryIterator_ & iter) {
  Log::trace() << "Increment<MODEL>::setLocal starting" << std::endl;
  util::Timer timer(classname(), "setLocal");
  increment_->setLocal(gp, iter.geometryiter());
  Log::trace() << "Increment<MODEL>::setLocal done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "Increment<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  increment_->read(conf);
  Log::trace() << "Increment<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "Increment<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  increment_->write(conf);
  Log::trace() << "Increment<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::norm() const {
  Log::trace() << "Increment<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = increment_->norm();
  Log::trace() << "Increment<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Geometry<MODEL> Increment<MODEL>::geometry() const {
  Log::trace() << "Increment<MODEL>::geometry starting" << std::endl;
  util::Timer timer(classname(), "geometry");
  Geometry<MODEL> geom(increment_->geometry());
  Log::trace() << "Increment<MODEL>::geometry done" << std::endl;
  return geom;
}
// -----------------------------------------------------------------------------
#if ATLASIFIED
// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::setAtlas(atlas::FieldSet * afieldset) const {
  Log::trace() << "Increment<MODEL>::setAtlas starting" << std::endl;
  util::Timer timer(classname(), "setAtlas");
  increment_->setAtlas(afieldset);
  Log::trace() << "Increment<MODEL>::setAtlas done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::toAtlas(atlas::FieldSet * afieldset) const {
  Log::trace() << "Increment<MODEL>::toAtlas starting" << std::endl;
  util::Timer timer(classname(), "toAtlas");
  increment_->toAtlas(afieldset);
  Log::trace() << "Increment<MODEL>::toAtlas done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::fromAtlas(atlas::FieldSet * afieldset) {
  Log::trace() << "Increment<MODEL>::fromAtlas starting" << std::endl;
  util::Timer timer(classname(), "fromAtlas");
  increment_->fromAtlas(afieldset);
  Log::trace() << "Increment<MODEL>::fromAtlas done" << std::endl;
}
// -----------------------------------------------------------------------------
#else
// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::ug_coord(UnstructuredGrid & ug) const {
  Log::trace() << "Increment<MODEL>::ug_coord starting" << std::endl;
  util::Timer timer(classname(), "ug_coord");
  increment_->ug_coord(ug);
  Log::trace() << "Increment<MODEL>::ug_coord done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::field_to_ug(UnstructuredGrid & ug, const int & its) const {
  Log::trace() << "Increment<MODEL>::field_to_ug starting" << std::endl;
  util::Timer timer(classname(), "field_to_ug");
  increment_->field_to_ug(ug, its);
  Log::trace() << "Increment<MODEL>::field_to_ug done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::field_from_ug(const UnstructuredGrid & ug, const int & its) {
  Log::trace() << "Increment<MODEL>::field_from_ug starting" << std::endl;
  util::Timer timer(classname(), "field_from_ug");
  increment_->field_from_ug(ug, its);
  Log::trace() << "Increment<MODEL>::field_from_ug done" << std::endl;
}

// -----------------------------------------------------------------------------
#endif
// -----------------------------------------------------------------------------

template<typename MODEL>
size_t Increment<MODEL>::serialSize() const {
  Log::trace() << "Increment<MODEL>::serialSize" << std::endl;
  util::Timer timer(classname(), "serialSize");
  return increment_->serialSize();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::serialize(std::vector<double> & vect) const {
  Log::trace() << "Increment<MODEL>::serialize starting" << std::endl;
  util::Timer timer(classname(), "serialize");
  increment_->serialize(vect);
  Log::trace() << "Increment<MODEL>::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::deserialize(const std::vector<double> & vect, size_t & current) {
  Log::trace() << "Increment<MODEL>::Increment deserialize starting" << std::endl;
  util::Timer timer(classname(), "deserialize");
  increment_->deserialize(vect, current);
  Log::trace() << "Increment<MODEL>::Increment deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Increment<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *increment_;
  Log::trace() << "Increment<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_INCREMENT_H_
