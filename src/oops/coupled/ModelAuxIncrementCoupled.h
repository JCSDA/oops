/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/coupled/GeometryCoupled.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace oops {

/// Coupled implementation of ModelAuxIncrement
template <typename MODEL1, typename MODEL2>
class ModelAuxIncrementCoupled : public util::Printable,
                                 public util::Serializable {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;
  typedef AuxCoupledModel<MODEL1, MODEL2>  ModelAuxControlCoupled_;
  typedef ModelAuxIncrement<MODEL1>        ModelAuxIncrement1;
  typedef ModelAuxIncrement<MODEL2>        ModelAuxIncrement2;

 public:
  static const std::string classname() {return "oops::ModelAuxIncrementCoupled";}

  /// Constructor for specified \p resol and \p conf
  ModelAuxIncrementCoupled(const GeometryCoupled_ & resol, const eckit::Configuration & conf);
  /// Copies \p other ModelAuxIncrementCoupled if \p copy is true,
  /// otherwise creates zero ModelAuxIncrementCoupled with same variables and geometry
  explicit ModelAuxIncrementCoupled(const ModelAuxIncrementCoupled & other, const bool copy = true);
  /// Copies \p other ModelAuxIncrementCoupled, reading extra information from \p conf
  ModelAuxIncrementCoupled(const ModelAuxIncrementCoupled & other,
      const eckit::Configuration & conf);
  /// Destructor (defined explicitly for timing and tracing)
  ~ModelAuxIncrementCoupled() = default;

  /// Sets this ModelAuxIncrement to the difference between two ModelAuxControlCoupled objects
  void diff(const ModelAuxControlCoupled_ &, const ModelAuxControlCoupled_ &);
  /// Zero out this ModelAuxIncrement
  void zero();
  /// Linear algebra operators
  ModelAuxIncrementCoupled & operator=(const ModelAuxIncrementCoupled &);
  ModelAuxIncrementCoupled & operator+=(const ModelAuxIncrementCoupled &);
  ModelAuxIncrementCoupled & operator-=(const ModelAuxIncrementCoupled &);
  ModelAuxIncrementCoupled & operator*=(const double &);
  void axpy(const double &, const ModelAuxIncrementCoupled &);
  /// dot product with the \p other ModelAuxIncrementCoupled
  double dot_product_with(const ModelAuxIncrementCoupled & other) const;

  /// Read this ModelAuxIncrement from file
  void read(const eckit::Configuration &);
  /// Write this ModelAuxIncrement out to file
  void write(const eckit::Configuration &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Serialize and deserialize (used in 4DEnVar, weak-constraint 4DVar and Block-Lanczos minimizer)
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  /// Accessors to the coupled  components
  ModelAuxIncrement1 & modelauxincrement1() {ASSERT(inc1_); return *inc1_;}
  ModelAuxIncrement2 & modelauxincrement2() {ASSERT(inc2_); return *inc2_;}
  const ModelAuxIncrement1 & modelauxincrement1() const {ASSERT(inc1_); return *inc1_;}
  const ModelAuxIncrement2 & modelauxincrement2() const {ASSERT(inc2_); return *inc2_;}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ModelAuxIncrement1> inc1_;
  std::unique_ptr<ModelAuxIncrement2> inc2_;
};

// -----------------------------------------------------------------------------
template <typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2> & operator+=(AuxCoupledModel<MODEL1, MODEL2> & xx,
                                    const ModelAuxIncrementCoupled<MODEL1, MODEL2> & dx) {
  Log::trace() << "operator+=(ModelAuxControl, ModelAuxIncrement) starting" << std::endl;
  xx.aux1() += dx.modelauxincrement1();
  xx.aux2() += dx.modelauxincrement2();
  Log::trace() << "operator+=(ModelAuxControl, ModelAuxIncrement) done" << std::endl;
  return xx;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxIncrementCoupled<MODEL1, MODEL2>::ModelAuxIncrementCoupled(
    const GeometryCoupled_ & resol,
    const eckit::Configuration & conf) : inc1_(), inc2_() {
  Log::trace() << "ModelAuxIncrementCoupled::ModelAuxIncrementCoupled starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrementCoupled");
  inc1_.reset(new ModelAuxIncrement1(resol.geometry1(), conf));
  inc2_.reset(new ModelAuxIncrement2(resol.geometry2(), conf));
  Log::trace() << "ModelAuxIncrementCoupled::ModelAuxIncrementCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxIncrementCoupled<MODEL1, MODEL2>::ModelAuxIncrementCoupled(
    const ModelAuxIncrementCoupled & other, const bool copy) : inc1_(), inc2_() {
  Log::trace() << "ModelAuxIncrementCoupled::ModelAuxIncrementCoupled copy starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrementCoupled");
  inc1_.reset(new ModelAuxIncrement1(*other.inc1_, copy));
  inc2_.reset(new ModelAuxIncrement2(*other.inc2_, copy));
  Log::trace() << "ModelAuxIncrementCoupled::ModelAuxIncrementCoupled copy done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxIncrementCoupled<MODEL1, MODEL2>::ModelAuxIncrementCoupled(
    const ModelAuxIncrementCoupled & other, const eckit::Configuration & conf) : inc1_(), inc2_() {
  Log::trace() << "ModelAuxIncrementCoupled::ModelAuxIncrementCoupled conf starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrementCoupled");
  inc1_.reset(new ModelAuxIncrement1(*other.inc1_, conf));
  inc2_.reset(new ModelAuxIncrement2(*other.inc2_, conf));
  Log::trace() << "ModelAuxIncrementCoupled::ModelAuxIncrementCoupled conf done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::diff(
    const ModelAuxControlCoupled_ & xx1, const ModelAuxControlCoupled_ & xx2) {
  Log::trace() << "ModelAuxIncrementCoupled::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  inc1_->diff(xx1.aux1(), xx2.aux1());
  inc2_->diff(xx1.aux2(), xx2.aux2());
  Log::trace() << "ModelAuxIncrementCoupled::diff done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::zero() {
  Log::trace() << "ModelAuxIncrementCoupled::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  inc1_->zero();
  inc2_->zero();
  Log::trace() << "ModelAuxIncrementCoupled::zero done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxIncrementCoupled<MODEL1, MODEL2> & ModelAuxIncrementCoupled<MODEL1, MODEL2>::operator=(
    const ModelAuxIncrementCoupled & rhs) {
  Log::trace() << "ModelAuxIncrementCoupled::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *inc1_ = *rhs.inc1_;
  *inc2_ = *rhs.inc2_;
  Log::trace() << "ModelAuxIncrementCoupled::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxIncrementCoupled<MODEL1, MODEL2> & ModelAuxIncrementCoupled<MODEL1, MODEL2>::operator+=(
    const ModelAuxIncrementCoupled & rhs) {
  Log::trace() << "ModelAuxIncrementCoupled::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *inc1_ += *rhs.inc1_;
  *inc2_ += *rhs.inc2_;
  Log::trace() << "ModelAuxIncrementCoupled::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxIncrementCoupled<MODEL1, MODEL2> & ModelAuxIncrementCoupled<MODEL1, MODEL2>::operator-=(
    const ModelAuxIncrementCoupled & rhs) {
  Log::trace() << "ModelAuxIncrementCoupled::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *inc1_ -= *rhs.inc1_;
  *inc2_ -= *rhs.inc2_;
  Log::trace() << "ModelAuxIncrementCoupled::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxIncrementCoupled<MODEL1, MODEL2> & ModelAuxIncrementCoupled<MODEL1, MODEL2>::operator*=(
    const double & zz) {
  Log::trace() << "ModelAuxIncrementCoupled::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *inc1_ *= zz;
  *inc2_ *= zz;
  Log::trace() << "ModelAuxIncrementCoupled::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::axpy(
    const double & zz, const ModelAuxIncrementCoupled & dx) {
  Log::trace() << "ModelAuxIncrementCoupled::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  inc1_->axpy(zz, *dx.inc1_);
  inc2_->axpy(zz, *dx.inc2_);
  Log::trace() << "ModelAuxIncrementCoupled::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
double ModelAuxIncrementCoupled<MODEL1, MODEL2>::dot_product_with(
    const ModelAuxIncrementCoupled & dx) const {
  Log::trace() << "ModelAuxIncrementCoupled::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = 0.0;
  if (inc1_) zz += inc1_->dot_product_with(*dx.inc1_);
  if (inc2_) zz += inc2_->dot_product_with(*dx.inc2_);
  Log::trace() << "ModelAuxIncrementCoupled::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::read(
    const eckit::Configuration & conf) {
  Log::trace() << "ModelAuxIncrementCoupled::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  inc1_->read(conf);
  inc2_->read(conf);
  Log::trace() << "ModelAuxIncrementCoupled::read done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::write(
    const eckit::Configuration & conf) const {
  Log::trace() << "ModelAuxIncrementCoupled::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  inc1_->write(conf);
  inc2_->write(conf);
  Log::trace() << "ModelAuxIncrementCoupled::write done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
double ModelAuxIncrementCoupled<MODEL1, MODEL2>::norm() const {
  Log::trace() << "ModelAuxIncrementCoupled::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = 0.0;
  if (inc1_) zz += inc1_->norm();
  if (inc2_) zz += inc2_->norm();
  Log::trace() << "ModelAuxIncrementCoupled::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
size_t ModelAuxIncrementCoupled<MODEL1, MODEL2>::serialSize() const {
  Log::trace() << "ModelAuxIncrementCoupled::serialSize starting" << std::endl;
  util::Timer timer(classname(), "serialSize");
  size_t sz = inc1_->serialSize() + inc2_->serialSize();
  Log::trace() << "ModelAuxIncrementCoupled::serialSize done" << std::endl;
  return sz;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::serialize(
    std::vector<double> & vect) const {
  Log::trace() << "ModelAuxIncrementCoupled::serialize starting" << std::endl;
  util::Timer timer(classname(), "serialize");
  inc1_->serialize(vect);
  inc2_->serialize(vect);
  Log::trace() << "ModelAuxIncrementCoupled::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::deserialize(
    const std::vector<double> & vect, size_t & ii) {
  Log::trace() << "ModelAuxIncrementCoupled::deserialize starting" << std::endl;
  util::Timer timer(classname(), "deserialize");
  inc1_->deserialize(vect, ii);
  inc2_->deserialize(vect, ii);
  Log::trace() << "ModelAuxIncrementCoupled::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxIncrementCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "ModelAuxIncrementCoupled::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *inc1_ << *inc2_;
  Log::trace() << "ModelAuxIncrementCoupled::print done" << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace oops
