/*
 * (C) Copyright 2023- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Increment.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Coupled model increment
template <typename MODEL1, typename MODEL2>
class IncrementCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;
  typedef StateCoupled<MODEL1, MODEL2>     StateCoupled_;

 public:
  /// Constructor, destructor
  IncrementCoupled(const GeometryCoupled_ &, const Variables &, const util::DateTime &);
  IncrementCoupled(const GeometryCoupled_ &, const IncrementCoupled &, const bool ad = false);
  IncrementCoupled(const IncrementCoupled &, const bool copy = true);
  virtual ~IncrementCoupled();

  /// Make this increment a difference of two states
  void diff(const StateCoupled_ &, const StateCoupled_ &);

  /// Accessors to the coupled increment components
  Increment<MODEL1> & increment1() {ASSERT(dx1_); return *dx1_;}
  Increment<MODEL2> & increment2() {ASSERT(dx2_); return *dx2_;}
  const Increment<MODEL1> & increment1() const {ASSERT(dx1_); return *dx1_;}
  const Increment<MODEL2> & increment2() const {ASSERT(dx2_); return *dx2_;}

  /// ATLAS
  void toFieldSet(atlas::FieldSet &) const {throw eckit::NotImplemented(Here());}
  void fromFieldSet(const atlas::FieldSet &) {throw eckit::NotImplemented(Here());}

  /// Accessor to the increment's valid time
  const util::DateTime validTime() const;
  /// Update increment's time (increase by \p dt), used in the PseudoModel
  void updateTime(const util::Duration & dt);

  /// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  void dirac(const eckit::Configuration &);
  void ones();
  double norm() const;

  /// Randomize both components of the increment (assuming they are
  /// uncorrelated).
  void random();

  /// Linear algebra operations
  IncrementCoupled & operator=(const IncrementCoupled &);
  IncrementCoupled & operator+=(const IncrementCoupled &);
  IncrementCoupled & operator-=(const IncrementCoupled &);
  IncrementCoupled & operator*=(const double &);
  double dot_product_with(const IncrementCoupled & other) const;
  void schur_product_with(const IncrementCoupled & other);
  void axpy(const double & w, const IncrementCoupled & dx, const bool check = true);

  /// Returns this increment's geometry
  std::shared_ptr<const GeometryCoupled_> geometry() const {return geom_;}
  /// Returns this increment variables
  const Variables & variables() const {return vars_;}

  // For accumulator
  void zero();
  void zero(const util::DateTime &);
  void accumul(const double & zz, const StateCoupled_ & other);

  /// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> & buf) const;
  void deserialize(const std::vector<double> & buf, size_t & ii);

 private:
  static std::string name1() {return MODEL1::name();}
  static std::string name2() {return MODEL2::name();}
  void print(std::ostream &) const;
  std::shared_ptr<const GeometryCoupled_> geom_;
  std::unique_ptr<Increment<MODEL1>> dx1_;
  std::unique_ptr<Increment<MODEL2>> dx2_;
  bool parallel_;
  Variables vars_;
};

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2>::IncrementCoupled(const GeometryCoupled_ & resol,
                                                   const Variables & vars,
                                                   const util::DateTime & time)
  : geom_(new GeometryCoupled_(resol)), dx1_(), dx2_(), parallel_(resol.isParallel()) {
  Log::trace() << "IncrementCoupled::IncrementCoupled starting" << std::endl;
  if (parallel_) {
    if (resol.modelNumber() == 1) {
      dx1_ = std::make_unique<Increment<MODEL1>>(resol.geometry1(), vars, time);
      vars_ = dx1_->variables();
    }
    if (resol.modelNumber() == 2) {
      dx2_ = std::make_unique<Increment<MODEL2>>(resol.geometry2(), vars, time);
      vars_ = dx2_->variables();
    }
  } else {
    dx1_ = std::make_unique<Increment<MODEL1>>(resol.geometry1(), vars, time);
    dx2_ = std::make_unique<Increment<MODEL2>>(resol.geometry2(), vars, time);
    vars_ = dx1_->variables();
    vars_ += dx2_->variables();
  }
  Log::trace() << "IncrementCoupled::IncrementCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2>::IncrementCoupled(const GeometryCoupled_ & resol,
                                                   const IncrementCoupled & other,
                                                   const bool ad)
  : geom_(new GeometryCoupled_(resol)), dx1_(), dx2_(), parallel_(resol.isParallel()),
    vars_(other.vars_) {
  Log::trace() << "IncrementCoupled::IncrementCoupled interpolated starting" << std::endl;
  ASSERT(parallel_ == other.parallel_);
  if (other.dx1_) dx1_ = std::make_unique<Increment<MODEL1>>(resol.geometry1(), *other.dx1_, ad);
  if (other.dx2_) dx2_ = std::make_unique<Increment<MODEL2>>(resol.geometry2(), *other.dx2_, ad);
  Log::trace() << "IncrementCoupled::IncrementCoupled interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2>::IncrementCoupled(const IncrementCoupled & other, const bool copy)
  : geom_(other.geom_), dx1_(), dx2_(), parallel_(other.parallel_), vars_(other.vars_) {
  Log::trace() << "IncrementCoupled::IncrementCoupled copy starting" << std::endl;
  if (other.dx1_) dx1_ = std::make_unique<Increment<MODEL1>>(*other.dx1_, copy);
  if (other.dx2_) dx2_ = std::make_unique<Increment<MODEL2>>(*other.dx2_, copy);
  Log::trace() << "IncrementCoupled::IncrementCoupled copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2>::~IncrementCoupled() {
  Log::trace() << "IncrementCoupled::~IncrementCoupled starting" << std::endl;
  dx1_.reset();
  dx2_.reset();
  Log::trace() << "IncrementCoupled::~IncrementCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::diff(const StateCoupled_ & xx1,
                                            const StateCoupled_ & xx2) {
  Log::trace() << "IncrementCoupled::diff starting" << std::endl;
  dx1_->diff(xx1.state1(), xx2.state1());
  dx2_->diff(xx1.state2(), xx2.state2());
  Log::trace() << "IncrementCoupled::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::read(const eckit::Configuration & config) {
  Log::trace() << "IncrementCoupled::read starting" << std::endl;
  if (dx1_) {
    const eckit::LocalConfiguration conf1(config, name1());
    dx1_->read(conf1);
  }
  if (dx2_) {
    const eckit::LocalConfiguration conf2(config, name2());
    dx2_->read(conf2);
  }
  Log::trace() << "IncrementCoupled::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::write(const eckit::Configuration & config) const {
  Log::trace() << "IncrementCoupled::write starting" << std::endl;
  if (dx1_) {
    const eckit::LocalConfiguration conf1(config, name1());
    dx1_->write(conf1);
  }
  if (dx2_) {
    const eckit::LocalConfiguration conf2(config, name2());
    dx2_->write(conf2);
  }
  Log::trace() << "IncrementCoupled::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::dirac(const eckit::Configuration & config) {
  Log::trace() << "IncrementCoupled::dirac starting" << std::endl;
  if (dx1_) {
    const eckit::LocalConfiguration conf1(config, name1());
    dx1_->dirac(conf1);
  }
  if (dx2_) {
    const eckit::LocalConfiguration conf2(config, name2());
    dx2_->dirac(conf2);
  }
  Log::trace() << "IncrementCoupled::dirac done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
double IncrementCoupled<MODEL1, MODEL2>::norm() const {
  Log::trace() << "IncrementCoupled::norm starting" << std::endl;
  double zz = 0.0;
  if (parallel_) {
    if (dx1_) zz = dx1_->norm();
    if (dx2_) zz = dx2_->norm();
    geom_->getCommPairRanks().allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  } else {
    zz = dx1_->norm() + dx2_->norm();
  }
  Log::trace() << "IncrementCoupled::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::random() {
  Log::trace() << "IncrementCoupled::random starting" << std::endl;
  if (dx1_) dx1_->random();
  if (dx2_) dx2_->random();
  Log::trace() << "IncrementCoupled::random done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
double IncrementCoupled<MODEL1, MODEL2>::dot_product_with(const IncrementCoupled & other) const {
  Log::trace() << "IncrementCoupled::dot_product_with starting" << std::endl;
  double zz = 0.0;
  if (parallel_) {
    if (dx1_) zz = dx1_->dot_product_with(other.increment1());
    if (dx2_) zz = dx2_->dot_product_with(other.increment2());
    geom_->getCommPairRanks().allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  } else {
    zz = dx1_->dot_product_with(other.increment1()) +
         dx2_->dot_product_with(other.increment2());
  }
  Log::trace() << "IncrementCoupled::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::schur_product_with(const IncrementCoupled & other) {
  Log::trace() << "IncrementCoupled::schur_product_with starting" << std::endl;
  if (dx1_) dx1_->schur_product_with(other.increment1());
  if (dx2_) dx2_->schur_product_with(other.increment2());
  Log::trace() << "IncrementCoupled::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::axpy(const double & w,
                     const IncrementCoupled & dx, const bool check) {
  if (dx1_) dx1_->axpy(w, *dx.dx1_, check);
  if (dx2_) dx2_->axpy(w, *dx.dx2_, check);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
const util::DateTime IncrementCoupled<MODEL1, MODEL2>::validTime() const {
  util::DateTime dt;
  if (dx1_) dt = dx1_->validTime();
  if (dx2_) dt = dx2_->validTime();
  return dt;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::updateTime(const util::Duration & dt) {
  if (dx1_) dx1_->updateTime(dt);
  if (dx2_) dx2_->updateTime(dt);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::zero() {
  if (dx1_) dx1_->zero();
  if (dx2_) dx2_->zero();
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::zero(const util::DateTime & time) {
  if (dx1_) dx1_->zero(time);
  if (dx2_) dx2_->zero(time);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::ones() {
  if (dx1_) dx1_->ones();
  if (dx2_) dx2_->ones();
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::accumul(const double & zz, const StateCoupled_ & other) {
  if (dx1_) dx1_->accumul(zz, other.state1());
  if (dx2_) dx2_->accumul(zz, other.state2());
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2> &
IncrementCoupled<MODEL1, MODEL2>::operator=(const IncrementCoupled & other) {
  if (dx1_) *dx1_ = *other.dx1_;
  if (dx2_) *dx2_ = *other.dx2_;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2> &
IncrementCoupled<MODEL1, MODEL2>::operator+=(const IncrementCoupled & other) {
  if (dx1_) *dx1_ += *other.dx1_;
  if (dx2_) *dx2_ += *other.dx2_;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2> &
IncrementCoupled<MODEL1, MODEL2>::operator-=(const IncrementCoupled & other) {
  if (dx1_) *dx1_ -= *other.dx1_;
  if (dx2_) *dx2_ -= *other.dx2_;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
IncrementCoupled<MODEL1, MODEL2> & IncrementCoupled<MODEL1, MODEL2>::operator*=(const double & k) {
  if (dx1_) *dx1_ *= k;
  if (dx2_) *dx2_ *= k;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
size_t IncrementCoupled<MODEL1, MODEL2>::serialSize() const {
  size_t ss = 0;
  if (dx1_) ss += dx1_->serialSize();
  if (dx2_) ss += dx2_->serialSize();
  return ss;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::serialize(std::vector<double> & buf) const {
  if (dx1_) dx1_->serialize(buf);
  if (dx2_) dx2_->serialize(buf);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::deserialize(const std::vector<double> & buf, size_t & ii) {
  if (dx1_) dx1_->deserialize(buf, ii);
  if (dx2_) dx2_->deserialize(buf, ii);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void IncrementCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "IncrementCoupled::print starting" << std::endl;

  if (parallel_) {
    std::stringstream ss;
    ss.setf(os.flags());
    ss.precision(os.precision());
    if (dx1_) {
      ss << std::endl << "IncrementCoupled: " << MODEL1::name() << std::endl;
      ss << *dx1_ << std::endl;
    }
    if (dx2_) {
      ss << std::endl << "IncrementCoupled: " << MODEL2::name() << std::endl;
      ss << *dx2_ << std::endl;
    }
    util::gatherPrint(os, ss.str(), geom_->getCommPairRanks());
  } else {
    os << std::endl << "IncrementCoupled: " << MODEL1::name() << std::endl;
    os << *dx1_ << std::endl;
    os << std::endl << "IncrementCoupled: " << MODEL2::name() << std::endl;
    os << *dx2_;
  }

  Log::trace() << "IncrementCoupled::print done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2> & operator+=(StateCoupled<MODEL1, MODEL2> & xx,
                                          const IncrementCoupled<MODEL1, MODEL2> & dx) {
  if (!xx.geometry()->isParallel() || xx.geometry()->modelNumber() == 1) {
    xx.state1() += dx.increment1();
  }
  if (!xx.geometry()->isParallel() || xx.geometry()->modelNumber() == 2) {
    xx.state2() += dx.increment2();
  }
  return xx;
}

}  // namespace oops
