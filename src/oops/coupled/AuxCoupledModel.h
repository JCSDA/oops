/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Printable.h"

#include "oops/coupled/GeometryCoupled.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Implementation of the ModelAuxControl interface for a coupled model.
template <typename MODEL1, typename MODEL2>
class AuxCoupledModel : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;

 public:
  AuxCoupledModel(const GeometryCoupled_ &, const eckit::Configuration &);
  AuxCoupledModel(const GeometryCoupled_ &, const AuxCoupledModel &);
  AuxCoupledModel(const AuxCoupledModel &, const bool);
  ~AuxCoupledModel();

  /// Accessors to the individual components of the coupled ModelAuxControl
  ModelAuxControl<MODEL1> & aux1() {ASSERT(aux1_); return *aux1_;}
  ModelAuxControl<MODEL2> & aux2() {ASSERT(aux2_); return *aux2_;}
  const ModelAuxControl<MODEL1> & aux1() const {ASSERT(aux1_); return *aux1_;}
  const ModelAuxControl<MODEL2> & aux2() const {ASSERT(aux2_); return *aux2_;}

  /// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

 private:
  void print(std::ostream &) const override;
  std::shared_ptr<const GeometryCoupled_> geom_;
  std::unique_ptr<ModelAuxControl<MODEL1>> aux1_;
  std::unique_ptr<ModelAuxControl<MODEL2>> aux2_;
  bool parallel_;
};

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2>::AuxCoupledModel(const GeometryCoupled_ & geom,
                                                 const eckit::Configuration & config)
  : geom_(new GeometryCoupled_(geom)), aux1_(), aux2_(), parallel_(geom.isParallel())
{
  Log::trace() << "AuxCoupledModel::AuxCoupledModel starting" << std::endl;
  const eckit::LocalConfiguration conf1 = config.getSubConfiguration(MODEL1::name());;
  const eckit::LocalConfiguration conf2 = config.getSubConfiguration(MODEL2::name());;
  if (parallel_) {
    if (geom.modelNumber() == 1) {
      aux1_ = std::make_unique<ModelAuxControl<MODEL1>>(geom.geometry1(), conf1);
    }
    if (geom.modelNumber() == 2) {
      aux2_ = std::make_unique<ModelAuxControl<MODEL2>>(geom.geometry2(), conf2);
    }
  } else {
    aux1_ = std::make_unique<ModelAuxControl<MODEL1>>(geom.geometry1(), conf1);
    aux2_ = std::make_unique<ModelAuxControl<MODEL2>>(geom.geometry2(), conf2);
  }

  Log::trace() << "AuxCoupledModel::AuxCoupledModel read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2>::AuxCoupledModel(const GeometryCoupled_ & geom,
                                                 const AuxCoupledModel & other)
  : geom_(new GeometryCoupled_(geom)), aux1_(), aux2_(), parallel_(other.parallel_)
{
  Log::trace() << "AuxCoupledModel::AuxCoupledModel interpolated starting" << std::endl;
  if (other.aux1_) aux1_.reset(new ModelAuxControl<MODEL1>(geom.geometry1(), *other.aux1_));
  if (other.aux2_) aux2_.reset(new ModelAuxControl<MODEL2>(geom.geometry2(), *other.aux2_));
  Log::trace() << "AuxCoupledModel::AuxCoupledModel interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2>::AuxCoupledModel(const AuxCoupledModel & other,
                                                 const bool copy)
  : geom_(other.geom_), aux1_(), aux2_(), parallel_(other.parallel_)
{
  Log::trace() << "AuxCoupledModel::AuxCoupledModel copy starting" << std::endl;
  if (other.aux1_) aux1_.reset(new ModelAuxControl<MODEL1>(*other.aux1_, copy));
  if (other.aux2_) aux2_.reset(new ModelAuxControl<MODEL2>(*other.aux2_, copy));
  Log::trace() << "AuxCoupledModel::AuxCoupledModel copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2>::~AuxCoupledModel() {
  Log::trace() << "AuxCoupledModel::~AuxCoupledModel starting" << std::endl;
  aux1_.reset();
  aux2_.reset();
  Log::trace() << "AuxCoupledModel::~AuxCoupledModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void AuxCoupledModel<MODEL1, MODEL2>::read(const eckit::Configuration & conf) {
  Log::trace() << "AuxCoupledModel::read starting" << std::endl;
  if (aux1_) aux1_->read(conf.getSubConfiguration(MODEL1::name()));
  if (aux2_) aux2_->read(conf.getSubConfiguration(MODEL2::name()));
  Log::trace() << "AuxCoupledModel::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void AuxCoupledModel<MODEL1, MODEL2>::write(const eckit::Configuration & conf) const {
  Log::trace() << "AuxCoupledModel::write starting" << std::endl;
  if (aux1_) aux1_->write(conf.getSubConfiguration(MODEL1::name()));
  if (aux2_) aux2_->write(conf.getSubConfiguration(MODEL2::name()));
  Log::trace() << "AuxCoupledModel::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
double AuxCoupledModel<MODEL1, MODEL2>::norm() const {
  Log::trace() << "AuxCoupledModel::norm starting" << std::endl;
  double zz = 0.0;
  if (parallel_) {
    if (aux1_) zz = aux1_->norm();
    if (aux2_) zz = aux2_->norm();
    geom_->getCommPairRanks().allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  } else {
    zz = aux1_->norm() + aux2_->norm();
  }
  Log::trace() << "AuxCoupledModel::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void AuxCoupledModel<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "AuxCoupledModel::print starting" << std::endl;

  if (parallel_) {
    std::stringstream ss;
    ss.setf(os.flags());
    ss.precision(os.precision());
    if (aux1_) {
      ss << std::endl << "AuxCoupledModel: " << MODEL1::name() << std::endl;
      ss << *aux1_ << std::endl;
    }
    if (aux2_) {
      ss << std::endl << "AuxCoupledModel: " << MODEL2::name() << std::endl;
      ss << *aux2_ << std::endl;
    }
    util::gatherPrint(os, ss.str(), geom_->getCommPairRanks());
  } else {
    os << std::endl << "AuxCoupledModel: " << MODEL1::name() << std::endl;
    os << *aux1_ << std::endl;
    os << std::endl << "AuxCoupledModel: " << MODEL2::name() << std::endl;
    os << *aux2_;
  }

  Log::trace() << "AuxCoupledModel::print done" << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace oops
