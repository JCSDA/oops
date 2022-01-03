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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

#include "oops/coupled/GeometryCoupled.h"

namespace oops {

/// Parameters for ModelAuxControl describing a coupled model bias
template <typename MODEL1, typename MODEL2>
class AuxCoupledModelParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(AuxCoupledModelParameters, Parameters)

  typedef typename ModelAuxControl<MODEL1>::Parameters_ Parameters1_;
  typedef typename ModelAuxControl<MODEL2>::Parameters_ Parameters2_;
 public:
  /// Parameters for ModelAuxControl of MODEL1 and ModelAuxControl of MODEL2
  Parameter<Parameters1_> modelaux1{MODEL1::name().c_str(), {}, this};
  Parameter<Parameters2_> modelaux2{MODEL2::name().c_str(), {}, this};
};


// -----------------------------------------------------------------------------
/// Implementation of the ModelAuxControl interface for a coupled model.
template <typename MODEL1, typename MODEL2>
class AuxCoupledModel : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;

 public:
  typedef AuxCoupledModelParameters<MODEL1, MODEL2> Parameters_;

  AuxCoupledModel(const GeometryCoupled_ &, const Parameters_ &);
  AuxCoupledModel(const GeometryCoupled_ &, const AuxCoupledModel &);
  AuxCoupledModel(const AuxCoupledModel &, const bool);
  ~AuxCoupledModel();

  /// Accessors to the individual components of the coupled ModelAuxControl
  ModelAuxControl<MODEL1> & aux1() {return *aux1_;}
  ModelAuxControl<MODEL2> & aux2() {return *aux2_;}
  const ModelAuxControl<MODEL1> & aux1() const {return *aux1_;}
  const ModelAuxControl<MODEL2> & aux2() const {return *aux2_;}

  /// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ModelAuxControl<MODEL1>> aux1_;
  std::unique_ptr<ModelAuxControl<MODEL2>> aux2_;
};

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2>::AuxCoupledModel(const GeometryCoupled_ & geom,
                                                 const Parameters_ & params)
  : aux1_(), aux2_()
{
  Log::trace() << "AuxCoupledModel::AuxCoupledModel read starting" << std::endl;

  aux1_ = std::make_unique<ModelAuxControl<MODEL1>>(geom.geometry1(), params.modelaux1);
  aux2_ = std::make_unique<ModelAuxControl<MODEL2>>(geom.geometry2(), params.modelaux2);

  Log::trace() << "AuxCoupledModel::AuxCoupledModel read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2>::AuxCoupledModel(const GeometryCoupled_ & geom,
                                                 const AuxCoupledModel & other)
  : aux1_(), aux2_()
{
  Log::trace() << "AuxCoupledModel::AuxCoupledModel interpolated starting" << std::endl;
  aux1_.reset(new ModelAuxControl<MODEL1>(geom.geometry1(), *other.aux1_));
  aux2_.reset(new ModelAuxControl<MODEL2>(geom.geometry2(), *other.aux2_));
  Log::trace() << "AuxCoupledModel::AuxCoupledModel interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
AuxCoupledModel<MODEL1, MODEL2>::AuxCoupledModel(const AuxCoupledModel & other, const bool copy)
  : aux1_(), aux2_()
{
  Log::trace() << "AuxCoupledModel::AuxCoupledModel starting copy" << std::endl;
  aux1_.reset(new ModelAuxControl<MODEL1>(*other.aux1_, copy));
  aux2_.reset(new ModelAuxControl<MODEL2>(*other.aux2_, copy));
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
  aux1_->read(conf.getSubConfiguration(MODEL1::name()));
  aux2_->read(conf.getSubConfiguration(MODEL2::name()));
  Log::trace() << "AuxCoupledModel::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void AuxCoupledModel<MODEL1, MODEL2>::write(const eckit::Configuration & conf) const {
  Log::trace() << "AuxCoupledModel::write starting" << std::endl;
  aux1_->write(conf.getSubConfiguration(MODEL1::name()));
  aux2_->write(conf.getSubConfiguration(MODEL2::name()));
  Log::trace() << "AuxCoupledModel::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
double AuxCoupledModel<MODEL1, MODEL2>::norm() const {
  Log::trace() << "AuxCoupledModel::norm starting" << std::endl;
  const double zz = aux1_->norm() + aux2_->norm();
  Log::trace() << "AuxCoupledModel::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void AuxCoupledModel<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "AuxCoupledModel::print starting" << std::endl;
  os << "AuxCoupledModel: " << MODEL1::name() << std::endl;
  os << *aux1_ << std::endl;
  os << "AuxCoupledModel: " << MODEL2::name() << std::endl;
  os << *aux2_;
  Log::trace() << "AuxCoupledModel::print done" << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace oops
