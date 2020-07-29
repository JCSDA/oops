/*
 * (C) Copyright 2018-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PSEUDOMODEL_H_
#define OOPS_GENERIC_PSEUDOMODEL_H_

#include <string>
#include <vector>

#include "oops/base/ModelBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

///  Generic implementation of the pseudo model (steps through time by reading states)
template <typename MODEL>
class PseudoModel : public ModelBase<MODEL> {
  typedef typename MODEL::Geometry          Geometry_;
  typedef typename MODEL::ModelAuxControl   ModelAux_;
  typedef typename MODEL::State             State_;

 public:
  static const std::string classname() {return "oops::PseudoModel";}

  PseudoModel(const Geometry_ &, const eckit::Configuration &);

/// initialize forecast
  void initialize(State_ &) const override;
/// one forecast step
  void step(State_ &, const ModelAux_ &) const override;
/// finalize forecast
  void finalize(State_ &) const override;

/// model time step
  const util::Duration & timeResolution() const override {return tstep_;}

/// model variables
  const oops::Variables & variables() const override {return vars_;}

 private:
  void print(std::ostream &) const override {}
  const util::Duration tstep_;
  const oops::Variables vars_;
  std::vector<eckit::LocalConfiguration> confs_;
  mutable size_t currentstate_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoModel<MODEL>::PseudoModel(const Geometry_ & resol, const eckit::Configuration & conf)
  : tstep_(util::Duration(conf.getString("tstep"))), vars_(conf, "state variables"),
    confs_(conf.getSubConfigurations("states")),
    currentstate_(0) {
  Log::trace() << "PseudoModel<MODEL>::PseudoModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::initialize(State_ & xx) const {
  currentstate_ = 0;
  Log::trace() << "PseudoModel<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "PseudoModel<MODEL>:step Starting " << std::endl;
  xx.updateTime(tstep_);
  xx.read(confs_[currentstate_++]);
  Log::trace() << "PseudoModel<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::finalize(State_ & xx) const {
  Log::trace() << "PseudoModel<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_PSEUDOMODEL_H_
