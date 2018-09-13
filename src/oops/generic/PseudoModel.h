/*
 * (C) Copyright 2018 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PSEUDOMODEL_H_
#define OOPS_GENERIC_PSEUDOMODEL_H_

#include <string>

#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelBase.h"
#include "oops/interface/State.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Encapsulates a pseudo forecast model.
/*!
 * Generic implementation of the pseudo model.
 */

// -----------------------------------------------------------------------------

template <typename MODEL>
class PseudoModel : public ModelBase<MODEL> {
  typedef typename MODEL::Geometry          Geometry_;
  typedef typename MODEL::ModelAuxControl   ModelAux_;
  typedef typename MODEL::State             State_;

 public:
  static const std::string classname() {return "oops::PseudoModel";}

  PseudoModel(const Geometry_ &, const eckit::Configuration &);
  ~PseudoModel();

// Run the Pseudo forecast
  void initialize(State_ &) const override;
  void step(State_ &, const ModelAux_ &) const override;
  void finalize(State_ &) const override;

// Information and diagnostics
  const util::Duration & timeResolution() const override {return tstep_;}
  void print(std::ostream &) const override {}

 private:
  const Geometry_ resol_;
  const util::Duration tstep_;
};

// =============================================================================

template<typename MODEL>
PseudoModel<MODEL>::PseudoModel(const Geometry_ & resol, const eckit::Configuration & tlConf)
  : resol_(resol), tstep_(util::Duration(tlConf.getString("tstep")))
{
  Log::trace() << "PseudoModel<MODEL>::PseudoModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoModel<MODEL>::~PseudoModel() {
  Log::trace() << "PseudoModel<MODEL>::~PseudoModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::initialize(State_ & xx) const {
  Log::info() << "PseudoModel<MODEL>:initialize Starting " << std::endl;
  Log::trace() << "PseudoModel<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::info() << "PseudoModel<MODEL>:step Starting " << std::endl;
  xx.updateTime(tstep_);
  Log::trace() << "PseudoModel<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::finalize(State_ & xx) const {
  Log::info() << "PseudoModel<MODEL>:finalize Starting " << std::endl;
  Log::trace() << "PseudoModel<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_PSEUDOMODEL_H_
