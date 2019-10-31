/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_CONTROLVARIABLE_H_
#define OOPS_ASSIMILATION_CONTROLVARIABLE_H_

#include <cmath>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

/// Control variable
/*!
 * The control variable acts as a container for the inputs of the variational
 * data assimilation cost functions in physical space.
 * That includes the states at the start the assimilation window or of each
 * sub-window but also additional variables such as model bias, VarBC
 * coefficients, or other control variables for algorithms that use them.
 * This is mostly a convenience class that is used to keep things together
 * and reduce the number of arguments to be passed around.
 */

template<typename MODEL> class ControlVariable;

// -----------------------------------------------------------------------------
template<typename MODEL>
class ControlVariable : public util::Printable,
                        private util::ObjectCounter<ControlVariable<MODEL> > {
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef State4D<MODEL>             State4D_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::ControlVariable";}

/// The arguments define the number of sub-windows and the resolution
  ControlVariable(const eckit::Configuration &, const Variables &, const Geometry_ &);
  ControlVariable(const eckit::Configuration &, const State_ &);
  explicit ControlVariable(const ControlVariable &);
  ~ControlVariable();

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

/// Get state control variable
  State4D_ & state() {return state4d_;}
  const State4D_ & state() const {return state4d_;}

/// Get augmented model control variable
  ModelAux_ & modVar() {return modbias_;}
  const ModelAux_ & modVar() const {return modbias_;}

/// Get augmented observation control variable
  ObsAuxCtrls_ & obsVar() {return obsbias_;}
  const ObsAuxCtrls_ & obsVar() const {return obsbias_;}

 private:
  ControlVariable & operator= (const ControlVariable &);  // No assignment
  void print(std::ostream &) const;

  State4D_ state4d_;
  ModelAux_ modbias_;     // not only for bias, better name?
  ObsAuxCtrls_ obsbias_;  // not only for bias, better name?
};

// =============================================================================

template<typename MODEL>
ControlVariable<MODEL>::ControlVariable(const eckit::Configuration & conf,
                                        const Variables & vars, const Geometry_ & resol)
  : state4d_(eckit::LocalConfiguration(conf, "Jb.Background"), vars, resol),
    modbias_(resol, conf.getSubConfiguration("Jb.Background.ModelBias")),
    obsbias_(eckit::LocalConfiguration(conf, "Jo"))
{
  Log::trace() << "ControlVariable contructed" << std::endl;
}

// =============================================================================

template<typename MODEL>
ControlVariable<MODEL>::ControlVariable(const eckit::Configuration & conf,
                                        const State_ & statein)
  : state4d_(statein),
    modbias_(statein.geometry(), conf.getSubConfiguration("Jb.Background.ModelBias")),
    obsbias_(eckit::LocalConfiguration(conf, "Jo"))
{
  Log::trace() << "ControlVariable contructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ControlVariable<MODEL>::ControlVariable(const ControlVariable & other)
  : state4d_(other.state4d_), modbias_(other.modbias_), obsbias_(other.obsbias_)
{
  Log::trace() << "ControlVariable copied" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ControlVariable<MODEL>::~ControlVariable() {
  Log::trace() << "ControlVariable destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ControlVariable<MODEL>::read(const eckit::Configuration & config) {
  state4d_.read(config);
  modbias_.read(config);
  obsbias_.read(config);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ControlVariable<MODEL>::write(const eckit::Configuration & config) const {
  state4d_.write(config);
  modbias_.write(config);
  obsbias_.write(config);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ControlVariable<MODEL>::print(std::ostream & outs) const {
  outs << state4d_;
  outs << modbias_;
  outs << obsbias_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double ControlVariable<MODEL>::norm() const {
  double zz = state4d_.norm();
  double zn = zz * zz;
  zz = modbias_.norm();
  zn += zz * zz;
  zz = obsbias_.norm();
  zn += zz * zz;
  return sqrt(zn);
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_CONTROLVARIABLE_H_
