/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_MODELL95_H_
#define LORENZ95_MODELL95_H_

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "lorenz95/L95Traits.h"
#include "lorenz95/Resolution.h"

namespace lorenz95 {
  class FieldL95;
  class ModelBias;
  class ModelTrajectory;
  class StateL95;

// -----------------------------------------------------------------------------

class ModelL95 : public oops::interface::ModelBase<L95Traits>,
                 private util::ObjectCounter<ModelL95> {
 public:
  static const std::string classname() {return "lorenz95::ModelL95";}

  ModelL95(const Resolution &, const eckit::Configuration &);
  ~ModelL95();

// Run the forecast
  void initialize(StateL95 &) const;
  void step(StateL95 &, const ModelBias &) const;
  void finalize(StateL95 &) const;
  void stepRK(FieldL95 &, const ModelBias &, ModelTrajectory &) const;

// Information and diagnostics
  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  void tendencies(const FieldL95 &, const double &, FieldL95 &) const;

// Data
  const Resolution resol_;
  const double f_;
  const util::Duration tstep_;
  const double dt_;
  const oops::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_MODELL95_H_
