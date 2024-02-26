/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_MODELQG_H_
#define QG_MODEL_MODELQG_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"

#include "oops/qg/GeometryQG.h"
#include "oops/qg/QgFortran.h"
#include "oops/qg/QgTraits.h"

namespace eckit {
  class Configuration;
}

namespace qg {
  class ModelBias;
  class FieldsQG;
  class StateQG;

// -----------------------------------------------------------------------------
/// QG nonlinear model definition.

class ModelQG: public oops::interface::ModelBase<QgTraits>,
               private util::ObjectCounter<ModelQG> {
 public:
  static const std::string classname() {return "qg::ModelQG";}

  ModelQG(const GeometryQG &, const eckit::Configuration &);
  ~ModelQG();

/// Prepare model integration
  void initialize(StateQG &) const;

/// Model integration
  void step(StateQG &, const ModelBias &) const;
  int saveTrajectory(StateQG &, const ModelBias &) const;

/// Finish model integration
  void finalize(StateQG &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  const util::Duration tstep_;
  const GeometryQG geom_;
  oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_MODELQG_H_
