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

#include "oops/base/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "oops/qg/GeometryQG.h"
#include "oops/qg/QgFortran.h"
#include "oops/qg/QgTraits.h"


namespace qg {
  class ModelBias;
  class FieldsQG;
  class StateQG;

// -----------------------------------------------------------------------------

class ModelQgParameters : public oops::ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(ModelQgParameters, ModelParametersBase)

 public:
  oops::RequiredParameter<util::Duration> tstep{"tstep", this};
};


// -----------------------------------------------------------------------------
/// QG model definition.
/*!
 *  QG nonlinear model definition and configuration parameters.
 */

class ModelQG: public oops::ModelBase<QgTraits>,
               private util::ObjectCounter<ModelQG> {
 public:
  typedef ModelQgParameters Parameters_;

  static const std::string classname() {return "qg::ModelQG";}

  ModelQG(const GeometryQG &, const ModelQgParameters &);
  ~ModelQG();

/// Prepare model integration
  void initialize(StateQG &) const;

/// Model integration
  void step(StateQG &, const ModelBias &) const;
  int saveTrajectory(StateQG &, const ModelBias &) const;

/// Finish model integration
  void finalize(StateQG &) const;

/// Utilities
  const util::Duration & timeResolution() const {return params_.tstep;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  ModelQgParameters params_;
  const GeometryQG geom_;
  const oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_MODELQG_H_
