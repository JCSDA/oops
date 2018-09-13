/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ModelQG.h"

#include <vector>

#include "eckit/config/Configuration.h"

#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/ModelBias.h"
#include "model/QgFortran.h"
#include "model/StateQG.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
static oops::ModelMaker<QgTraits, ModelQG> makermodel_("QG");
// -----------------------------------------------------------------------------
ModelQG::ModelQG(const GeometryQG & resol, const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol),
    vars_(std::vector<std::string>{"x", "q", "u", "v", "bc"})
{
  oops::Log::trace() << "ModelQG::ModelQG" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  const eckit::Configuration * configc = &model;
  qg_setup_f90(&configc, geom_.toFortran(), keyConfig_);
  oops::Log::trace() << "ModelQG created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelQG::~ModelQG() {
  qg_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelQG::initialize(StateQG & xx) const {
  ASSERT(xx.fields().isForModel(true));
  qg_prepare_integration_f90(keyConfig_, xx.fields().toFortran());
  oops::Log::debug() << "ModelQG::initialize" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void ModelQG::step(StateQG & xx, const ModelBias &) const {
  ASSERT(xx.fields().isForModel(true));
  oops::Log::debug() << "ModelQG::step fields in" << xx.fields() << std::endl;
  qg_propagate_f90(keyConfig_, xx.fields().toFortran());
  xx.validTime() += tstep_;
  oops::Log::debug() << "ModelQG::step fields out" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void ModelQG::finalize(StateQG & xx) const {
  ASSERT(xx.fields().isForModel(true));
  oops::Log::debug() << "ModelQG::finalize" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
int ModelQG::saveTrajectory(StateQG & xx, const ModelBias &) const {
  ASSERT(xx.fields().isForModel(true));
  int ftraj = 0;
  oops::Log::debug() << "ModelQG::saveTrajectory fields in" << xx.fields() << std::endl;
  qg_prop_traj_f90(keyConfig_, xx.fields().toFortran(), ftraj);
  ASSERT(ftraj != 0);
  oops::Log::debug() << "ModelQG::saveTrajectory fields out" << xx.fields() << std::endl;
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelQG::print(std::ostream & os) const {
  os << "ModelQG::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace qg
