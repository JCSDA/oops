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
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/ModelBias.h"
#include "model/QgFortran.h"
#include "model/StateQG.h"


namespace qg {
// -----------------------------------------------------------------------------
static oops::ModelMaker<QgTraits, ModelQG> makermodel_("QG");
// -----------------------------------------------------------------------------
ModelQG::ModelQG(const GeometryQG & resol, const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(model)
{
  oops::Log::trace() << "ModelQG::ModelQG" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  qg_model_setup_f90(keyConfig_, model);
  oops::Log::trace() << "ModelQG created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelQG::~ModelQG() {
  qg_model_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelQG::initialize(StateQG & xx) const {
  ASSERT(xx.fields().isForModel(true));
  oops::Log::debug() << "ModelQG::initialize" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void ModelQG::step(StateQG & xx, const ModelBias &) const {
  ASSERT(xx.fields().isForModel(true));
  oops::Log::debug() << "ModelQG::step fields in" << xx.fields() << std::endl;
  qg_model_propagate_f90(keyConfig_, xx.fields().toFortran());
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
  qg_fields_create_from_other_f90(ftraj, xx.fields().toFortran());
  qg_fields_copy_f90(ftraj, xx.fields().toFortran());
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
