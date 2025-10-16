/*
 * (C) Copyright 2023- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "model/ModelQG.h"
#include "model/TlmQG.h"

#include "model/ChangeVarQG.h"
#include "model/ChangeVarTLADQG.h"
#include "model/ErrorCovarianceIdQG.h"
#include "model/GeometryQG.h"
#include "model/GeometryQGIterator.h"
#include "model/IncrementQG.h"
#include "model/InterpolatorQG.h"
#include "model/ModelBias.h"
#include "model/ModelBiasCovariance.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelData.h"
#include "model/StateQG.h"

namespace qg {

// copy of QgTraits with a different name for use with qg-qg
// coupled applications
struct QgTraits2 {
  static std::string name() {return "QG 2";}
  static std::string nameCovar() {return "QgErrorId";}
  static std::string nameCovar4D() {return "QgErrorId";}

  typedef qg::GeometryQG            Geometry;

  typedef qg::GeometryQGIterator    GeometryIterator;

  typedef qg::ModelQG               Model;
  typedef qg::TlmQG                 LinearModel;

  typedef qg::ChangeVarQG           VariableChange;
  typedef qg::ChangeVarTLADQG       LinearVariableChange;

  typedef qg::StateQG               State;
  typedef qg::IncrementQG           Increment;
  typedef qg::ErrorCovarianceIdQG   Covariance;
  typedef qg::InterpolatorQG        LocalInterpolator;

  typedef qg::ModelBias             ModelAuxControl;
  typedef qg::ModelBiasIncrement    ModelAuxIncrement;
  typedef qg::ModelBiasCovariance   ModelAuxCovariance;
  typedef qg::ModelData             ModelData;
};

}  // namespace qg
