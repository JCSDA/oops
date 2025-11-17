/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "model/ChangeVarQG.h"
#include "model/ChangeVarTLADQG.h"
#include "model/ErrorCovarianceQG.h"
#include "model/GeometryQG.h"
#include "model/GeometryQGIterator.h"
#include "model/IncrementQG.h"
#include "model/InterpolatorQG.h"
#include "model/LocalizationMatrixQG.h"
#include "model/ModelBias.h"
#include "model/ModelBiasCovariance.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelData.h"
#include "model/ModelQG.h"
#include "model/ObsLocQG.h"
#include "model/StateQG.h"
#include "model/TlmQG.h"

namespace qg {

struct QgTraits {
  static std::string name() {return "QG";}
  static std::string nameCovar() {return "QgError";}
  static std::string nameCovar4D() {return "QgError";}

  typedef qg::GeometryQG            Geometry;

  typedef qg::GeometryQGIterator    GeometryIterator;

  typedef qg::ModelQG               Model;
  typedef qg::TlmQG                 LinearModel;

  typedef qg::ChangeVarQG           VariableChange;
  typedef qg::ChangeVarTLADQG       LinearVariableChange;

  typedef qg::StateQG               State;
  typedef qg::IncrementQG           Increment;
  typedef qg::ErrorCovarianceQG     Covariance;
  typedef qg::InterpolatorQG        LocalInterpolator;
  typedef qg::LocalizationMatrixQG  Localization;

  typedef qg::ModelBias             ModelAuxControl;
  typedef qg::ModelBiasIncrement    ModelAuxIncrement;
  typedef qg::ModelBiasCovariance   ModelAuxCovariance;
  typedef qg::ModelData             ModelData;

  typedef qg::ObsLocQG              ObsLocalization;
};

}  // namespace qg
