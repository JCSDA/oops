/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_QGTRAITS_H_
#define QG_MODEL_QGTRAITS_H_

#include <string>

#include "model/ErrorCovarianceQG.h"
#include "model/GeometryQG.h"
#include "model/GomQG.h"
#include "model/IncrementQG.h"
#include "model/LinearObsOp.h"
#include "model/LocalizationMatrixQG.h"
#include "model/LocationsQG.h"
#include "model/ModelQG.h"
#include "model/ModelBias.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelBiasCovariance.h"
#include "model/ObservationsQG.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsBiasCovariance.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/StateQG.h"
#include "model/VariablesQG.h"

namespace qg {

struct QgTraits {
  static std::string name() {return "QG";}

  typedef qg::GeometryQG            Geometry;
  typedef qg::VariablesQG           Variables;

  typedef qg::StateQG               State;
  typedef qg::ModelQG               Model;
  typedef qg::IncrementQG           Increment;
  typedef qg::ErrorCovarianceQG     Covariance;

  typedef qg::ModelBias             ModelAuxControl;
  typedef qg::ModelBiasIncrement    ModelAuxIncrement;
  typedef qg::ModelBiasCovariance   ModelAuxCovariance;

  typedef qg::ObsSpaceQG            ObsSpace;
  typedef qg::ObservationsQG        ObsOperator;
  typedef qg::LinearObsOp           LinearObsOperator;
  typedef qg::ObsVecQG              ObsVector;

  typedef qg::ObsBias               ObsAuxControl;
  typedef qg::ObsBiasIncrement      ObsAuxIncrement;
  typedef qg::ObsBiasCovariance     ObsAuxCovariance;

  typedef qg::GomQG                 GeoVaLs;
  typedef qg::LocationsQG           Locations;

  typedef qg::LocalizationMatrixQG  LocalizationMatrix;
};

}  // namespace qg

#endif  // QG_MODEL_QGTRAITS_H_
