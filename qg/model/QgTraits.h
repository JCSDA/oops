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

#include "model/ChangeVarQG.h"
#include "model/ErrorCovarianceQG.h"
#include "model/GeometryQG.h"
#include "model/GeometryQGIterator.h"
#include "model/GetValuesQG.h"
#include "model/GetValuesTLAD.h"
#include "model/GomQG.h"
#include "model/IncrementQG.h"
#include "model/LocationsQG.h"
#include "model/ModelBias.h"
#include "model/ModelBiasCovariance.h"
#include "model/ModelBiasIncrement.h"
#include "model/ObsBias.h"
#include "model/ObsBiasCovariance.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsDataQG.h"
#include "model/ObsDiagsQG.h"
#include "model/ObsOperatorQG.h"
#include "model/ObsOperatorTLAD.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/StateQG.h"

namespace qg {

struct QgTraits {
  static std::string name() {return "QG";}
  static std::string nameCovar() {return "QgError";}
  static std::string nameCovar4D() {return "QgError";}

  typedef qg::GeometryQG            Geometry;

  typedef qg::GeometryQGIterator    GeometryIterator;

  typedef qg::GetValuesQG           GetValues;
  typedef qg::GetValuesTLAD         LinearGetValues;

  typedef qg::StateQG               State;
  typedef qg::IncrementQG           Increment;
  typedef qg::ErrorCovarianceQG     Covariance;

  typedef qg::ModelBias             ModelAuxControl;
  typedef qg::ModelBiasIncrement    ModelAuxIncrement;
  typedef qg::ModelBiasCovariance   ModelAuxCovariance;
};

struct QgObsTraits {
  static std::string name() {return "QG obs";}

  typedef qg::ObsSpaceQG            ObsSpace;
  typedef qg::ObsVecQG              ObsVector;
  typedef qg::ObsOperatorQG         ObsOperator;
  typedef qg::ObsOperatorTLAD       LinearObsOperator;
  template <typename DATATYPE> using ObsDataVector = qg::ObsDataQG<DATATYPE>;

  typedef qg::ObsBias               ObsAuxControl;
  typedef qg::ObsBiasIncrement      ObsAuxIncrement;
  typedef qg::ObsBiasCovariance     ObsAuxCovariance;

  typedef qg::ObsDiagsQG            ObsDiagnostics;

  typedef qg::GomQG                 GeoVaLs;
  typedef qg::LocationsQG           Locations;
};

}  // namespace qg

#endif  // QG_MODEL_QGTRAITS_H_
