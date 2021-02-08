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

#include "oops/qg/AnalyticInit.h"
#include "oops/qg/ErrorCovarianceQG.h"
#include "oops/qg/GeometryQG.h"
#include "oops/qg/GeometryQGIterator.h"
#include "oops/qg/GetValuesQG.h"
#include "oops/qg/GetValuesTLAD.h"
#include "oops/qg/GomQG.h"
#include "oops/qg/IncrementQG.h"
#include "oops/qg/LocationsQG.h"
#include "oops/qg/ModelBias.h"
#include "oops/qg/ModelBiasCovariance.h"
#include "oops/qg/ModelBiasIncrement.h"
#include "oops/qg/ObsBias.h"
#include "oops/qg/ObsBiasCovariance.h"
#include "oops/qg/ObsBiasIncrement.h"
#include "oops/qg/ObsDataQG.h"
#include "oops/qg/ObsDiagsQG.h"
#include "oops/qg/ObsOperatorQG.h"
#include "oops/qg/ObsOperatorTLAD.h"
#include "oops/qg/ObsSpaceQG.h"
#include "oops/qg/ObsVecQG.h"
#include "oops/qg/StateQG.h"

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
  typedef qg::AnalyticInit          AnalyticInit;
};

}  // namespace qg

#endif  // QG_MODEL_QGTRAITS_H_
