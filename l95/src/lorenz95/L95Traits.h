/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_L95TRAITS_H_
#define LORENZ95_L95TRAITS_H_

#include <string>

#include "lorenz95/ErrorCovarianceL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/LocalizationMatrixL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/NoVariables.h"
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "lorenz95/ObsBiasCovariance.h"
#include "lorenz95/ObservationL95.h"
#include "lorenz95/ObservationTLAD.h"
#include "lorenz95/ObsTable.h"
#include "lorenz95/ObsVec1D.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"

namespace lorenz95 {

struct L95Traits {
  static std::string name() {return "Lorenz 95";}

  typedef lorenz95::Resolution             Geometry;
  typedef lorenz95::NoVariables            Variables;

  typedef lorenz95::StateL95               State;
  typedef lorenz95::ModelL95               Model;
  typedef lorenz95::IncrementL95           Increment;
  typedef lorenz95::ErrorCovarianceL95     Covariance;

  typedef lorenz95::ModelBias              ModelAuxControl;
  typedef lorenz95::ModelBiasCorrection    ModelAuxIncrement;
  typedef lorenz95::ModelBiasCovariance    ModelAuxCovariance;

  typedef lorenz95::ObsTable               ObsSpace;
  typedef lorenz95::ObservationL95         ObsOperator;
  typedef lorenz95::ObservationTLAD        LinearObsOperator;
  typedef lorenz95::ObsVec1D               ObsVector;

  typedef lorenz95::ObsBias                ObsAuxControl;
  typedef lorenz95::ObsBiasCorrection      ObsAuxIncrement;
  typedef lorenz95::ObsBiasCovariance      ObsAuxCovariance;

  typedef lorenz95::GomL95                 ModelAtLocations;
  typedef lorenz95::LocsL95                Locations;

  typedef lorenz95::LocalizationMatrixL95  LocalizationMatrix;
};

}  // namespace lorenz95

#endif  // LORENZ95_L95TRAITS_H_
