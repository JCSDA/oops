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

#include "lorenz95/ErrorCovarianceL95.h"
#include "lorenz95/IdChangeVariable.h"
#include "lorenz95/IdChangeVarTLADL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/InterpolatorL95.h"
#include "lorenz95/Iterator.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "lorenz95/ModelData.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/NormGradientL95.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "lorenz95/TLML95.h"

namespace lorenz95 {

struct L95Traits {
  static std::string name() {return "Lorenz 95";}
  static std::string nameCovar() {return "L95Error";}
  static std::string nameCovar4D() {return "L95Error";}

  typedef lorenz95::Resolution             Geometry;
  typedef lorenz95::Iterator               GeometryIterator;

  typedef lorenz95::StateL95               State;
  typedef lorenz95::IncrementL95           Increment;
  typedef lorenz95::ErrorCovarianceL95     Covariance;
  typedef lorenz95::InterpolatorL95        LocalInterpolator;

  typedef lorenz95::ModelL95               Model;
  typedef lorenz95::TLML95                 LinearModel;

  typedef lorenz95::IdChangeVariable       VariableChange;
  typedef lorenz95::IdChangeVarTLADL95     LinearVariableChange;

  typedef lorenz95::NormGradientL95        NormGradient;

  typedef lorenz95::ModelBias              ModelAuxControl;
  typedef lorenz95::ModelBiasCorrection    ModelAuxIncrement;
  typedef lorenz95::ModelBiasCovariance    ModelAuxCovariance;
  typedef lorenz95::ModelData              ModelData;
};

}  // namespace lorenz95

