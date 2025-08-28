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

// The L95Traits and L95ObsTraits classes are defined in L95TraitsFwd.h, which, however,
// contains only forward declarations of the lorenz95 implementations of oops interfaces.
// This file includes headers in which all these implementations are defined.

#include <string>

#include "lorenz95/ErrorCovarianceL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IdChangeVariable.h"
#include "lorenz95/IdChangeVarTLADL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/InterpolatorL95.h"
#include "lorenz95/Iterator.h"
#include "lorenz95/L95TraitsFwd.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "lorenz95/ModelData.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/NormGradientL95.h"
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "lorenz95/ObsBiasCovariance.h"
#include "lorenz95/ObsBiasPreconditioner.h"
#include "lorenz95/ObsData1D.h"
#include "lorenz95/ObsDiags1D.h"
#include "lorenz95/ObservationL95.h"
#include "lorenz95/ObservationTLAD.h"
#include "lorenz95/ObsFilter.h"
#include "lorenz95/ObsIterator.h"
#include "lorenz95/ObsTable.h"
#include "lorenz95/ObsVec1D.h"
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

#endif  // LORENZ95_L95TRAITS_H_
