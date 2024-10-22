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

#endif  // LORENZ95_L95TRAITS_H_
