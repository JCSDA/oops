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

// QgObsTraits is defined in QgObsTraitsFwd.h, which, however,
// contains only forward declarations of the QG implementations of oops interfaces.
// This file includes headers in which all these implementations are defined.

#include "model/GomQG.h"
#include "model/LocationsQG.h"
#include "model/ObsBias.h"
#include "model/ObsBiasCovariance.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsBiasPreconditioner.h"
#include "model/ObsDataQG.h"
#include "model/ObsDiagsQG.h"
#include "model/ObsErrorDiagQG.h"
#include "model/ObsFilter.h"
#include "model/ObsIteratorQG.h"
#include "model/ObsOperatorQG.h"
#include "model/ObsOperatorTLAD.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgObsTraitsFwd.h"

