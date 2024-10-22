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

// The QgTraits and QgObsTraits classes are defined in QgTraitsFwd.h, which, however,
// contains only forward declarations of the QG implementations of oops interfaces.
// This file includes headers in which all these implementations are defined.

#include "oops/qg/ChangeVarQG.h"
#include "oops/qg/ChangeVarTLADQG.h"
#include "oops/qg/ErrorCovarianceQG.h"
#include "oops/qg/GeometryQG.h"
#include "oops/qg/GeometryQGIterator.h"
#include "oops/qg/GomQG.h"
#include "oops/qg/IncrementQG.h"
#include "oops/qg/InterpolatorQG.h"
#include "oops/qg/LocationsQG.h"
#include "oops/qg/ModelBias.h"
#include "oops/qg/ModelBiasCovariance.h"
#include "oops/qg/ModelBiasIncrement.h"
#include "oops/qg/ModelData.h"
#include "oops/qg/ObsBias.h"
#include "oops/qg/ObsBiasCovariance.h"
#include "oops/qg/ObsBiasIncrement.h"
#include "oops/qg/ObsBiasPreconditioner.h"
#include "oops/qg/ObsDataQG.h"
#include "oops/qg/ObsDiagsQG.h"
#include "oops/qg/ObsFilter.h"
#include "oops/qg/ObsIteratorQG.h"
#include "oops/qg/ObsOperatorQG.h"
#include "oops/qg/ObsOperatorTLAD.h"
#include "oops/qg/ObsSpaceQG.h"
#include "oops/qg/ObsVecQG.h"
#include "oops/qg/QgTraitsFwd.h"
#include "oops/qg/StateQG.h"

#endif  // QG_MODEL_QGTRAITS_H_
