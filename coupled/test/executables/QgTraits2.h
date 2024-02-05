/*
 * (C) Copyright 2023- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

namespace qg {

class GeometryQG;
class GeometryQGIterator;

class StateQG;
class IncrementQG;
class ErrorCovarianceQG;
class InterpolatorQG;

class ChangeVarQG;
class ChangeVarTLADQG;

class ModelBias;
class ModelBiasIncrement;
class ModelBiasCovariance;
class ModelData;

// copy of QgTraits with a different name for use with qg-qg
// coupled applications
struct QgTraits2 {
  static std::string name() {return "QG 2";}
  static std::string nameCovar() {return "QgError";}
  static std::string nameCovar4D() {return "QgError";}

  typedef qg::GeometryQG            Geometry;

  typedef qg::GeometryQGIterator    GeometryIterator;

  typedef qg::ChangeVarQG           VariableChange;
  typedef qg::ChangeVarTLADQG       LinearVariableChange;

  typedef qg::StateQG               State;
  typedef qg::IncrementQG           Increment;
  typedef qg::ErrorCovarianceQG     Covariance;
  typedef qg::InterpolatorQG        LocalInterpolator;

  typedef qg::ModelBias             ModelAuxControl;
  typedef qg::ModelBiasIncrement    ModelAuxIncrement;
  typedef qg::ModelBiasCovariance   ModelAuxCovariance;
  typedef qg::ModelData             ModelData;
};

}  // namespace qg
