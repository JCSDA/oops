/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_QGTRAITSFWD_H_
#define QG_MODEL_QGTRAITSFWD_H_

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

class ObsSpaceQG;
class ObsVecQG;
template <typename DATATYPE> class ObsDataQG;
class ObsIteratorQG;

class ObsOperatorQG;
class ObsOperatorTLAD;
class ObsBias;
class ObsBiasIncrement;
class ObsBiasCovariance;
class ObsBiasPreconditioner;
class ObsDiagsQG;
class ObsFilter;

class GomQG;
class LocationsQG;

struct QgTraits {
  static std::string name() {return "QG";}
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

struct QgObsTraits {
  static std::string name() {return "QG obs";}

  typedef qg::ObsSpaceQG            ObsSpace;
  typedef qg::ObsVecQG              ObsVector;
  typedef qg::ObsOperatorQG         ObsOperator;
  typedef qg::ObsOperatorTLAD       LinearObsOperator;
  template <typename DATATYPE> using ObsDataVector = qg::ObsDataQG<DATATYPE>;
  typedef qg::ObsIteratorQG         GeometryIterator;

  typedef qg::ObsBias               ObsAuxControl;
  typedef qg::ObsBiasIncrement      ObsAuxIncrement;
  typedef qg::ObsBiasCovariance     ObsAuxCovariance;
  typedef qg::ObsBiasPreconditioner ObsAuxPreconditioner;

  typedef qg::ObsDiagsQG            ObsDiagnostics;
  typedef qg::ObsFilter             ObsFilter;

  typedef qg::GomQG                 GeoVaLs;
  typedef qg::LocationsQG           SampledLocations;
};

}  // namespace qg

#endif  // QG_MODEL_QGTRAITSFWD_H_
