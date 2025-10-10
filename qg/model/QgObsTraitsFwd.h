/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <string>

namespace qg {

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

