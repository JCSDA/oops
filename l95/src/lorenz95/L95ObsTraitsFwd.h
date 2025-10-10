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

namespace lorenz95 {

class ObsTable;
class ObsVec1D;
template <typename DATATYPE> class ObsData1D;
class ObsIterator;

class ObservationL95;
class ObservationTLAD;
class ObsBias;
class ObsBiasCorrection;
class ObsBiasCovariance;
class ObsBiasPreconditioner;
class ObsDiags1D;
class ObsFilter;

class GomL95;
class LocsL95;

struct L95ObsTraits {
  static std::string name() {return "Lorenz 95 Obs";}

  typedef lorenz95::ObsTable               ObsSpace;
  typedef lorenz95::ObsVec1D               ObsVector;
  template <typename DATATYPE> using ObsDataVector = lorenz95::ObsData1D<DATATYPE>;
  typedef lorenz95::ObsIterator            GeometryIterator;

  typedef lorenz95::ObservationL95         ObsOperator;
  typedef lorenz95::ObservationTLAD        LinearObsOperator;
  typedef lorenz95::ObsBias                ObsAuxControl;
  typedef lorenz95::ObsBiasCorrection      ObsAuxIncrement;
  typedef lorenz95::ObsBiasCovariance      ObsAuxCovariance;
  typedef lorenz95::ObsBiasPreconditioner  ObsAuxPreconditioner;
  typedef lorenz95::ObsDiags1D             ObsDiagnostics;
  typedef lorenz95::ObsFilter              ObsFilter;

  typedef lorenz95::GomL95                 GeoVaLs;
  typedef lorenz95::LocsL95                SampledLocations;
};

}  // namespace lorenz95

