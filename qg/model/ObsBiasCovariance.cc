/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsBiasCovariance.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>

#include "model/ObsBiasIncrement.h"
#include "model/ObsBiasPreconditioner.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"


// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
ObsBiasCovariance::ObsBiasCovariance(const ObsSpaceQG &, const Parameters_ & params)
{
  std::array<double, ObsBias::ntypes> zz;
  zz.fill(0.0);
  if (params.covariance.value() != boost::none) {
    const ObsBiasCovarianceParameters& covparams = *params.covariance.value();
    if (covparams.stream.value() != boost::none) zz[0] = *covparams.stream.value();
    if (covparams.uwind.value() != boost::none)  zz[1] = *covparams.uwind.value();
    if (covparams.vwind.value() != boost::none)  zz[2] = *covparams.vwind.value();
    if (covparams.wspeed.value() != boost::none) zz[3] = *covparams.wspeed.value();
  }
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (std::abs(zz[jj]) > 1.0e-8) {
      variance_[jj] = zz[jj] * zz[jj];
      std::ostringstream strs;
      strs << variance_[jj];
      strn += strs.str();
    } else {
      variance_[jj] = 0.0;
      strn += "0.0";
    }
  }
  oops::Log::info() << "ObsBiasCovariance created, variances = " << strn << std::endl;
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::multiply(const ObsBiasIncrement & dxin,
                                 ObsBiasIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] * variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & dxin,
                                        ObsBiasIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] / variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::randomize(ObsBiasIncrement & dx) const {
  static util::NormalDistribution<double> dist(ObsBias::ntypes, 0.0, 1.0, 4);
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dx[jj] = dist[jj] * std::sqrt(variance_[jj]);
    } else {
      dx[jj] = 0.0;
    }
  }
}

// -----------------------------------------------------------------------------
std::unique_ptr<ObsBiasPreconditioner> ObsBiasCovariance::preconditioner() const {
  return std::make_unique<ObsBiasPreconditioner> (variance_);
}

// -----------------------------------------------------------------------------
void ObsBiasCovariance::print(std::ostream & os) const {
  os << "ObsBiasCovariance::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace qg
