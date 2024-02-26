/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ErrorCovarianceL95.h"

#include <algorithm>
#include <cmath>

#include "eckit/config/Configuration.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "oops/util/Logger.h"

namespace oops {
class Variables;
}

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ErrorCovarianceL95::ErrorCovarianceL95(const Resolution & geom, const oops::Variables &,
                                       const eckit::Configuration & config,
                                       const StateL95 &, const StateL95 &) :
  time_(util::DateTime(config.getString("date"))),
  sigmab_(config.getDouble("standard_deviation")),
  rscale_(1.0/config.getDouble("length_scale"))
{
// Gaussian structure function
  resol_ = geom.npoints();
  size_ = resol_/2+1;
  std::vector<double> structFct(resol_);
  for (unsigned int jj = 0; jj < resol_; ++jj) {
    double zz = rscale_ * std::min(jj, resol_-jj);
    structFct[jj] = std::exp(-0.5*zz*zz);
  }
// Go to Fourier space
  std::vector<std::complex<double> > coefs(size_);
  fft_.fwd(coefs, structFct);
// Save Fourier coefficients
  bcoefs_.resize(size_);
  for (unsigned int jj = 0; jj < size_; ++jj) {
    bcoefs_[jj] = std::real(coefs[jj]);
  }
}
// -----------------------------------------------------------------------------
ErrorCovarianceL95::~ErrorCovarianceL95() {}
// -----------------------------------------------------------------------------
void ErrorCovarianceL95::multiply(const IncrementL95 & dxin,
                                  IncrementL95 & dxout) const {
  std::vector<std::complex<double> > coefs(resol_);
  fft_.fwd(coefs, dxin.asVector());
  for (unsigned int jj = 0; jj < size_; ++jj) {
    coefs[jj] *= bcoefs_[jj];
  }
  fft_.inv(dxout.asVector(), coefs);
  double var = sigmab_ * sigmab_;
  dxout *= var;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceL95::inverseMultiply(const IncrementL95 & dxin,
                                         IncrementL95 & dxout) const {
  std::vector<std::complex<double> > coefs(resol_);
  fft_.fwd(coefs, dxin.asVector());
  for (unsigned int jj = 0; jj < size_; ++jj) {
    coefs[jj] /= bcoefs_[jj];
  }
  fft_.inv(dxout.asVector(), coefs);
  double vari = 1.0 / (sigmab_ * sigmab_);
  dxout *= vari;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceL95::randomize(IncrementL95 & dx) const {
  dx.random();
  std::vector<std::complex<double> > coefs(resol_);
  fft_.fwd(coefs, dx.asVector());
  for (unsigned int jj = 0; jj < size_; ++jj) {
    coefs[jj] *= std::sqrt(bcoefs_[jj]);
  }
  fft_.inv(dx.asVector(), coefs);
  dx *= sigmab_;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceL95::print(std::ostream & os) const {
  os << "ErrorCovarianceL95: time = " << time_ << ", std dev = " << sigmab_
     << ", length scale = " << 1.0/rscale_;
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
