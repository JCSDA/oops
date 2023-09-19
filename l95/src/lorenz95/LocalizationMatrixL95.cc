/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/LocalizationMatrixL95.h"

#include <unsupported/Eigen/FFT>
#include <algorithm>
#include <cmath>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "lorenz95/IncrementL95.h"
#include "lorenz95/Resolution.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
LocalizationMatrixL95::LocalizationMatrixL95(const Resolution & resol,
                                             const oops::Variables & vars,
                                             const eckit::Configuration & config)
  : resol_(resol.npoints()),
    rscale_(1.0/config.getDouble("length_scale"))
{
// Gaussian structure function
  unsigned int size = resol_/2+1;
  std::vector<double> locfct(resol_);
  for (unsigned int jj = 0; jj < resol_; ++jj) {
    double zz = rscale_ * std::min(jj, resol_-jj);
    locfct[jj] = std::exp(-0.5*zz*zz);
  }
// Go to Fourier space
  Eigen::FFT<double> fft;
  std::vector<std::complex<double> > four(size);
  fft.fwd(four, locfct);
// Save Fourier coefficients
  coefs_.resize(size);
  for (unsigned int jj = 0; jj < size; ++jj) {
    coefs_[jj] = std::real(four[jj]);
  }
}
// -----------------------------------------------------------------------------
void LocalizationMatrixL95::randomize(IncrementL95 & dx) const {
  dx.random();
  unsigned int size = resol_/2+1;
  Eigen::FFT<double> fft;
  std::vector<std::complex<double> > four(size);
  fft.fwd(four, dx.asVector());
  for (unsigned int jj = 0; jj < size; ++jj) {
    four[jj] *= std::sqrt(coefs_[jj]);
  }
  fft.inv(dx.asVector(), four);
}
// -----------------------------------------------------------------------------
void LocalizationMatrixL95::multiply(IncrementL95 & dx) const {
  unsigned int size = resol_/2+1;
  Eigen::FFT<double> fft;
  std::vector<std::complex<double> > four(size);
  fft.fwd(four, dx.asVector());
  for (unsigned int jj = 0; jj < size; ++jj) {
    four[jj] *= coefs_[jj];
  }
  fft.inv(dx.asVector(), four);
}
// -----------------------------------------------------------------------------
void LocalizationMatrixL95::print(std::ostream & os) const {
  os << "Localization with Gaussian, lengthscale = " << 1.0/rscale_;
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
