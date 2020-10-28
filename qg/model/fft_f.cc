/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/fft_f.h"

#include <unsupported/Eigen/FFT>
#include <cmath>
#include <vector>

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------

void fft_fwd_f(const int & kk, const double * xx, double * ff) {
  static Eigen::FFT<double> fft;
  std::vector<double> grid(kk);
  const unsigned int size = kk/2+1;
  std::vector<std::complex<double> > coefs(size);

  for (int jj = 0; jj < kk; ++jj) grid[jj] = xx[jj];

  fft.fwd(coefs, grid);

  for (unsigned int jj = 0; jj < size; ++jj) {
    ff[2*jj]   = coefs[jj].real();
    ff[2*jj+1] = coefs[jj].imag();
  }
}

// -----------------------------------------------------------------------------

void fft_inv_f(const int & kk, const double * ff, double * xx) {
  static Eigen::FFT<double> fft;
  std::vector<double> grid(kk);
  std::vector<std::complex<double> > coefs(kk);

  const std::complex<double> zero(0.0, 0.0);
  const int size = kk/2+1;
  for (int jj = 0; jj < size; ++jj) {
    std::complex<double> zz(ff[2*jj], ff[2*jj+1]);
    coefs[jj] = zz;
  }
  for (int jj = size; jj < kk; ++jj) coefs[jj] = zero;

  fft.inv(grid, coefs);

  for (int jj = 0; jj < kk; ++jj) xx[jj] = grid[jj];
}

// -----------------------------------------------------------------------------

}  // namespace qg
