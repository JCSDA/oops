/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_FFT_F_H_
#define QG_MODEL_FFT_F_H_

namespace qg {
extern "C" {
  void fft_fwd_f(const int &, const double *, double *);
  void fft_inv_f(const int &, const double *, double *);
}
}  // namespace qg

#endif  // QG_MODEL_FFT_F_H_
