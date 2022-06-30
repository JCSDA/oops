
/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


//--


#ifndef OOPS_GENERIC_FFT_INTERFACE_F_H_
#define OOPS_GENERIC_FFT_INTERFACE_F_H_


extern "C" {

  void fft_gp2spe(float *ds_mfft, const int *size_mfft, const int *no_seq, const int *no_el_seq);
  void fft_spe2gp(float *ds_mfft, const int *size_mfft, const int *no_seq, const int *no_el_seq);

}


#endif  // OOPS_GENERIC_FFT_INTERFACE_F_H_
