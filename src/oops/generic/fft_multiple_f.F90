
!
! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!--


subroutine fft_gp2spe(ds_mfft, size_mfft, no_seq, no_el_seq) bind(c)

  use iso_c_binding

  use fft_init_f

  !--

  implicit none

  !--

  integer(c_int) :: size_mfft
  ! number of sequences (number of FFTs to be carried out)
  integer(c_int) :: no_seq
  ! number of elements per sequence
  integer(c_int) :: no_el_seq
  ! increment between the start of each sequence
  integer :: jump
  ! data structure containing the data to be processed/tranformed
  real(c_float) :: ds_mfft(size_mfft)

  jump = (no_el_seq+2)

  ! initializing the procedure for multiple FFTs
  call initialize_FFT(no_el_seq, factors, trigs)

  ! performing multiple FFTs (from gridpoint to spectral)
  call fft_gpoint2spectral_f(ds_mfft, trigs, factors, inc, jump, no_el_seq, no_seq)


end subroutine fft_gp2spe


subroutine fft_spe2gp(ds_mfft, size_mfft, no_seq, no_el_seq) bind(c)

  use iso_c_binding

  use fft_init_f

  !--

  implicit none

  !--

  integer(c_int) :: size_mfft
  ! number of sequences (number of FFTs to be carried out)
  integer(c_int) :: no_seq
  ! number of elements per sequence
  integer(c_int) :: no_el_seq
  ! increment between the start of each sequence
  integer :: jump
  ! data structure containing the data to be processed/tranformed
  real(c_float) :: ds_mfft(size_mfft)

  jump = (no_el_seq+2)

  ! initializing the procedure for multiple FFTs
  call initialize_FFT(no_el_seq, factors, trigs)

  ! performing multiple FFTs (from spectral to gridpoint)
  call fft_spectral2gpoint_f(ds_mfft, trigs, factors, inc, jump, no_el_seq, no_seq)


end subroutine fft_spe2gp
