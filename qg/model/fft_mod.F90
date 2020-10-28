! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for Eigen FFTs
module fft_mod

use, intrinsic :: iso_c_binding
use kinds

implicit none
private
public :: fft_fwd, fft_inv

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------
! The lines below:
! real(kind=c_double), intent(in)  :: pgrid(:)
! real(kind=c_double), intent(out) :: pfour(:)
! would not work because then fortran passes an array descriptor
! (somebody could pass a non-contiguous array using array syntax)
! instead of the address of the first element of the array. YT
subroutine fft_fwd_c(kk, pgrid, pfour) bind(C,name='fft_fwd_f')
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: kk
  real(c_double), intent(in)  :: pgrid
  real(c_double), intent(inout) :: pfour
end subroutine fft_fwd_c
!-------------------------------------------------------------------------------
subroutine fft_inv_c(kk, pfour, pgrid) bind(C,name='fft_inv_f')
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: kk
  real(c_double), intent(in)  :: pfour
  real(c_double), intent(inout) :: pgrid
end subroutine fft_inv_c
!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine fft_fwd(kk, pg, pf)
implicit none
integer, intent(in) :: kk
real(kind_real), intent(in)  :: pg(kk)
real(kind_real), intent(out) :: pf(kk+2)
real(c_double) :: zz(kk), ww(kk+2)
integer(c_int) :: nn

nn=kk
zz(:) = pg(:)
call fft_fwd_c(nn, zz(1), ww(1))
pf(:) = ww(:)

end subroutine fft_fwd

! ------------------------------------------------------------------------------

subroutine fft_inv(kk, pf, pg)
implicit none
integer, intent(in) :: kk
real(kind_real), intent(in)  :: pf(kk+2)
real(kind_real), intent(out) :: pg(kk)
real(c_double) :: zz(kk), ww(kk+2)
integer(c_int) :: nn

nn=kk
ww(:) = pf(:)
call fft_inv_c(nn, ww(1), zz(1))
pg(:) = zz(:)

end subroutine fft_inv

! ------------------------------------------------------------------------------

end module fft_mod
