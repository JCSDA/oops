! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Interfaces to be called from C++ for Fortran computation of eigenvalues
!> and eigenvectors of symetric tri-diagonal matrix.

! ------------------------------------------------------------------------------

subroutine tridiagev(nn, diag, subd, eval, evec) bind(c,name='FtnTriDiagSpectrum')

use iso_c_binding
use kinds

implicit none
integer(c_int), intent(in) :: nn              !< Size of matrix
real(c_double), intent(in) :: diag(nn)        !< Diagonal elements
real(c_double), intent(in) :: subd(nn-1)      !< Sub-diagonal elements
real(c_double), intent(inout) :: eval(nn)     !< Eigenvalues
real(c_double), intent(inout) :: evec(nn,nn)  !< Eigenvectors

real(kind=kind_real) :: zdiag(nn)
real(kind=kind_real) :: zsubd(nn-1)
real(kind=kind_real) :: zvecs(nn, nn)
real(kind=kind_real) :: zwork(2*nn)
integer :: ii, info

ii = nn
zdiag(:) = diag(:)
zsubd(:) = subd(:)

call dsteqr('I', ii, zdiag, zsubd, zvecs, ii, zwork, info)

if (info/=0) call abor1_ftn("tridiagev: Error in dsteqr")

eval(:) = zdiag(:)
evec(:,:) = zvecs(:,:)

end subroutine tridiagev

! ------------------------------------------------------------------------------

