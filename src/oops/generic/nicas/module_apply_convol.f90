!----------------------------------------------------------------------
! Module: module_apply_convol.f90
!> Purpose: convolution routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_convol

use module_apply_com, only: alpha_com_AB,alpha_com_BA,alpha_com_CA,alpha_copy_AC,alpha_copy_BC
use module_namelist, only: nam
use omp_lib
use tools_missing, only: msr
use type_fields, only: alphatype
use type_linop, only: apply_linop_sym
use type_mpl, only: mpl
use type_sdata, only: sdatatype,sdatampitype
implicit none

private
public :: convol

contains

!----------------------------------------------------------------------
! Subroutine: convol
!> Purpose: convolution
!----------------------------------------------------------------------
subroutine convol(sdata,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                 !< Sampling data
type(alphatype),intent(inout) :: alpha(sdata%nproc) !< Subgrid variable

! Local variables
integer :: iproc

if (nam%nproc==0) then
   ! Apply global convolution
   call convol_global(sdata,alpha(1))
elseif (nam%nproc>0) then
   ! Communication
   if (sdata%mpicom==1) then
      ! Copy zone B into zone C
      do iproc=1,sdata%nproc
         call alpha_copy_BC(sdata%mpi(iproc),alpha(iproc))
      end do
   elseif (sdata%mpicom==2) then
      ! Halo reduction from zone B to zone A
      call alpha_com_BA(sdata,alpha)

      ! Copy zone A into zone C
      do iproc=1,sdata%nproc
         call alpha_copy_AC(sdata%mpi(iproc),alpha(iproc))
      end do
   end if

   ! Apply local convolution
   do iproc=1,sdata%nproc
      call convol_local(sdata%mpi(iproc),alpha(iproc))
   end do
end if

! Communication
if (nam%nproc>0) then
   ! Halo reduction from zone C to zone A
   call alpha_com_CA(sdata,alpha)

   ! Halo extension from zone A to zone B
   call alpha_com_AB(sdata,alpha)
end if

end subroutine convol

!----------------------------------------------------------------------
! Subroutine: convol_global
!> Purpose: convolution, global
!----------------------------------------------------------------------
subroutine convol_global(sdata,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata    !< Sampling data
type(alphatype),intent(inout) :: alpha !< Subgrid variable

! Apply linear operator, symmetric
call apply_linop_sym(sdata%c,alpha%val)

end subroutine convol_global

!----------------------------------------------------------------------
! Subroutine: convol_local
!> Purpose: convolution, local
!----------------------------------------------------------------------
subroutine convol_local(sdatampi,alpha)

implicit none

! Passed variables
type(sdatampitype),intent(in) :: sdatampi !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Apply linear operator, symmetric
call apply_linop_sym(sdatampi%c,alpha%valc)

end subroutine convol_local

end module module_apply_convol
