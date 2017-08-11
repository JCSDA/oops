!----------------------------------------------------------------------
! Module: type_randgen
!> Purpose: random numbers generator derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_randgen

use iso_c_binding, only: c_ptr,c_int,c_double
use module_namelist, only: nam
use tools_kinds, only: kind_real
use type_mpl, only: mpl,mpl_bcast
implicit none

type randgentype
   type(c_ptr) :: ptr !< Pointer to the C++ class
end type randgentype

! C++ interface
interface
   function create_randgen_c(default_seed) bind(C,name="create_randgen")
   use iso_c_binding
   implicit none
   type(c_ptr) :: create_randgen_c
   integer(c_int),value :: default_seed
   end function create_randgen_c
end interface
interface
   subroutine delete_randgen_c(randgen) bind(C,name="delete_randgen")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: randgen
   end subroutine delete_randgen_c
end interface
interface
   subroutine rand_integer_c(randgen,binf,bsup,ir) bind(C,name="rand_integer")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: randgen
   integer(c_int),value :: binf
   integer(c_int),value :: bsup
   integer(c_int) :: ir
   end subroutine rand_integer_c
end interface
interface
   subroutine initialize_sampling_c(randgen,n,lon,lat,mask,L,ntry,nrep,ns,nfor,ifor,ihor) bind(C,name="initialize_sampling")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: randgen
   integer(c_int),value :: n
   real(c_double) :: lon(*)
   real(c_double) :: lat(*)
   integer(c_int) :: mask(*)
   real(c_double) :: L(*)
   integer(c_int),value :: ntry
   integer(c_int),value :: nrep
   integer(c_int),value :: ns
   integer(c_int),value :: nfor
   integer(c_int) :: ifor(*)
   integer(c_int) :: ihor(*)
   end subroutine initialize_sampling_c
end interface

private
public :: randgentype,create_randgen,delete_randgen,rand_integer,initialize_sampling

contains

!----------------------------------------------------------------------
! Subroutine: create_randgen
!> Purpose: create a random number generator
!----------------------------------------------------------------------
function create_randgen()

implicit none

! Passed variables
type(randgentype) :: create_randgen !< Random number generator object

! Local variable
integer :: default_seed

if (mpl%main) then
   ! Set default seed key to integer
   if (nam%sam_default_seed) then
      default_seed = 1
   else
      default_seed = 0
   end if

   ! Call C++ function
   create_randgen%ptr = create_randgen_c(default_seed)
end if

end function create_randgen

!----------------------------------------------------------------------
! Subroutine: delete_randgen
!> Purpose: delete a random number generator
!----------------------------------------------------------------------
subroutine delete_randgen(this)

implicit none

! Passed variables
type(randgentype),intent(inout) :: this !< Random number generator object

if (mpl%main) then
   ! Call C++ function
   call delete_randgen_c(this%ptr)
end if

end subroutine delete_randgen

!----------------------------------------------------------------------
! Subroutine: rand_integer
!> Purpose: generate a random integer
!----------------------------------------------------------------------
subroutine rand_integer(this,binf,bsup,ir)

implicit none

! Passed variables
class(randgentype),intent(in) :: this !< Random number generator object
integer,intent(in) :: binf            !< Lower bound
integer,intent(in) :: bsup            !< Upper bound
integer,intent(out) :: ir             !< Random integer

if (mpl%main) then
   ! Call C++ function
   call rand_integer_c(this%ptr,binf,bsup,ir)
end if

! Broadcast
call mpl_bcast(ir,mpl%ioproc)

end subroutine rand_integer

!----------------------------------------------------------------------
! Subroutine: initialize_sampling
!> Purpose: intialize sampling
!----------------------------------------------------------------------
subroutine initialize_sampling(this,n,lon,lat,mask,L,ntry,nrep,ns,nfor,ifor,ihor)

implicit none

! Passed variables
class(randgentype),intent(in) :: this !< Random number generator object
integer,intent(in) :: n               !< Number of grid points
real(kind_real),intent(in) :: lon(n)  !< Grid points longitude
real(kind_real),intent(in) :: lat(n)  !< Grid points latitude
integer,intent(in) :: mask(n)         !< Grid mask
real(kind_real),intent(in) :: L(n)    !< Grid length-scale
integer,intent(in) :: ntry            !< Number of tries
integer,intent(in) :: nrep            !< Number of replacements
integer,intent(in) :: ns              !< Number of samplings points
integer,intent(in) :: nfor            !< Number of boundary points
integer,intent(in) :: ifor(nfor)      !< Boudary points
integer,intent(out) :: ihor(ns)       !< Horizontal sampling index

if (mpl%main) then
   ! Call C++ function
   call initialize_sampling_c(this%ptr,n,lon,lat,mask,L,ntry,nrep,ns,nfor,ifor,ihor)
end if

! Broadcast
call mpl_bcast(ihor,mpl%ioproc)

end subroutine initialize_sampling

end module type_randgen
