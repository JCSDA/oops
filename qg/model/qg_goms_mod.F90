! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling interpolated (to obs locations) model variables

module qg_goms_mod

use iso_c_binding
use qg_locs_mod
use qg_vars_mod
use kinds

implicit none
private
public :: qg_goms, gom_setup
public :: qg_goms_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold interpolated fields required by the obs operators
type :: qg_goms
  integer :: nobs
  integer :: nvar
  integer :: used
  integer, allocatable :: indx(:)
  real(kind=kind_real), allocatable :: values(:,:)
  character(len=1), allocatable :: variables(:)
  logical :: lalloc
end type qg_goms

#define LISTED_TYPE qg_goms

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_goms_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_qg_gom_setup(c_key_self, c_key_locs, c_vars) bind(c,name='qg_gom_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(inout) :: c_key_locs
integer(c_int), dimension(*), intent(in) :: c_vars     !< List of variables

type(qg_goms), pointer :: self
type(qg_locs), pointer :: locs
type(qg_vars) :: vars

call qg_goms_registry%init()
call qg_goms_registry%add(c_key_self)
call qg_goms_registry%get(c_key_self, self)
call qg_locs_registry%get(c_key_locs, locs)
call qg_vars_create(vars, c_vars)

call gom_setup(self, vars, locs%indx)

end subroutine c_qg_gom_setup

! ------------------------------------------------------------------------------

subroutine c_qg_gom_create(c_key_self) bind(c,name='qg_gom_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(qg_goms), pointer :: self
call qg_goms_registry%init()
call qg_goms_registry%add(c_key_self)
call qg_goms_registry%get(c_key_self, self)

self%lalloc = .false.

end subroutine c_qg_gom_create

! ------------------------------------------------------------------------------

subroutine gom_setup(self, vars, kobs)
implicit none
type(qg_goms), intent(inout) :: self
type(qg_vars), intent(in) :: vars
integer, intent(in) :: kobs(:)

self%nobs=size(kobs)
self%nvar=vars%nv
self%used=0

allocate(self%indx(self%nobs))
self%indx(:)=kobs(:)

allocate(self%variables(self%nvar))
self%variables(:)=vars%fldnames(:)

allocate(self%values(self%nvar,self%nobs))

self%lalloc = .true.

end subroutine gom_setup

! ------------------------------------------------------------------------------

subroutine c_qg_gom_delete(c_key_self) bind(c,name='qg_gom_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(qg_goms), pointer :: self

call qg_goms_registry%get(c_key_self, self)
if (self%lalloc) then
  deallocate(self%values)
  deallocate(self%indx)
  deallocate(self%variables)
endif
call qg_goms_registry%remove(c_key_self)

end subroutine c_qg_gom_delete

! ------------------------------------------------------------------------------

subroutine c_qg_gom_assign(c_key_self, c_key_other) bind(c,name='qg_gom_assign_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(qg_goms), pointer :: self
type(qg_goms), pointer :: other
integer :: jo, jv

call qg_goms_registry%get(c_key_self, self)
call qg_goms_registry%get(c_key_other, other)

self%nobs = other%nobs
self%nvar = other%nvar

if (.not. self%lalloc) then
   allocate(self%variables(self%nvar))
   allocate(self%values(self%nvar,self%nobs))
   allocate(self%indx(self%nobs))
   self%lalloc = .true.
endif
   
self%variables = other%variables
self%values = other%values
self%indx = other%indx
self%used = other%used

end subroutine c_qg_gom_assign

! ------------------------------------------------------------------------------

subroutine c_qg_gom_zero(c_key_self) bind(c,name='qg_gom_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(qg_goms), pointer :: self
call qg_goms_registry%get(c_key_self, self)
self%values(:,:)=0.0_kind_real
end subroutine c_qg_gom_zero

! ------------------------------------------------------------------------------

subroutine c_qg_gom_abs(c_key_self) bind(c,name='qg_gom_abs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(qg_goms), pointer :: self
call qg_goms_registry%get(c_key_self, self)
self%values(:,:)=abs(self%values(:,:))
end subroutine c_qg_gom_abs

! ------------------------------------------------------------------------------

subroutine c_qg_gom_random(c_key_self) bind(c,name='qg_gom_random_f90')
use random_vectors_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(qg_goms), pointer :: self
call qg_goms_registry%get(c_key_self, self)
call random_vector(self%values(:,:))
end subroutine c_qg_gom_random

! ------------------------------------------------------------------------------

subroutine c_qg_gom_mult(c_key_self, zz) bind(c,name='qg_gom_mult_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: zz
type(qg_goms), pointer :: self
integer :: jo, jv

call qg_goms_registry%get(c_key_self, self)
do jo=1,self%nobs
  do jv=1,self%nvar
    self%values(jv,jo) = zz * self%values(jv,jo)
  enddo
enddo

end subroutine c_qg_gom_mult

! ------------------------------------------------------------------------------

subroutine c_qg_gom_add(c_key_self, c_key_other) bind(c,name='qg_gom_add_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(qg_goms), pointer :: self
type(qg_goms), pointer :: other
integer :: jo, jv

call qg_goms_registry%get(c_key_self, self)
call qg_goms_registry%get(c_key_other, other)
do jo=1,self%nobs
  do jv=1,self%nvar
    self%values(jv,jo) = self%values(jv,jo) + other%values(jv,jo)
  enddo
enddo

end subroutine c_qg_gom_add

! ------------------------------------------------------------------------------

subroutine c_qg_gom_diff(c_key_self, c_key_other) bind(c,name='qg_gom_diff_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(qg_goms), pointer :: self
type(qg_goms), pointer :: other
integer :: jo, jv

call qg_goms_registry%get(c_key_self, self)
call qg_goms_registry%get(c_key_other, other)
do jo=1,self%nobs
  write(*,*) "qg_gom_diff u: ",jo,self%values(3,jo),other%values(3,jo)
  write(*,*) "qg_gom_diff v: ",jo,self%values(4,jo),other%values(4,jo)
  do jv=1,self%nvar
    self%values(jv,jo) = self%values(jv,jo) - other%values(jv,jo)
  enddo
enddo

end subroutine c_qg_gom_diff

! ------------------------------------------------------------------------------
!> QG GeoVaLs division Operator

subroutine c_qg_gom_divide(c_key_self, c_key_other) bind(c,name='qg_gom_divide_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(qg_goms), pointer :: self
type(qg_goms), pointer :: other
real(kind_real) :: tol
integer :: jloc, jvar, ii

call qg_goms_registry%get(c_key_self, self)
call qg_goms_registry%get(c_key_other, other)

tol = epsilon(tol)
ii=0

do jvar = 1, self%nvar
  write(*,*)'qg_gom_divide variable =', trim(self%variables(jvar))
  do jloc = 1, self%nobs
    if (abs(other%values(jvar,jloc)) > tol) then
      write(*,*)'qg_gom_divide self, other = ',self%values(jvar,jloc), other%values(jvar,jloc), &
                                             & self%values(jvar,jloc) / other%values(jvar,jloc)
      self%values(jvar,jloc) = self%values(jvar,jloc) / other%values(jvar,jloc)
    else
      write(*,*)'qg_gom_divide self, other = ',self%values(jvar,jloc), other%values(jvar,jloc)
      self%values(jvar,jloc) = 0.0
      ii=ii+1
    endif
  enddo
enddo

end subroutine c_qg_gom_divide

! ------------------------------------------------------------------------------

subroutine c_qg_gom_rms(c_key_self, rms) bind(c,name='qg_gom_rms_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: rms
type(qg_goms), pointer :: self
integer :: jo, jv

call qg_goms_registry%get(c_key_self, self)

rms=0.0_kind_real
do jo=1,self%nobs
  do jv=1,self%nvar
write(*,*) "MSM: ",jv,jo,self%values(jv,jo)
    rms=rms+self%values(jv,jo)**2
  enddo
enddo

rms = sqrt(rms/(self%nobs*self%nvar))

end subroutine c_qg_gom_rms

! ------------------------------------------------------------------------------

subroutine c_qg_gom_dotprod(c_key_self, c_key_other, prod) bind(c,name='qg_gom_dotprod_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other
real(c_double), intent(inout) :: prod
type(qg_goms), pointer :: self, other
integer :: jo, jv

call qg_goms_registry%get(c_key_self, self)
call qg_goms_registry%get(c_key_other, other)
prod=0.0_kind_real
do jo=1,self%nobs
  do jv=1,self%nvar
    prod=prod+self%values(jv,jo)*other%values(jv,jo)
  enddo
enddo

end subroutine c_qg_gom_dotprod

! ------------------------------------------------------------------------------

subroutine c_qg_gom_minmaxavg(c_key_self, kobs, pmin, pmax, prms) bind(c,name='qg_gom_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
real(c_double), intent(inout) :: pmin, pmax, prms
type(qg_goms), pointer :: self

call qg_goms_registry%get(c_key_self, self)

kobs = self%nobs
pmin=minval(self%values(:,:))
pmax=maxval(self%values(:,:))
prms=sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

end subroutine c_qg_gom_minmaxavg

! ------------------------------------------------------------------------------

subroutine c_qg_gom_maxloc(c_key_self, mxval, iloc, ivar) bind(c,name='qg_gom_maxloc_f90')
implicit none
integer(c_int), intent(in) :: c_key_self !< key to GeoVals object
real(c_double), intent(inout) :: mxval   !< maximum value
integer(c_int), intent(inout) :: iloc    !< location of maximum value
integer(c_int), intent(inout) :: ivar    !< variable with maximum value

type(qg_goms), pointer :: self
integer :: mxloc(2)

call qg_goms_registry%get(c_key_self, self)

mxval=maxval(self%values(:,:))
mxloc=maxloc(self%values(:,:))

ivar = mxloc(1)
iloc = mxloc(2)

end subroutine c_qg_gom_maxloc

! ------------------------------------------------------------------------------

subroutine qg_gom_read_file_c(c_key_self, c_conf) bind(c,name='qg_gom_read_file_f90')
use config_mod
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
type(qg_goms), pointer :: self

integer, parameter :: iunit=10
integer, parameter :: max_string_length=250 ! Yuk!
character(len=max_string_length) :: filename, record
character(len=4)  :: cnx
character(len=17) :: fmtn
character(len=11) :: fmt1='(X,ES24.16)'
integer :: jj, jo, jv

call qg_goms_registry%get(c_key_self, self)
if (self%lalloc) call abor1_ftn("qg_gom_read_file gom alredy allocated")

filename = config_get_string(c_conf,len(filename),"filename")
write(record,*)'qg_gom_read_file: opening '//trim(filename)
call fckit_log%info(record)
open(unit=iunit, file=trim(filename), form='formatted', action='read')

read(iunit,*) self%nobs, self%nvar, self%used
allocate(self%indx(self%nobs))
allocate(self%variables(self%nvar))
allocate(self%values(self%nvar,self%nobs))

read(iunit,*) self%indx(:)
do jv=1,self%nvar
  read(iunit,*) self%variables(jv)
enddo

if (self%nvar>9999)  call abor1_ftn("Format too small")
write(cnx,'(I4)')self%nvar
fmtn='('//trim(cnx)//fmt1//')'

do jo=1,self%nobs
  read(iunit,fmtn) (self%values(jj,jo), jj=1,self%nvar)
enddo

close(iunit)
self%lalloc = .true.

end subroutine qg_gom_read_file_c

! ------------------------------------------------------------------------------

subroutine qg_gom_write_file_c(c_key_self, c_conf) bind(c,name='qg_gom_write_file_f90')
use config_mod
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_conf
type(qg_goms), pointer :: self

integer, parameter :: iunit=10
integer, parameter :: max_string_length=250 ! Yuk!
character(len=max_string_length) :: filename, record
character(len=4)  :: cnx
character(len=17) :: fmtn
character(len=11) :: fmt1='(X,ES24.16)'
integer :: jj, jo, jv

call qg_goms_registry%get(c_key_self, self)
if (.not.self%lalloc) call abor1_ftn("qg_gom_write_file gom not allocated")

filename = config_get_string(c_conf,len(filename),"filename")
write(record,*)'qg_gom_write_file: opening '//trim(filename)
call fckit_log%info(record)
open(unit=iunit, file=trim(filename), form='formatted', action='write')

write(iunit,*) self%nobs, self%nvar, self%used
write(iunit,*) self%indx(:)
do jv=1,self%nvar
  write(iunit,*) self%variables(jv)
enddo

if (self%nvar>9999) call abor1_ftn("Format too small")
write(cnx,'(I4)')self%nvar
fmtn='('//trim(cnx)//fmt1//')'

do jo=1,self%nobs
  write(iunit,fmtn) (self%values(jj,jo), jj=1,self%nvar)
enddo

close(iunit)

end subroutine qg_gom_write_file_c

! ------------------------------------------------------------------------------
!> Initialize a QG GeoVaLs object based on an analytic state
!!
!! \details **qg_gom_analytic_init_c()** creates a QG GeoVaLs object and fills in 
!! values based on one of several analytic solutions.  For a description of the
!! solutions that are currently implemented, see the analytic_init() routine in
!! qg_fields.F90.  This initialization is intended to be used with the
!! **TestStateInterpolation()** test.
!!
!! \date April, 2018: Created by M. Miesch (JCSDA)
!!
!! \sa qg::GomQG analytic_init()
!!
subroutine qg_gom_analytic_init_c(c_key_self, c_key_locs, c_conf) bind(c,name='qg_gom_analytic_init_f90')
use config_mod
use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, test1_advection_hadley  
use qg_constants, only : ubar 
use fckit_log_module, only : fckit_log

implicit none
integer(c_int), intent(in) :: c_key_self           !< key to GeoVaLs object we are creating
integer(c_int), intent(in) :: c_key_locs           !< key to the F90 Locations object
type(c_ptr), intent(in)    :: c_conf               !< "StateGenerate" configuration

type(qg_goms), pointer :: self  ! GeoVals object we are creating         
type(qg_locs), pointer :: locs

character(len=30) :: ic
character(len=1024) :: buf
real(kind_real) :: d1, d2, height(2), klon, klat, kz
real(kind_real) :: pi = acos(-1.0_kind_real)
real(kind_real) :: p0,u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4
integer :: iloc, ivar

! Get F90 Locations object
call qg_locs_registry%get(c_key_locs, locs)

! The GeoVaLs object should already be allocated and defined.
! We just have to replace the values with the appropriate analytic values

call qg_goms_registry%get(c_key_self, self)

if (.not. self%lalloc) &
   call abor1_ftn("qg_gom_analytic init: gom not allocated")

ic = config_get_string(c_conf,len(ic),"analytic_init")
write(buf,*)'qg_gom_analytic_init: ic '//trim(ic)
call fckit_log%info(buf)

! The height is stored in the locations object as the layer number because
! this is what the interp_tl() function currently requires.  Here we have
! to convert this to meters in order to use the analytic formulae.
! So, as in the fields analytic_init() routine, set the height to be the
! middle of the corresponding layer.
d1 = config_get_real(c_conf,"top_layer_depth")
d2 = config_get_real(c_conf,"bottom_layer_depth")
height(2) = 0.5_kind_real*d2
height(1) = d2 + 0.5_kind_real*d1

! Now loop over locations

do iloc = 1, locs%nloc

   ! First get the locations and convert to latitude (rad), longiude (rad) and height (meters)

   klon = (locs%xyz(1,iloc)*2.0_kind_real - 1.0_kind_real)*pi
   klat = (locs%xyz(2,iloc) - 0.5_kind_real)*pi
   kz   = height(nint(locs%xyz(3,iloc)))

   init_option: select case (ic)

   case ("dcmip-test-1-1")
      
      call test1_advection_deformation(klon,klat,p0,kz,1,u0,v0,w0,&
                                       t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4)

      ! Make Nondimensional
      u0 = u0 / ubar
      v0 = v0 / ubar
      
   case ("dcmip-test-1-2")

      call test1_advection_hadley(klon,klat,p0,kz,1,u0,v0,w0,&
                                  t0,phis0,ps0,rho0,hum0,q1)

      ! Make Nondimensional
      u0 = u0 / ubar
      v0 = v0 / ubar

   case default

      write(buf,*) "qg_gom_analytic_init: unknown init = " // ic
      call fckit_log%error(buf)
      call abor1_ftn(buf)

   end select init_option

   do ivar = 1,self%nvar

      variable: select case(trim(self%variables(ivar)))
      case ("u")
         self%values(ivar,iloc) = u0
      case ("v")
         self%values(ivar,iloc) = v0
      case default
         self%values(ivar,iloc) = 0.0_kind_real
      end select variable

   enddo
      
enddo

end subroutine qg_gom_analytic_init_c

! ------------------------------------------------------------------------------

end module qg_goms_mod
