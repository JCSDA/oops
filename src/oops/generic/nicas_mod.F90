! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic nicas localization

module nicas_mod

use iso_c_binding
use kinds
use config_mod
use unstructured_grid_mod
use module_nicas, only: nicas_driver
use model_oops, only: model_oops_coord
use module_namelist, only: nam,namcheck
use tools_display, only: listing_setup
use type_esmf, only: esmf_start,esmf_end
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype
use fckit_log_module, only : log

implicit none
private
public nicas, create_nicas, delete_nicas, nicas_multiply

! ------------------------------------------------------------------------------

!>  Derived type containing the data

type nicas
  type(ndataloctype) :: ndataloc
end type nicas

! ------------------------------------------------------------------------------

#define LISTED_TYPE nicas

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: nicas_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------
!  C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_nicas_c(key, c_conf, cndims, cdims, cnh, clats, clons, cnv, clevs, carea, cmask) bind(c, name='create_nicas_f90')
implicit none
integer(c_int), intent(inout) :: key
type(c_ptr), intent(in) :: c_conf
integer(c_int), intent(in) :: cndims, cnh, cnv
real(c_double), intent(in) :: clats(cnh), clons(cnh), clevs(cnv), carea(cnv)
integer(c_int), intent(in) :: cdims(cndims), cmask(cnh*cnv)
type(nicas), pointer :: self
real(kind=kind_real) :: lats(cnh), lons(cnh), levs(cnv), area(cnv)
integer :: dims(cndims), mask(cnh*cnv)
call nicas_registry%init()
call nicas_registry%add(key)
call nicas_registry%get(key,self)
dims(:)=cdims(:)
lats(:)=clats(:)
lons(:)=clons(:)
levs(:)=clevs(:)
area(:)=carea(:)
mask(:)=cmask(:)
call create_nicas(self, c_conf, dims, lats, lons, levs, area, mask)
end subroutine create_nicas_c

! ------------------------------------------------------------------------------

subroutine delete_nicas_c(key) bind(c, name='delete_nicas_f90')
implicit none
integer(c_int), intent(inout) :: key
type(nicas), pointer :: self
call nicas_registry%get(key,self)
call delete_nicas(self)
call nicas_registry%remove(key)
end subroutine delete_nicas_c

! ------------------------------------------------------------------------------

subroutine nicas_multiply_c(key, idx) bind(c, name='nicas_multiply_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx
type(nicas), pointer :: self
type(unstructured_grid), pointer :: udx
call nicas_registry%get(key,self)
call unstructured_grid_registry%get(idx, udx)
call nicas_multiply(self, udx)
end subroutine nicas_multiply_c

! ------------------------------------------------------------------------------
!  End C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_nicas(self, c_conf, dims, lats, lons, levs, area, mask)
implicit none
type(nicas), intent(inout) :: self
type(c_ptr), intent(in) :: c_conf
real(kind=kind_real), intent(in) :: lats(:), lons(:), levs(:), area(:)
integer, intent(in) :: dims(:),mask(:)
integer :: ndims,nc0,nlev
character(len=4) :: myprocchar,nprocchar,nthreadchar
type(ndatatype) :: ndata

! NICAS setup
call log%info("NICAS setup")

! Read JSON
call log%info("Read JSON")
call nicas_read_conf(c_conf, self)

! Setup display
call log%info("Listing setup")
call listing_setup(nam%colorlog)

! Check namelist parameters
call namcheck

! Write parallel setup
write(nprocchar,'(i4)') mpl%nproc
write(nthreadchar,'(i4)') mpl%nthread
call log%info("Parallel setup: "//nprocchar//" MPI tasks and "//nthreadchar//" OpenMP threads")

! Initialize ESMF
call log%info("Initialize ESMF")
call esmf_start

! Initialize coordinates
call log%info("Initialize coordinates")
call model_oops_coord(dims, lats, lons, levs, area, mask, ndata)

! Call driver
call nicas_driver(ndata,self%ndataloc)

! Close listing files
if ((mpl%main.and..not.nam%colorlog).or..not.mpl%main) close(unit=mpl%unit)

write(*,*) 'NICAS setup done'

end subroutine create_nicas

!-------------------------------------------------------------------------------

subroutine nicas_read_conf(c_conf,self)
implicit none
type(c_ptr), intent(in) :: c_conf
type(nicas), intent(inout) :: self
integer :: il,idir
character(len=3) :: ilchar,idirchar

! general_param
nam%datadir = config_get_string(c_conf,1024,"datadir")
nam%prefix = config_get_string(c_conf,1024,"prefix")
nam%colorlog = integer_to_logical(config_get_int(c_conf,"colorlog"))
nam%model = config_get_string(c_conf,1024,"model")
nam%nl = config_get_int(c_conf,"nl")
do il=1,nam%nl
   write(ilchar,'(i3)') il
   nam%levs(il) = config_get_int(c_conf,"levs("//trim(adjustl(ilchar))//")")
end do
nam%new_param = integer_to_logical(config_get_int(c_conf,"new_param"))
nam%new_mpi = integer_to_logical(config_get_int(c_conf,"new_mpi"))
nam%check_adjoints = integer_to_logical(config_get_int(c_conf,"check_adjoints"))
nam%check_pos_def = integer_to_logical(config_get_int(c_conf,"check_pos_def"))
nam%check_mpi = integer_to_logical(config_get_int(c_conf,"check_mpi"))
nam%check_dirac = integer_to_logical(config_get_int(c_conf,"check_dirac"))
nam%check_perf = integer_to_logical(config_get_int(c_conf,"check_perf"))
nam%ndir = config_get_int(c_conf,"ndir")
nam%dirlev = config_get_int(c_conf,"dirlev")
do idir=1,nam%ndir
   write(idirchar,'(i3)') idir
   nam%dirlon(idir) = config_get_real(c_conf,"dirlon("//trim(adjustl(idirchar))//")")
   nam%dirlat(idir) = config_get_real(c_conf,"dirlat("//trim(adjustl(idirchar))//")")
end do

! sampling_param
nam%sam_default_seed = integer_to_logical(config_get_int(c_conf,"sam_default_seed"))
nam%mask_check = integer_to_logical(config_get_int(c_conf,"mask_check"))
nam%ntry = config_get_int(c_conf,"ntry")
nam%nrep = config_get_int(c_conf,"nrep")
nam%logpres = integer_to_logical(config_get_int(c_conf,"logpres"))

! nicas_param
nam%lsqrt = integer_to_logical(config_get_int(c_conf,"lsqrt"))
nam%Lbh_file = config_get_string(c_conf,1024,"Lbh_file")
do il=1,nam%nl
   write(ilchar,'(i3)') il
   nam%Lbh(il) = config_get_real(c_conf,"Lbh("//trim(adjustl(ilchar))//")")
end do
nam%Lbv_file = config_get_string(c_conf,1024,"Lbv_file")
do il=1,nam%nl
   write(ilchar,'(i3)') il
   nam%Lbv(il) = config_get_real(c_conf,"Lbv("//trim(adjustl(ilchar))//")")
end do
nam%resol = config_get_real(c_conf,"resol")
nam%network = integer_to_logical(config_get_int(c_conf,"network"))
nam%nproc = config_get_int(c_conf,"nproc")
nam%mpicom = config_get_int(c_conf,"mpicom")

end subroutine nicas_read_conf

!-------------------------------------------------------------------------------

logical function integer_to_logical(i)
implicit none
integer,intent(in) :: i

if (i==0) then
   integer_to_logical = .false.
elseif (i==1) then
   integer_to_logical = .true.
else
   call abor1_ftn('wrong integer in integer_to_logical')
end if

end function integer_to_logical

!-------------------------------------------------------------------------------

subroutine delete_nicas(self)
implicit none
type(nicas), intent(inout) :: self

end subroutine delete_nicas

!-------------------------------------------------------------------------------

subroutine nicas_multiply(self,dx)
implicit none
type(nicas), intent(in) :: self
type(unstructured_grid), intent(inout) :: dx

!integer :: ncol
!type(column_element), pointer :: current, prev

! Multiply with NICAS
call log%info("NICAS multiply")

! Count columns
!current = dx%head
!ncol = 0
!do while (associated(current))
!   ncol = ncol+1
!   current => current%next
!end do
!write(*,*) "Number of columns:",ncol
! Allocation


! Loop over 3D variables
!do ivars=1,dx%head%column%nvars
!   do while (associated(current))
      !
!      current%column%cols

      ! Next column
!      current => current%next
!   end do
!enddo

! Loop over 2D variables
!do ivars=1,dx%head%column%nsurfs


!end do

end subroutine nicas_multiply

!-------------------------------------------------------------------------------

end module nicas_mod
