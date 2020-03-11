! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module oobump_mod

use atlas_module
use datetime_mod
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use kinds
use missing_values_mod
use oops_variables_mod
use type_bump, only: bump_type
use type_nam, only: nvmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax 

implicit none

private
public :: oobump_type
public :: oobump_registry
public :: oobump_create, oobump_delete, oobump_get_cv_size, oobump_add_member, oobump_remove_member, oobump_run_drivers, &
        & oobump_multiply_vbal, oobump_multiply_vbal_inv, oobump_multiply_vbal_ad, oobump_multiply_vbal_inv_ad, &
        & oobump_multiply_nicas, oobump_multiply_nicas_sqrt, oobump_multiply_nicas_sqrt_ad, &
        & oobump_randomize_nicas, oobump_get_param, oobump_set_param
! ------------------------------------------------------------------------------
type oobump_type
   integer :: separate_log !> Separate log files for BUMP
   type(bump_type) :: bump !> Instances of BUMP
end type oobump_type

#define LISTED_TYPE oobump_type

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: oobump_registry
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
!-------------------------------------------------------------------------------
!> Create OOBUMP
subroutine oobump_create(self, f_comm, afunctionspace, afieldset, fconf, fgrid, &
 & ens1_ne, ens1_nsub, ens2_ne, ens2_nsub)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self                !< OOBUMP
type(fckit_mpi_comm),intent(in) :: f_comm               !< FCKIT MPI communicator wrapper
type(atlas_functionspace), intent(in) :: afunctionspace !< ATLAS function space
type(atlas_fieldset), intent(in) :: afieldset           !< ATLAS fieldset
type(fckit_configuration),intent(in) :: fconf           !< FCKIT configuration
type(fckit_configuration),intent(in) :: fgrid           !< FCKIT grid configuration
integer, intent(in) :: ens1_ne                          !< First ensemble size
integer, intent(in) :: ens1_nsub                        !< Number of sub-ensembles in the first ensemble
integer, intent(in) :: ens2_ne                          !< Second ensemble size
integer, intent(in) :: ens2_nsub                        !< Number of sub-ensembles in the second ensemble

! Local variables
integer :: lunit, grid_index, iproc, ifileunit
real(kind_real) :: msvalr
integer(c_size_t),parameter :: csize = 1024
character(len=1024) :: filename
character(len=:),allocatable :: str
character(kind=c_char,len=1024),allocatable :: char_array(:)

! Initialization
self%separate_log = 0
lunit = -999
if (fconf%has("separate_log")) call fconf%get_or_die("separate_log",self%separate_log)

! Initialize namelist
call self%bump%nam%init(f_comm%size())

! Read configuration
call self%bump%nam%from_conf(fconf)

! Get missing value
msvalr = missing_value(1.0_kind_real)

! Get grid index
call fgrid%get_or_die("grid_index", grid_index)

! Get number of levels
call fgrid%get_or_die("nl", self%bump%nam%nl)

! Get variables
call fgrid%get_or_die("variables",csize,char_array)
self%bump%nam%nv = size(char_array)
self%bump%nam%varname(1:self%bump%nam%nv) = char_array(1:self%bump%nam%nv)

! Get timeslots
call fgrid%get_or_die("timeslots",csize,char_array)
self%bump%nam%nts = size(char_array)
self%bump%nam%timeslot(1:self%bump%nam%nts) = char_array(1:self%bump%nam%nts)

! Get 2D level
call fgrid%get_or_die("lev2d",str)
self%bump%nam%lev2d = str

! Add suffix for multiple grid case
write(self%bump%nam%prefix,'(a,a,i2.2)') trim(self%bump%nam%prefix),'_',grid_index

! Open separate log files for BUMP
if (self%separate_log==1) then
   ! Initialize MPI
   call self%bump%mpl%init(f_comm)

   do iproc=1,self%bump%mpl%nproc
      if ((trim(self%bump%nam%verbosity)=='all').or.((trim(self%bump%nam%verbosity)=='main') &
   & .and.(iproc==self%bump%mpl%rootproc))) then
         if (iproc==self%bump%mpl%myproc) then
            ! Find a free unit
            call self%bump%mpl%newunit(lunit)

            ! Open listing file
            write(filename,'(a,i4.4,a)') trim(self%bump%nam%prefix)//'.',self%bump%mpl%myproc-1,'.out'
            inquire(file=filename,number=ifileunit)
            if (ifileunit<0) then
               open(unit=lunit,file=trim(filename),action='write',status='replace')
            else
               close(ifileunit)
               open(unit=lunit,file=trim(filename),action='write',status='replace')
            end if
         end if
         call self%bump%mpl%f_comm%barrier
      end if
   end do
end if

! Setup BUMP
call self%bump%setup(f_comm,afunctionspace,afieldset,ens1_ne=ens1_ne,ens1_nsub=ens1_nsub,ens2_ne=ens2_ne,ens2_nsub=ens2_nsub, &
 & lunit=lunit,msvalr=msvalr)

end subroutine oobump_create
!-------------------------------------------------------------------------------
!> Delete OOBUMP
subroutine oobump_delete(self)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP

! Close log files
if ((self%separate_log==1).and.(self%bump%mpl%lunit/=-999)) close(unit=self%bump%mpl%lunit)

! Release memory      
call self%bump%dealloc

end subroutine oobump_delete
!-------------------------------------------------------------------------------
!> Get control variable size
subroutine oobump_get_cv_size(self,n)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP
integer, intent(out) :: n                !< Control variable size

! Get control variable size
call self%bump%get_cv_size(n)

end subroutine oobump_get_cv_size
!-------------------------------------------------------------------------------
!> Add ensemble member
subroutine oobump_add_member(self,afieldset,ie,iens)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset
integer, intent(in) :: ie                        !< Ensemble member index
integer, intent(in) :: iens                      !< Ensemble index

! Add member
call self%bump%add_member(afieldset,ie,iens)

end subroutine oobump_add_member
!-------------------------------------------------------------------------------
!> Remove ensemble member
subroutine oobump_remove_member(self,afieldset,ie,iens)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset
integer, intent(in) :: ie                        !< Ensemble member index
integer, intent(in) :: iens                      !< Ensemble index

! Remove member
call self%bump%remove_member(afieldset,ie,iens)

end subroutine oobump_remove_member
!-------------------------------------------------------------------------------
!> Run OOBUMP drivers
subroutine oobump_run_drivers(self)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP

! Run BUMP drivers
call self%bump%run_drivers

end subroutine oobump_run_drivers
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator
subroutine oobump_multiply_vbal(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance
call self%bump%apply_vbal(afieldset)

end subroutine oobump_multiply_vbal
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator inverse
subroutine oobump_multiply_vbal_inv(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance, inverse
call self%bump%apply_vbal_inv(afieldset)

end subroutine oobump_multiply_vbal_inv
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint
subroutine oobump_multiply_vbal_ad(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance, adjoint
call self%bump%apply_vbal_ad(afieldset)

end subroutine oobump_multiply_vbal_ad
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint inverse
subroutine oobump_multiply_vbal_inv_ad(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply vertical balance, inverse adjoint
call self%bump%apply_vbal_inv_ad(afieldset)

end subroutine oobump_multiply_vbal_inv_ad
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator
subroutine oobump_multiply_nicas(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply NICAS
call self%bump%apply_nicas(afieldset)

end subroutine oobump_multiply_nicas
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root
subroutine oobump_multiply_nicas_sqrt(self,cv,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
real(kind_real), intent(in) :: cv(:)             !< Control variable
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Apply NICAS square-root
call self%bump%apply_nicas_sqrt(cv, afieldset)

end subroutine oobump_multiply_nicas_sqrt
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root adjoint
subroutine oobump_multiply_nicas_sqrt_ad(self,afieldset,cv)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset
real(kind_real), intent(inout) :: cv(:)          !< Control variable

! Apply NICAS square-root
call self%bump%apply_nicas_sqrt_ad(afieldset, cv)

end subroutine oobump_multiply_nicas_sqrt_ad
!-------------------------------------------------------------------------------
!> Randomize the BUMP NICAS operator
subroutine oobump_randomize_nicas(self,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Randomize NICAS
call self%bump%randomize(afieldset)

end subroutine oobump_randomize_nicas
!-------------------------------------------------------------------------------
!> Get BUMP parameter
subroutine oobump_get_param(self,param,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
character(len=*), intent(in) :: param            !< Parameter name
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Get parameter
call self%bump%get_parameter(param,afieldset)

end subroutine oobump_get_param
!-------------------------------------------------------------------------------
!> Set BUMP parameter
subroutine oobump_set_param(self,param,afieldset)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self         !< OOBUMP
character(len=*), intent(in) :: param            !< Parameter name
type(atlas_fieldset), intent(inout) :: afieldset !< ATLAS fieldset

! Set parameter
call self%bump%set_parameter(param,afieldset)

end subroutine oobump_set_param
!-------------------------------------------------------------------------------
end module oobump_mod
