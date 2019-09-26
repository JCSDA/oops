! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module oobump_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use missing_values_mod
use type_bump, only: bump_type
use type_nam, only: nvmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax 
use unstructured_grid_mod

implicit none

private
public :: oobump_type
public :: oobump_registry
public :: create_oobump, delete_oobump, get_oobump_colocated, get_oobump_nts, get_oobump_cv_size, add_oobump_member, &
        & remove_oobump_member,run_oobump_drivers,multiply_oobump_vbal, multiply_oobump_vbal_inv, multiply_oobump_vbal_ad,  &
        & multiply_oobump_vbal_inv_ad,multiply_oobump_nicas, multiply_oobump_nicas_sqrt, multiply_oobump_nicas_sqrt_ad, &
        & randomize_oobump_nicas, get_oobump_param, set_oobump_param, bump_read_conf
! ------------------------------------------------------------------------------
type oobump_type
   integer :: colocated                   !> Colocated flag
   integer :: separate_log                !> Separate log files for BUMP
   integer :: ngrid                       !> Number of instances of BUMP
   type(bump_type),allocatable :: bump(:) !> Instances of BUMP
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
subroutine create_oobump(self, ug, f_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub, mpi_comm)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self       !< OOBUMP
type(unstructured_grid), intent(in) :: ug      !< Unstructured grid
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
integer, intent(in) :: ens1_ne                 !< First ensemble size
integer, intent(in) :: ens1_nsub               !< Number of sub-ensembles in the first ensemble
integer, intent(in) :: ens2_ne                 !< Second ensemble size
integer, intent(in) :: ens2_nsub               !< Number of sub-ensembles in the second ensemble
integer,intent(in) :: mpi_comm                 ! MPI communicator

! Local variables
integer :: igrid, lunit, iproc, ifileunit
real(kind_real) :: msvalr
character(len=1024) :: filename

! Initialization
self%colocated = ug%colocated
self%separate_log = 0
lunit = -999
if (f_conf%has("separate_log")) call f_conf%get_or_die("separate_log",self%separate_log)
self%ngrid = ug%ngrid

! Allocation
allocate(self%bump(self%ngrid))

do igrid=1,self%ngrid
   ! Initialize namelist
   call self%bump(igrid)%nam%init

   ! Read JSON
   call bump_read_conf(f_conf,self%bump(igrid))

   ! Add suffix for multiple grid case
   if (self%ngrid>1) write(self%bump(igrid)%nam%prefix,'(a,a,i2.2)') trim(self%bump(igrid)%nam%prefix),'_',igrid

   ! Get missing value
   msvalr = missing_value(1.0_kind_real)

   ! Open separate log files for BUMP
   if (self%separate_log==1) then
      ! Initialize MPI
      call self%bump(igrid)%mpl%init(mpi_comm)

      do iproc=1,self%bump(igrid)%mpl%nproc
         if ((trim(self%bump(igrid)%nam%verbosity)=='all').or.((trim(self%bump(igrid)%nam%verbosity)=='main') &
      & .and.(iproc==self%bump(igrid)%mpl%rootproc))) then
            if (iproc==self%bump(igrid)%mpl%myproc) then
               ! Find a free unit
               call self%bump(igrid)%mpl%newunit(lunit)

               ! Open listing file
               write(filename,'(a,i4.4,a)') trim(self%bump(igrid)%nam%prefix)//'.',self%bump(igrid)%mpl%myproc-1,'.out'
               inquire(file=filename,number=ifileunit)
               if (ifileunit<0) then
                  open(unit=lunit,file=trim(filename),action='write',status='replace')
               else
                  close(ifileunit)
                  open(unit=lunit,file=trim(filename),action='write',status='replace')
               end if
            end if
            call self%bump(igrid)%mpl%f_comm%barrier
         end if
      end do
   end if

   ! Online setup
   call self%bump(igrid)%setup_online(ug%grid(igrid)%nmga,ug%grid(igrid)%nl0,ug%grid(igrid)%nv,ug%grid(igrid)%nts, &
 & ug%grid(igrid)%lon,ug%grid(igrid)%lat,ug%grid(igrid)%area,ug%grid(igrid)%vunit,ug%grid(igrid)%lmask,ens1_ne=ens1_ne, &
 & ens1_nsub=ens1_nsub,ens2_ne=ens2_ne,ens2_nsub=ens2_nsub,lunit=lunit,msvalr=msvalr)
end do

end subroutine create_oobump
!-------------------------------------------------------------------------------
!> Delete OOBUMP
subroutine delete_oobump(self)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP

! Local variables
integer :: igrid

if (allocated(self%bump)) then
   do igrid=1,self%ngrid
      ! Close log files
      if ((self%separate_log==1).and.(self%bump(igrid)%mpl%lunit/=-999)) close(unit=self%bump(igrid)%mpl%lunit)

      ! Release memory      
      call self%bump(igrid)%dealloc
   end do

   ! Release memory
   deallocate(self%bump)
end if

end subroutine delete_oobump
!-------------------------------------------------------------------------------
!> Get colocated flag
subroutine get_oobump_colocated(self,colocated)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP
integer, intent(out) :: colocated        !< Colocated flag

! Copy colocated
colocated = self%colocated

end subroutine get_oobump_colocated
!-------------------------------------------------------------------------------
!> Get number of timeslots
subroutine get_oobump_nts(self,nts)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP
integer, intent(out) :: nts              !< Number of timeslots

! Copy nts
nts = self%bump(1)%nam%nts

end subroutine get_oobump_nts
!-------------------------------------------------------------------------------
!> Get control variable size
subroutine get_oobump_cv_size(self,n)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP
integer, intent(out) :: n                !< Control variable size

! Local variables
integer :: igrid,nn

! Add control variable sizes for each grid
n = 0
do igrid=1,self%ngrid
   call self%bump(igrid)%get_cv_size(nn)
   n = n+nn
end do

end subroutine get_oobump_cv_size
!-------------------------------------------------------------------------------
!> Add ensemble member
subroutine add_oobump_member(self,ug,ie,iens)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid
integer, intent(in) :: ie                    !< Ensemble member index
integer, intent(in) :: iens                  !< Ensemble index

! Local variables
integer :: igrid

! Add member
do igrid=1,self%ngrid
   call self%bump(igrid)%add_member(ug%grid(igrid)%fld,ie,iens)
end do

end subroutine add_oobump_member
!-------------------------------------------------------------------------------
!> Remove ensemble member
subroutine remove_oobump_member(self,ug,ie,iens)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid
integer, intent(in) :: ie                    !< Ensemble member index
integer, intent(in) :: iens                  !< Ensemble index

! Local variables
integer :: igrid

! Remove member
do igrid=1,self%ngrid
   call self%bump(igrid)%remove_member(ug%grid(igrid)%fld,ie,iens)
end do

end subroutine remove_oobump_member
!-------------------------------------------------------------------------------
!> Run OOBUMP drivers
subroutine run_oobump_drivers(self)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self !< OOBUMP

! Local variables
integer :: igrid

! Run BUMP drivers
do igrid=1,self%ngrid
   call self%bump(igrid)%run_drivers
end do

end subroutine run_oobump_drivers
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator
subroutine multiply_oobump_vbal(self,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Apply vertical balance
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator inverse
subroutine multiply_oobump_vbal_inv(self,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Apply vertical balance, inverse
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal_inv(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal_inv
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint
subroutine multiply_oobump_vbal_ad(self,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Apply vertical balance, adjoint
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal_ad(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal_ad
!-------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint inverse
subroutine multiply_oobump_vbal_inv_ad(self,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Apply vertical balance, inverse adjoint
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal_inv_ad(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal_inv_ad
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator
subroutine multiply_oobump_nicas(self,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Apply NICAS
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_nicas(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_nicas
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root
subroutine multiply_oobump_nicas_sqrt(self,cv,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
real(kind_real), intent(in) :: cv(:)         !< Control variable
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: offset,igrid,nn

! Allocation
call allocate_unstructured_grid_field(ug)

! Initialization
offset = 0

do igrid=1,self%ngrid
   ! Get control variable size for this grid
   call self%bump(igrid)%get_cv_size(nn)

   ! Apply NICAS square-root
   call self%bump(igrid)%apply_nicas_sqrt(cv(offset+1:offset+nn), ug%grid(igrid)%fld)

   ! Update
   offset = offset+nn
end do

end subroutine multiply_oobump_nicas_sqrt
!-------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root adjoint
subroutine multiply_oobump_nicas_sqrt_ad(self,ug,cv)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self  !< OOBUMP
type(unstructured_grid), intent(in) :: ug !< Unstructured grid
real(kind_real), intent(inout) :: cv(:)   !< Control variable

! Local variables
integer :: offset,igrid,nn

! Initialization
offset = 0

do igrid=1,self%ngrid
   ! Get control variable size for this grid
   call self%bump(igrid)%get_cv_size(nn)

   ! Apply NICAS square-root
   call self%bump(igrid)%apply_nicas_sqrt_ad(ug%grid(igrid)%fld, cv(offset+1:offset+nn))

   ! Update
   offset = offset+nn
end do

end subroutine multiply_oobump_nicas_sqrt_ad
!-------------------------------------------------------------------------------
!> Randomize the BUMP NICAS operator
subroutine randomize_oobump_nicas(self,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Allocation
call allocate_unstructured_grid_field(ug)

! Randomize NICAS
do igrid=1,self%ngrid
   call self%bump(igrid)%randomize(ug%grid(igrid)%fld)
end do

end subroutine randomize_oobump_nicas
!-------------------------------------------------------------------------------
!> Get BUMP parameter
subroutine get_oobump_param(self,param,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self     !< OOBUMP
character(len=*), intent(in) :: param        !< Parameter name
type(unstructured_grid), intent(inout) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Allocation
call allocate_unstructured_grid_field(ug)

! Get parameter
do igrid=1,self%ngrid
   call self%bump(igrid)%get_parameter(param,ug%grid(igrid)%fld)
end do

end subroutine get_oobump_param
!-------------------------------------------------------------------------------
!> Set BUMP parameter
subroutine set_oobump_param(self,param,ug)

implicit none

! Passed variables
type(oobump_type), intent(inout) :: self  !< OOBUMP
character(len=*),intent(in) :: param      !< Parameter name
type(unstructured_grid), intent(in) :: ug !< Unstructured grid

! Local variables
integer :: igrid

! Set parameter
do igrid=1,self%ngrid
   call self%bump(igrid)%set_parameter(param,ug%grid(igrid)%fld)
end do

end subroutine set_oobump_param
!-------------------------------------------------------------------------------
!> Read BUMP configuration
subroutine bump_read_conf(f_conf,bump)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(bump_type), intent(inout) :: bump         !< BUMP

! Set BUMP namelist
call bump%nam%from_conf(f_conf)

end subroutine bump_read_conf
!-------------------------------------------------------------------------------
end module oobump_mod
