!----------------------------------------------------------------------
! Module: type_com
!> Purpose: communication derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_com

use netcdf
use tools_missing, only: msi,msr
use tools_nc, only: ncfloat,ncerr
implicit none

! Communication derived type
type comtype
   character(len=1024) :: prefix         !< Communication prefix
   integer :: nproc                      !< Number of tasks
   integer,allocatable :: jhalocounts(:) !< Halo counts
   integer,allocatable :: jexclcounts(:) !< Exclusive interior counts
   integer,allocatable :: jhalodispl(:)  !< Halo displacement
   integer,allocatable :: jexcldispl(:)  !< Exclusive interior displacement
   integer :: nhalo                      !< Halo buffer size
   integer :: nexcl                      !< Exclusive interior buffer size
   integer,allocatable :: halo(:)        !< Halo buffer
   integer,allocatable :: excl(:)        !< Exclusive interior buffer
end type comtype

private
public :: comtype
public :: com_read,com_write

contains

!----------------------------------------------------------------------
! Subroutine: com_read
!> Purpose: read communications from a NetCDF file
!----------------------------------------------------------------------
subroutine com_read(ncid,prefix,com)

implicit none

! Passed variables
integer,intent(in) :: ncid            !< NetCDF file id
character(len=*),intent(in) :: prefix !< Communication prefix
type(comtype),intent(inout) :: com    !< Communication

! Local variables
integer :: info
integer :: nproc_id,nhalo_id,nexcl_id,jhalocounts_id,jexclcounts_id,jhalodispl_id,jexcldispl_id,halo_id,excl_id
character(len=1024) :: subr = 'com_read'

! Copy prefix
com%prefix = prefix

! Get dimensions
call ncerr(subr,nf90_inq_dimid(ncid,'nproc',nproc_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nproc_id,len=com%nproc))
info = nf90_inq_dimid(ncid,trim(prefix)//'_nhalo',nhalo_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,nhalo_id,len=com%nhalo))
else
   com%nhalo = 0
end if
info = nf90_inq_dimid(ncid,trim(prefix)//'_nexcl',nexcl_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,nexcl_id,len=com%nexcl))
else
   com%nexcl = 0
end if

! Allocation
allocate(com%jhalocounts(com%nproc))
allocate(com%jexclcounts(com%nproc))
allocate(com%jhalodispl(com%nproc))
allocate(com%jexcldispl(com%nproc))
if (com%nhalo>0) allocate(com%halo(com%nhalo))
if (com%nexcl>0) allocate(com%excl(com%nexcl))

! Get variables id
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jhalocounts',jhalocounts_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jexclcounts',jexclcounts_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jhalodispl',jhalodispl_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jexcldispl',jexcldispl_id))
if (com%nhalo>0) call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_halo',halo_id))
if (com%nexcl>0) call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_excl',excl_id))

! Get variable
call ncerr(subr,nf90_get_var(ncid,jhalocounts_id,com%jhalocounts))
call ncerr(subr,nf90_get_var(ncid,jexclcounts_id,com%jexclcounts))
call ncerr(subr,nf90_get_var(ncid,jhalodispl_id,com%jhalodispl))
call ncerr(subr,nf90_get_var(ncid,jexcldispl_id,com%jexcldispl))
if (com%nhalo>0) call ncerr(subr,nf90_get_var(ncid,halo_id,com%halo))
if (com%nexcl>0) call ncerr(subr,nf90_get_var(ncid,excl_id,com%excl))

end subroutine com_read

!----------------------------------------------------------------------
! Subroutine: com_write
!> Purpose: write communications to a NetCDF file
!----------------------------------------------------------------------
subroutine com_write(ncid,com)

implicit none

! Passed variables
integer,intent(in) :: ncid      !< NetCDF file id
type(comtype),intent(in) :: com !< Communication

! Local variables
integer :: info
integer :: nproc_id,nhalo_id,nexcl_id,jhalocounts_id,jexclcounts_id,jhalodispl_id,jexcldispl_id,halo_id,excl_id
character(len=1024) :: subr = 'com_write'

! Start definition mode
call ncerr(subr,nf90_redef(ncid))

! Define dimensions
info = nf90_inq_dimid(ncid,'nproc',nproc_id)
if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nproc',com%nproc,nproc_id))
if (com%nhalo>0) call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nhalo',com%nhalo,nhalo_id))
if (com%nexcl>0) call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nexcl',com%nexcl,nexcl_id))

! Define variables
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalocounts',nf90_int,(/nproc_id/),jhalocounts_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexclcounts',nf90_int,(/nproc_id/),jexclcounts_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalodispl',nf90_int,(/nproc_id/),jhalodispl_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexcldispl',nf90_int,(/nproc_id/),jexcldispl_id))
if (com%nhalo>0) call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_halo',nf90_int,(/nhalo_id/),halo_id))
if (com%nexcl>0) call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_excl',nf90_int,(/nexcl_id/),excl_id))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Put variables
call ncerr(subr,nf90_put_var(ncid,jhalocounts_id,com%jhalocounts))
call ncerr(subr,nf90_put_var(ncid,jexclcounts_id,com%jexclcounts))
call ncerr(subr,nf90_put_var(ncid,jhalodispl_id,com%jhalodispl))
call ncerr(subr,nf90_put_var(ncid,jexcldispl_id,com%jexcldispl))
if (com%nhalo>0) call ncerr(subr,nf90_put_var(ncid,halo_id,com%halo))
if (com%nexcl>0) call ncerr(subr,nf90_put_var(ncid,excl_id,com%excl))

end subroutine com_write

end module type_com
