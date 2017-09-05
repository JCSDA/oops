!----------------------------------------------------------------------
! Module: type_com
!> Purpose: communications derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_com

use module_namelist, only: namtype
use netcdf
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr
use type_geom, only: geomtype
use type_mpl, only: mpl_alltoallv
use tools_nc, only: ncfloat,ncerr

implicit none

! Communication derived type
type comtype
   ! Setup data
   integer,allocatable :: iext_to_iproc(:)
   integer,allocatable :: iext_to_ired(:)

   ! Communication data
   character(len=1024) :: prefix         !< Communication prefix
   integer :: nred
   integer :: next
   integer,allocatable :: ired_to_iext(:) !< Indices conversion
   integer :: nhalo                      !< Halo buffer size
   integer :: nexcl                      !< Exclusive interior buffer size
   integer,allocatable :: jhalocounts(:) !< Halo counts
   integer,allocatable :: jexclcounts(:) !< Exclusive interior counts
   integer,allocatable :: jhalodispl(:)  !< Halo displacement
   integer,allocatable :: jexcldispl(:)  !< Exclusive interior displacement
   integer,allocatable :: halo(:)        !< Halo buffer
   integer,allocatable :: excl(:)        !< Exclusive interior buffer
end type comtype

private
public :: comtype
public :: com_dealloc,com_copy,com_setup,com_ext,com_red,com_read,com_write

contains

!----------------------------------------------------------------------
! Subroutine: com_dealloc
!> Purpose: communications object deallocation
!----------------------------------------------------------------------
subroutine com_dealloc(com)

implicit none

! Passed variables
type(comtype),intent(inout) :: com !< Linear operator

! Release memory
if (allocated(com%ired_to_iext)) deallocate(com%ired_to_iext)
if (allocated(com%jhalocounts)) deallocate(com%jhalocounts)
if (allocated(com%jexclcounts)) deallocate(com%jexclcounts)
if (allocated(com%jhalodispl)) deallocate(com%jhalodispl)
if (allocated(com%jexcldispl)) deallocate(com%jexcldispl)
if (allocated(com%halo)) deallocate(com%halo)
if (allocated(com%excl)) deallocate(com%excl)

end subroutine com_dealloc

!----------------------------------------------------------------------
! Subroutine: com_copy
!> Purpose: communications object copyment
!----------------------------------------------------------------------
subroutine com_copy(nam,com_in,com_out)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(comtype),intent(in) :: com_in     !< Input linear operator
type(comtype),intent(inout) :: com_out !< Output linear operator

! Copy attributes
com_out%prefix = trim(com_in%prefix)
com_out%nred = com_in%nred
com_out%next = com_in%next
com_out%nhalo = com_in%nhalo
com_out%nexcl = com_in%nexcl

! Deallocation
call com_dealloc(com_out)

! Allocation
allocate(com_out%ired_to_iext(nam%nproc))
allocate(com_out%jhalocounts(nam%nproc))
allocate(com_out%jexclcounts(nam%nproc))
allocate(com_out%jhalodispl(nam%nproc))
allocate(com_out%jexcldispl(nam%nproc))
allocate(com_out%halo(com_out%nhalo))
allocate(com_out%excl(com_out%nexcl))

! Copy data
com_out%ired_to_iext = com_in%ired_to_iext
com_out%jhalocounts = com_in%jhalocounts
com_out%jexclcounts = com_in%jexclcounts
com_out%jhalodispl = com_in%jhalodispl
com_out%jexcldispl = com_in%jexcldispl
com_out%halo = com_in%halo
com_out%excl = com_in%excl

end subroutine com_copy

!----------------------------------------------------------------------
! Subroutine: com_setup
!> Purpose: setup communications
!----------------------------------------------------------------------
subroutine com_setup(nam,com)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(comtype),intent(inout) :: com(nam%nproc) !< Communication

! Local variables
integer :: iproc,jproc,i

! Allocation
do iproc=1,nam%nproc
   allocate(com(iproc)%jhalocounts(nam%nproc))
   allocate(com(iproc)%jexclcounts(nam%nproc))
   allocate(com(iproc)%jhalodispl(nam%nproc))
   allocate(com(iproc)%jexcldispl(nam%nproc))
end do

! Initialization
do iproc=1,nam%nproc
   com(iproc)%jhalocounts = 0
   com(iproc)%jexclcounts = 0
   com(iproc)%jhalodispl = 0
   com(iproc)%jexcldispl = 0
end do

! Compute counts
do iproc=1,nam%nproc
   do i=1,com(iproc)%next
      jproc = com(iproc)%iext_to_iproc(i)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         com(iproc)%jhalocounts(jproc) = com(iproc)%jhalocounts(jproc)+1

         ! Count of points received on JPROC from IPROC
         com(jproc)%jexclcounts(iproc) = com(jproc)%jexclcounts(iproc)+1
      end if
   end do
end do

! Compute displacement
do iproc=1,nam%nproc
   com(iproc)%jhalodispl(1) = 0
   com(iproc)%jexcldispl(1) = 0
   do jproc=2,nam%nproc
      com(iproc)%jhalodispl(jproc) = com(iproc)%jhalodispl(jproc-1)+com(iproc)%jhalocounts(jproc-1)
      com(iproc)%jexcldispl(jproc) = com(iproc)%jexcldispl(jproc-1)+com(iproc)%jexclcounts(jproc-1)
   end do
end do

! Allocation
do iproc=1,nam%nproc
   com(iproc)%nhalo = sum(com(iproc)%jhalocounts)
   com(iproc)%nexcl = sum(com(iproc)%jexclcounts)
   allocate(com(iproc)%halo(com(iproc)%nhalo))
   allocate(com(iproc)%excl(com(iproc)%nexcl))
end do

! Fill halo array
do iproc=1,nam%nproc
   com(iproc)%jhalocounts = 0
end do
do iproc=1,nam%nproc
   do i=1,com(iproc)%next
      ! Check for halo points
      jproc = com(iproc)%iext_to_iproc(i)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         com(iproc)%jhalocounts(jproc) = com(iproc)%jhalocounts(jproc)+1

         ! Local index of points sent from IPROC to JPROC
         com(iproc)%halo(com(iproc)%jhalodispl(jproc)+com(iproc)%jhalocounts(jproc)) = i
      end if
   end do
end do

! Fill excl array
do jproc=1,nam%nproc
   ! Loop over processors sending data to JPROC
   do iproc=1,nam%nproc
      do i=1,com(iproc)%jhalocounts(jproc)
         ! Local index of points received on JPROC from IPROC
         com(jproc)%excl(com(jproc)%jexcldispl(iproc)+i) = &
       & com(iproc)%iext_to_ired(com(iproc)%halo(com(iproc)%jhalodispl(jproc)+i))
      end do
   end do
end do

end subroutine com_setup

!----------------------------------------------------------------------
! Subroutine: com_ext
!> Purpose: communicate field to halo (extension)
!----------------------------------------------------------------------
subroutine com_ext(com,vec)

implicit none

! Passed variables
type(comtype),intent(in) :: com !< Sampling data
real(kind_real),allocatable,intent(inout) :: vec(:)    !< Subgrid variable

! Local variables
real(kind_real) :: sbuf(com%nexcl),rbuf(com%nhalo),vec_tmp(com%nred)

! Check input vector size
if (size(vec)/=com%nred) call msgerror('vector size inconsistent in com_ext')

! Prepare buffers to send
sbuf = vec(com%excl)

! Communication
call mpl_alltoallv(com%nexcl,sbuf,com%jexclcounts,com%jexcldispl,com%nhalo,rbuf,com%jhalocounts,com%jhalodispl)

! Copy
vec_tmp = vec

! Reallocation
deallocate(vec)
allocate(vec(com%next))

! Copy zone A into zone B
vec(com%ired_to_iext) = vec_tmp

! Copy halo into zone B
vec(com%halo) = rbuf

end subroutine com_ext

!----------------------------------------------------------------------
! Subroutine: com_red
!> Purpose: communicate vector from halo (reduction)
!----------------------------------------------------------------------
subroutine com_red(com,vec)

implicit none

! Passed variables
type(comtype),intent(in) :: com !< Sampling data
real(kind_real),allocatable,intent(inout) :: vec(:)    !< Subgrid variable

! Local variables
real(kind_real) :: sbuf(com%nhalo),rbuf(com%nexcl),vec_tmp(com%next)

! Check input vector size
if (size(vec)/=com%next) call msgerror('vector size inconsistent in com_red')

! Prepare buffers to send
sbuf = vec(com%halo)

! Communication
call mpl_alltoallv(com%nhalo,sbuf,com%jhalocounts,com%jhalodispl,com%nexcl,rbuf,com%jexclcounts,com%jexcldispl)

! Copy
vec_tmp = vec

! Reallocation
deallocate(vec)
allocate(vec(com%nred))

! Copy zone B into zone A
vec = vec_tmp(com%ired_to_iext)

! Copy halo into zone B
vec(com%excl) = vec(com%excl)+rbuf

end subroutine com_red

!----------------------------------------------------------------------
! Subroutine: com_read
!> Purpose: read communications from a NetCDF file
!----------------------------------------------------------------------
subroutine com_read(nam,ncid,prefix,com)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
integer,intent(in) :: ncid            !< NetCDF file id
character(len=*),intent(in) :: prefix !< Communication prefix
type(comtype),intent(inout) :: com    !< Communication

! Local variables
integer :: info
integer :: nred_id,next_id,ired_to_iext_id,nhalo_id,nexcl_id
integer :: jhalocounts_id,jexclcounts_id,jhalodispl_id,jexcldispl_id,halo_id,excl_id
character(len=1024) :: subr = 'com_read'

! Copy prefix
com%prefix = trim(prefix)

! Get dimensions
info = nf90_inq_dimid(ncid,trim(prefix)//'_nred',nred_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,nred_id,len=com%nred))
else
   com%nred = 0
end if
info = nf90_inq_dimid(ncid,trim(prefix)//'_next',next_id)
if (info==nf90_noerr) then
   call ncerr(subr,nf90_inquire_dimension(ncid,next_id,len=com%next))
else
   com%next = 0
end if
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
allocate(com%ired_to_iext(com%nred))
allocate(com%jhalocounts(nam%nproc))
allocate(com%jexclcounts(nam%nproc))
allocate(com%jhalodispl(nam%nproc))
allocate(com%jexcldispl(nam%nproc))
if (com%nhalo>0) allocate(com%halo(com%nhalo))
if (com%nexcl>0) allocate(com%excl(com%nexcl))

! Get variables id
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_ired_to_iext',ired_to_iext_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jhalocounts',jhalocounts_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jexclcounts',jexclcounts_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jhalodispl',jhalodispl_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jexcldispl',jexcldispl_id))
if (com%nhalo>0) call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_halo',halo_id))
if (com%nexcl>0) call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_excl',excl_id))

! Get variable
call ncerr(subr,nf90_get_var(ncid,ired_to_iext_id,com%ired_to_iext))
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
subroutine com_write(nam,ncid,com)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
integer,intent(in) :: ncid      !< NetCDF file id
type(comtype),intent(in) :: com !< Communication

! Local variables
integer :: info
integer :: nproc_id,nred_id,next_id,ired_to_iext_id,nhalo_id,nexcl_id
integer :: jhalocounts_id,jexclcounts_id,jhalodispl_id,jexcldispl_id,halo_id,excl_id
character(len=1024) :: subr = 'com_write'

! Start definition mode
call ncerr(subr,nf90_redef(ncid))

! Define dimensions
info = nf90_inq_dimid(ncid,'nproc',nproc_id)
if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nproc',nam%nproc,nproc_id))
call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nred',com%nred,nred_id))
call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_next',com%next,next_id))
if (com%nhalo>0) call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nhalo',com%nhalo,nhalo_id))
if (com%nexcl>0) call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nexcl',com%nexcl,nexcl_id))

! Define variables
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_ired_to_iext',nf90_int,(/nred_id/),ired_to_iext_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalocounts',nf90_int,(/nproc_id/),jhalocounts_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexclcounts',nf90_int,(/nproc_id/),jexclcounts_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalodispl',nf90_int,(/nproc_id/),jhalodispl_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexcldispl',nf90_int,(/nproc_id/),jexcldispl_id))
if (com%nhalo>0) call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_halo',nf90_int,(/nhalo_id/),halo_id))
if (com%nexcl>0) call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_excl',nf90_int,(/nexcl_id/),excl_id))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Put variables
call ncerr(subr,nf90_put_var(ncid,ired_to_iext_id,com%ired_to_iext))
call ncerr(subr,nf90_put_var(ncid,jhalocounts_id,com%jhalocounts))
call ncerr(subr,nf90_put_var(ncid,jexclcounts_id,com%jexclcounts))
call ncerr(subr,nf90_put_var(ncid,jhalodispl_id,com%jhalodispl))
call ncerr(subr,nf90_put_var(ncid,jexcldispl_id,com%jexcldispl))
if (com%nhalo>0) call ncerr(subr,nf90_put_var(ncid,halo_id,com%halo))
if (com%nexcl>0) call ncerr(subr,nf90_put_var(ncid,excl_id,com%excl))

end subroutine com_write

end module type_com
