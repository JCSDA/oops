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

use netcdf
use omp_lib
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr
use type_geom, only: geomtype
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_alltoallv
use tools_nc, only: ncfloat,ncerr
use yomhook, only: lhook,dr_hook

implicit none

! Communication derived type
type comtype
   ! Setup data
   integer,allocatable :: ext_to_proc(:) !< Extended index to processor
   integer,allocatable :: ext_to_red(:)  !< Extended index to reduced index

   ! Communication data
   character(len=1024) :: prefix         !< Communication prefix
   integer :: nred                       !< Reduction size
   integer :: next                       !< Extension size
   integer,allocatable :: red_to_ext(:)  !< Indices conversion
   integer :: nhalo                      !< Halo buffer size
   integer :: nexcl                      !< Exclusive interior buffer size
   integer,allocatable :: jhalocounts(:) !< Halo counts
   integer,allocatable :: jexclcounts(:) !< Exclusive interior counts
   integer,allocatable :: jhalodispl(:)  !< Halo displacement
   integer,allocatable :: jexcldispl(:)  !< Exclusive interior displacement
   integer,allocatable :: halo(:)        !< Halo buffer
   integer,allocatable :: excl(:)        !< Exclusive interior buffer
end type comtype

interface com_ext
  module procedure com_ext_1d
  module procedure com_ext_2d
end interface

interface com_red
  module procedure com_red_1d
  module procedure com_red_2d
end interface

private
public :: comtype
public :: com_dealloc,com_bcast,com_setup,com_ext,com_red,com_read,com_write

contains

!----------------------------------------------------------------------
! Subroutine: com_dealloc
!> Purpose: communications object deallocation
!----------------------------------------------------------------------
subroutine com_dealloc(com)

implicit none

! Passed variables
type(comtype),intent(inout) :: com !< Communication data

! Release memory
if (allocated(com%red_to_ext)) deallocate(com%red_to_ext)
if (allocated(com%jhalocounts)) deallocate(com%jhalocounts)
if (allocated(com%jexclcounts)) deallocate(com%jexclcounts)
if (allocated(com%jhalodispl)) deallocate(com%jhalodispl)
if (allocated(com%jexcldispl)) deallocate(com%jexcldispl)
if (allocated(com%halo)) deallocate(com%halo)
if (allocated(com%excl)) deallocate(com%excl)

end subroutine com_dealloc

!----------------------------------------------------------------------
! Subroutine: com_setup
!> Purpose: setup communications
!----------------------------------------------------------------------
subroutine com_setup(com)

implicit none

! Passed variables
type(comtype),intent(inout) :: com(mpl%nproc) !< Communication data

! Local variables
integer :: iproc,jproc,iext,icount

! Allocation
do iproc=1,mpl%nproc
   allocate(com(iproc)%jhalocounts(mpl%nproc))
   allocate(com(iproc)%jexclcounts(mpl%nproc))
   allocate(com(iproc)%jhalodispl(mpl%nproc))
   allocate(com(iproc)%jexcldispl(mpl%nproc))
end do

! Initialization
do iproc=1,mpl%nproc
   com(iproc)%jhalocounts = 0
   com(iproc)%jexclcounts = 0
   com(iproc)%jhalodispl = 0
   com(iproc)%jexcldispl = 0
end do

! Compute counts
do iproc=1,mpl%nproc
   do iext=1,com(iproc)%next
      jproc = com(iproc)%ext_to_proc(iext)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         com(iproc)%jhalocounts(jproc) = com(iproc)%jhalocounts(jproc)+1

         ! Count of points received on JPROC from IPROC
         com(jproc)%jexclcounts(iproc) = com(jproc)%jexclcounts(iproc)+1
      end if
   end do
end do

! Compute displacement
do iproc=1,mpl%nproc
   com(iproc)%jhalodispl(1) = 0
   com(iproc)%jexcldispl(1) = 0
   do jproc=2,mpl%nproc
      com(iproc)%jhalodispl(jproc) = com(iproc)%jhalodispl(jproc-1)+com(iproc)%jhalocounts(jproc-1)
      com(iproc)%jexcldispl(jproc) = com(iproc)%jexcldispl(jproc-1)+com(iproc)%jexclcounts(jproc-1)
   end do
end do

! Allocation
do iproc=1,mpl%nproc
   com(iproc)%nhalo = sum(com(iproc)%jhalocounts)
   com(iproc)%nexcl = sum(com(iproc)%jexclcounts)
   allocate(com(iproc)%halo(com(iproc)%nhalo))
   allocate(com(iproc)%excl(com(iproc)%nexcl))
end do

! Fill halo array
do iproc=1,mpl%nproc
   com(iproc)%jhalocounts = 0
end do
do iproc=1,mpl%nproc
   do iext=1,com(iproc)%next
      ! Check for halo points
      jproc = com(iproc)%ext_to_proc(iext)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         com(iproc)%jhalocounts(jproc) = com(iproc)%jhalocounts(jproc)+1

         ! Local index of points sent from IPROC to JPROC
         com(iproc)%halo(com(iproc)%jhalodispl(jproc)+com(iproc)%jhalocounts(jproc)) = iext
      end if
   end do
end do

! Fill excl array
do jproc=1,mpl%nproc
   ! Loop over processors sending data to JPROC
   do iproc=1,mpl%nproc
      do icount=1,com(iproc)%jhalocounts(jproc)
         ! Local index of points received on JPROC from IPROC
         com(jproc)%excl(com(jproc)%jexcldispl(iproc)+icount) = &
       & com(iproc)%ext_to_red(com(iproc)%halo(com(iproc)%jhalodispl(jproc)+icount))
      end do
   end do
end do

end subroutine com_setup

!----------------------------------------------------------------------
! Subroutine: com_bcast
!> Purpose: communications object broadcast
!----------------------------------------------------------------------
subroutine com_bcast(com_in,com_out)

implicit none

! Passed variables
type(comtype),intent(in) :: com_in(mpl%nproc) !< Input communication data
type(comtype),intent(inout) :: com_out        !< Output communication data

! Local variables
integer :: iproc

! Communicate dimensions
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         com_out%nred = com_in(iproc)%nred
         com_out%next = com_in(iproc)%next
         com_out%nhalo = com_in(iproc)%nhalo
         com_out%nexcl = com_in(iproc)%nexcl
      else
         ! Send dimensions to iproc
         call mpl_send(com_in(iproc)%nred,iproc,mpl%tag)
         call mpl_send(com_in(iproc)%next,iproc,mpl%tag+1)
         call mpl_send(com_in(iproc)%nhalo,iproc,mpl%tag+2)
         call mpl_send(com_in(iproc)%nexcl,iproc,mpl%tag+3)
      end if
   end do
else
   ! Receive dimensions from ioproc
   call mpl_recv(com_out%nred,mpl%ioproc,mpl%tag)
   call mpl_recv(com_out%next,mpl%ioproc,mpl%tag+1)
   call mpl_recv(com_out%nhalo,mpl%ioproc,mpl%tag+2)
   call mpl_recv(com_out%nexcl,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4

! Allocation
allocate(com_out%red_to_ext(com_out%nred))
allocate(com_out%jhalocounts(mpl%nproc))
allocate(com_out%jexclcounts(mpl%nproc))
allocate(com_out%jhalodispl(mpl%nproc))
allocate(com_out%jexcldispl(mpl%nproc))
allocate(com_out%halo(com_out%nhalo))
allocate(com_out%excl(com_out%nexcl))

! Communicate data
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         com_out%red_to_ext = com_in(iproc)%red_to_ext
         com_out%jhalocounts = com_in(iproc)%jhalocounts
         com_out%jexclcounts = com_in(iproc)%jexclcounts
         com_out%jhalodispl = com_in(iproc)%jhalodispl
         com_out%jexcldispl = com_in(iproc)%jexcldispl
         com_out%halo = com_in(iproc)%halo
         com_out%excl = com_in(iproc)%excl
      else
         ! Send dimensions to iproc
         call mpl_send(com_in(iproc)%nred,com_in(iproc)%red_to_ext,iproc,mpl%tag)
         call mpl_send(mpl%nproc,com_in(iproc)%jhalocounts,iproc,mpl%tag+1)
         call mpl_send(mpl%nproc,com_in(iproc)%jexclcounts,iproc,mpl%tag+2)
         call mpl_send(mpl%nproc,com_in(iproc)%jhalodispl,iproc,mpl%tag+3)
         call mpl_send(mpl%nproc,com_in(iproc)%jexcldispl,iproc,mpl%tag+4)
         if (com_in(iproc)%nhalo>0) call mpl_send(com_in(iproc)%nhalo,com_in(iproc)%halo,iproc,mpl%tag+5)
         if (com_in(iproc)%nexcl>0) call mpl_send(com_in(iproc)%nexcl,com_in(iproc)%excl,iproc,mpl%tag+6)
      end if
   end do
else
   ! Receive dimensions from ioproc
   call mpl_recv(com_out%nred,com_out%red_to_ext,mpl%ioproc,mpl%tag)
   call mpl_recv(mpl%nproc,com_out%jhalocounts,mpl%ioproc,mpl%tag+1)
   call mpl_recv(mpl%nproc,com_out%jexclcounts,mpl%ioproc,mpl%tag+2)
   call mpl_recv(mpl%nproc,com_out%jhalodispl,mpl%ioproc,mpl%tag+3)
   call mpl_recv(mpl%nproc,com_out%jexcldispl,mpl%ioproc,mpl%tag+4)
   if (com_out%nhalo>0) call mpl_recv(com_out%nhalo,com_out%halo,mpl%ioproc,mpl%tag+5)
   if (com_out%nexcl>0) call mpl_recv(com_out%nexcl,com_out%excl,mpl%ioproc,mpl%tag+6)
end if
mpl%tag = mpl%tag+7

end subroutine com_bcast

!----------------------------------------------------------------------
! Subroutine: com_ext_1d
!> Purpose: communicate field to halo (extension), 1d
!----------------------------------------------------------------------
subroutine com_ext_1d(com,vec)

implicit none

! Passed variables
type(comtype),intent(in) :: com                     !< Communication data
real(kind_real),allocatable,intent(inout) :: vec(:) !< Vector

! Local variables
integer :: iexcl,ired,ihalo
real(kind_real) :: sbuf(com%nexcl),rbuf(com%nhalo),vec_tmp(com%nred)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('com_ext_1d',0,zhook_handle)

! Check input vector size
if (size(vec)/=com%nred) call msgerror('vector size inconsistent in com_ext')

! Prepare buffers to send
!$omp parallel do schedule(static) private(iexcl)
do iexcl=1,com%nexcl
   sbuf(iexcl) = vec(com%excl(iexcl))
end do
!$omp end parallel do

! Communication
call mpl_alltoallv(com%nexcl,sbuf,com%jexclcounts,com%jexcldispl,com%nhalo,rbuf,com%jhalocounts,com%jhalodispl)

! Copy
vec_tmp = vec

! Reallocation
deallocate(vec)
allocate(vec(com%next))

! Copy interior
!$omp parallel do schedule(static) private(ired)
do ired=1,com%nred
   vec(com%red_to_ext(ired)) = vec_tmp(ired)
end do
!$omp end parallel do

! Copy halo
!$omp parallel do schedule(static) private(ihalo)
do ihalo=1,com%nhalo
   vec(com%halo(ihalo)) = rbuf(ihalo)
end do
!$omp end parallel do

if (lhook) call dr_hook('com_ext_1d',1,zhook_handle)

end subroutine com_ext_1d

!----------------------------------------------------------------------
! Subroutine: com_ext_2d
!> Purpose: communicate field to halo (extension), 2d
!----------------------------------------------------------------------
subroutine com_ext_2d(com,vec)

implicit none

! Passed variables
type(comtype),intent(in) :: com                       !< Communication data
real(kind_real),allocatable,intent(inout) :: vec(:,:) !< Vector

! Local variables
integer :: nl,il,iexcl,ired,ihalo
real(kind_real),allocatable :: sbuf(:),rbuf(:),vec_tmp(:,:)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('com_ext_2d',0,zhook_handle)

! Check input vector size
if (size(vec,1)/=com%nred) call msgerror('vector size inconsistent in com_ext')

! Second dimension
nl = size(vec,2)

! Allocation
allocate(sbuf(com%nexcl*nl))
allocate(rbuf(com%nhalo*nl))
allocate(vec_tmp(com%nred,nl))

! Prepare buffers to send
!$omp parallel do schedule(static) private(il,iexcl)
do il=1,nl
   do iexcl=1,com%nexcl
      sbuf((iexcl-1)*nl+il) = vec(com%excl(iexcl),il)
   end do
end do
!$omp end parallel do

! Communication
call mpl_alltoallv(com%nexcl*nl,sbuf,com%jexclcounts*nl,com%jexcldispl*nl,com%nhalo*nl,rbuf,com%jhalocounts*nl,com%jhalodispl*nl)

! Copy
vec_tmp = vec

! Reallocation
deallocate(vec)
allocate(vec(com%next,nl))

! Copy interior
!$omp parallel do schedule(static) private(il,ired)
do il=1,nl
   do ired=1,com%nred
      vec(com%red_to_ext(ired),il) = vec_tmp(ired,il)
   end do
end do
!$omp end parallel do

! Copy halo
!$omp parallel do schedule(static) private(il,ihalo)
do il=1,nl
   do ihalo=1,com%nhalo
      vec(com%halo(ihalo),il) = rbuf((ihalo-1)*nl+il)
   end do
end do
!$omp end parallel do

if (lhook) call dr_hook('com_ext_2d',1,zhook_handle)

end subroutine com_ext_2d

!----------------------------------------------------------------------
! Subroutine: com_red_1d
!> Purpose: communicate vector from halo (reduction)
!----------------------------------------------------------------------
subroutine com_red_1d(com,vec)

implicit none

! Passed variables
type(comtype),intent(in) :: com                     !< Communication data
real(kind_real),allocatable,intent(inout) :: vec(:) !< Vector

! Local variables
integer :: ihalo,ired,iexcl,ithread
real(kind_real) :: sbuf(com%nhalo),rbuf(com%nexcl),vec_tmp(com%next),vec_arr(com%nred,mpl%nthread)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('com_red_1d',0,zhook_handle)

! Check input vector size
if (size(vec)/=com%next) call msgerror('vector size inconsistent in com_red_1d')

! Prepare buffers to send
!$omp parallel do schedule(static) private(ihalo)
do ihalo=1,com%nhalo
   sbuf(ihalo) = vec(com%halo(ihalo))
end do
!$omp end parallel do

! Communication
call mpl_alltoallv(com%nhalo,sbuf,com%jhalocounts,com%jhalodispl,com%nexcl,rbuf,com%jexclcounts,com%jexcldispl)

! Copy
vec_tmp = vec

! Reallocation
deallocate(vec)
allocate(vec(com%nred))

! Copy interior
!$omp parallel do schedule(static) private(ired)
do ired=1,com%nred
   vec(ired) = vec_tmp(com%red_to_ext(ired))
end do
!$omp end parallel do

! Copy halo
vec_arr = 0.0
!$omp parallel do schedule(static) private(iexcl,ithread)
do iexcl=1,com%nexcl
   ithread = omp_get_thread_num()+1
   vec_arr(com%excl(iexcl),ithread) = vec_arr(com%excl(iexcl),ithread)+rbuf(iexcl)
end do
!$omp end parallel do

! Sum over threads
do ithread=1,mpl%nthread
   vec = vec+vec_arr(:,ithread)
end do

if (lhook) call dr_hook('com_red_1d',1,zhook_handle)

end subroutine com_red_1d

!----------------------------------------------------------------------
! Subroutine: com_red_2d
!> Purpose: communicate vector from halo (reduction)
!----------------------------------------------------------------------
subroutine com_red_2d(com,vec)

implicit none

! Passed variables
type(comtype),intent(in) :: com                       !< Communication data
real(kind_real),allocatable,intent(inout) :: vec(:,:) !< Vector

! Local variables
integer :: nl,il,ihalo,ired,iexcl,ithread
real(kind_real),allocatable :: sbuf(:),rbuf(:),vec_tmp(:,:),vec_arr(:,:,:)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('com_red_2d',0,zhook_handle)

! Check input vector size
if (size(vec,1)/=com%next) call msgerror('vector size inconsistent in com_red_2d')

! Second dimension
nl = size(vec,2)

! Allocation
allocate(sbuf(com%nhalo*nl))
allocate(rbuf(com%nexcl*nl))
allocate(vec_tmp(com%next,nl))
allocate(vec_arr(com%nred,nl,mpl%nthread))

! Prepare buffers to send
!$omp parallel do schedule(static) private(il,ihalo)
do il=1,nl
   do ihalo=1,com%nhalo
      sbuf((ihalo-1)*nl+il) = vec(com%halo(ihalo),il)
   end do
end do
!$omp end parallel do

! Communication
call mpl_alltoallv(com%nhalo*nl,sbuf,com%jhalocounts*nl,com%jhalodispl*nl,com%nexcl*nl,rbuf,com%jexclcounts*nl,com%jexcldispl*nl)

! Copy
vec_tmp = vec

! Reallocation
deallocate(vec)
allocate(vec(com%nred,nl))

! Copy interior
!$omp parallel do schedule(static) private(il,ired)
do il=1,nl
   do ired=1,com%nred
      vec(ired,il) = vec_tmp(com%red_to_ext(ired),il)
   end do
end do
!$omp end parallel do

! Copy halo
vec_arr = 0.0
!$omp parallel do schedule(static) private(il,iexcl,ithread)
do il=1,nl
   do iexcl=1,com%nexcl
      ithread = omp_get_thread_num()+1
      vec_arr(com%excl(iexcl),il,ithread) = vec_arr(com%excl(iexcl),il,ithread)+rbuf((iexcl-1)*nl+il)
   end do
end do
!$omp end parallel do

! Sum over threads
do ithread=1,mpl%nthread
   vec = vec+vec_arr(:,:,ithread)
end do

if (lhook) call dr_hook('com_red_2d',1,zhook_handle)

end subroutine com_red_2d

!----------------------------------------------------------------------
! Subroutine: com_read
!> Purpose: read communications from a NetCDF file
!----------------------------------------------------------------------
subroutine com_read(ncid,prefix,com)

implicit none

! Passed variables
integer,intent(in) :: ncid            !< NetCDF file id
character(len=*),intent(in) :: prefix !< Communication prefix
type(comtype),intent(inout) :: com    !< Communication data

! Local variables
integer :: info
integer :: nred_id,next_id,red_to_ext_id,nhalo_id,nexcl_id
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
allocate(com%red_to_ext(com%nred))
allocate(com%jhalocounts(mpl%nproc))
allocate(com%jexclcounts(mpl%nproc))
allocate(com%jhalodispl(mpl%nproc))
allocate(com%jexcldispl(mpl%nproc))
if (com%nhalo>0) allocate(com%halo(com%nhalo))
if (com%nexcl>0) allocate(com%excl(com%nexcl))

! Get variables id
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_red_to_ext',red_to_ext_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jhalocounts',jhalocounts_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jexclcounts',jexclcounts_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jhalodispl',jhalodispl_id))
call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_jexcldispl',jexcldispl_id))
if (com%nhalo>0) call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_halo',halo_id))
if (com%nexcl>0) call ncerr(subr,nf90_inq_varid(ncid,trim(prefix)//'_excl',excl_id))

! Get variable
call ncerr(subr,nf90_get_var(ncid,red_to_ext_id,com%red_to_ext))
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
type(comtype),intent(in) :: com !< Communication data

! Local variables
integer :: info
integer :: nproc_id,nred_id,next_id,red_to_ext_id,nhalo_id,nexcl_id
integer :: jhalocounts_id,jexclcounts_id,jhalodispl_id,jexcldispl_id,halo_id,excl_id
character(len=1024) :: subr = 'com_write'

! Start definition mode
call ncerr(subr,nf90_redef(ncid))

! Define dimensions
info = nf90_inq_dimid(ncid,'nproc',nproc_id)
if (info/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nproc',mpl%nproc,nproc_id))
call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nred',com%nred,nred_id))
call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_next',com%next,next_id))
if (com%nhalo>0) call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nhalo',com%nhalo,nhalo_id))
if (com%nexcl>0) call ncerr(subr,nf90_def_dim(ncid,trim(com%prefix)//'_nexcl',com%nexcl,nexcl_id))

! Define variables
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_red_to_ext',nf90_int,(/nred_id/),red_to_ext_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalocounts',nf90_int,(/nproc_id/),jhalocounts_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexclcounts',nf90_int,(/nproc_id/),jexclcounts_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jhalodispl',nf90_int,(/nproc_id/),jhalodispl_id))
call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_jexcldispl',nf90_int,(/nproc_id/),jexcldispl_id))
if (com%nhalo>0) call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_halo',nf90_int,(/nhalo_id/),halo_id))
if (com%nexcl>0) call ncerr(subr,nf90_def_var(ncid,trim(com%prefix)//'_excl',nf90_int,(/nexcl_id/),excl_id))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Put variables
call ncerr(subr,nf90_put_var(ncid,red_to_ext_id,com%red_to_ext))
call ncerr(subr,nf90_put_var(ncid,jhalocounts_id,com%jhalocounts))
call ncerr(subr,nf90_put_var(ncid,jexclcounts_id,com%jexclcounts))
call ncerr(subr,nf90_put_var(ncid,jhalodispl_id,com%jhalodispl))
call ncerr(subr,nf90_put_var(ncid,jexcldispl_id,com%jexcldispl))
if (com%nhalo>0) call ncerr(subr,nf90_put_var(ncid,halo_id,com%halo))
if (com%nexcl>0) call ncerr(subr,nf90_put_var(ncid,excl_id,com%excl))

end subroutine com_write

end module type_com
