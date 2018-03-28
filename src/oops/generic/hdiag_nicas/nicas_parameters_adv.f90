!----------------------------------------------------------------------
! Module: nicas_parameters_adv.f90
!> Purpose: compute NICAS parameters (advection)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_parameters_adv

use tools_display, only: msgerror
use tools_interp, only: compute_interp
use tools_kinds,only: kind_real
use tools_missing, only: msi
use type_bdata, only: bdatatype
use type_com, only: comtype,com_dealloc,com_setup,com_bcast
use type_linop, only: linoptype,linop_alloc,linop_reorder
use type_mpl, only: mpl,mpl_send,mpl_recv
use type_nam, only: namtype
use type_ndata, only: ndatatype

implicit none

private
public :: compute_adv

contains

!----------------------------------------------------------------------
! Subroutine: compute_adv
!> Purpose: compute advection
!----------------------------------------------------------------------
subroutine compute_adv(bdata,ndata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata    !< B data
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: its,il0,ic0,ic0a,nc0a,nc0d,i_s,i_s_loc,ic0d,iproc,jc0,d_n_s_max,d_n_s_max_loc
integer,allocatable :: c0d_to_c0(:),c0_to_c0d(:),c0a_to_c0d(:),interpd_lg(:,:,:)
integer,allocatable :: c0d_to_c0_copy(:),c0a_to_c0d_copy(:)
logical :: mask_c0(ndata%geom%nc0)
logical,allocatable :: lcheck_c0d(:),lcheck_d(:,:,:)
type(linoptype),allocatable :: dfull(:,:)
type(comtype) :: comAD(mpl%nproc)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

write(mpl%unit,'(a7,a)') '','Compute advection'

! Allocation
allocate(dfull(geom%nl0,nam%nts-1))
allocate(ndata%d(geom%nl0,nam%nts-1))
allocate(lcheck_c0d(geom%nc0))

! Initialization
mask_c0 = .true.

do its=2,nam%nts
   do il0=1,geom%nl0
      ! Compute interpolation
      call compute_interp(geom%nc0,bdata%lon_c0_flt(:,il0,its),bdata%lat_c0_flt(:,il0,its),mask_c0,geom%nc0,geom%lon,geom%lat, &
    & mask_c0,nam%diag_interp,dfull(il0,its-1))
   end do
end do

! Allocation
d_n_s_max = 0
do its=2,nam%nts
   do il0=1,geom%nl0
      d_n_s_max = max(d_n_s_max,dfull(il0,its-1)%n_s)
   end do
end do
allocate(lcheck_d(d_n_s_max,geom%nl0,2:nam%nts))

! Halo definitions

! Halo A
lcheck_c0d = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   lcheck_c0d(ic0) = .true.
end do

! Halo B
lcheck_d = .false.
do its=2,nam%nts
   do il0=1,geom%nl0
      do i_s=1,dfull(il0,its-1)%n_s
         ic0 = dfull(il0,its-1)%row(i_s)
         iproc = geom%c0_to_proc(ic0)
         if (iproc==mpl%myproc) then
            jc0 = dfull(il0,its-1)%col(i_s)
            lcheck_d(i_s,il0,its) = .true.
            lcheck_c0d(jc0) = .true.
         end if
      end do
   end do
end do

! Halo sizes
do its=2,nam%nts
   do il0=1,geom%nl0
      ndata%d(il0,its-1)%n_s = count(lcheck_d(:,il0,its))
   end do
end do
ndata%nc0d = count(lcheck_c0d)

! Global <-> local conversions for fields
allocate(c0d_to_c0(ndata%nc0d))
allocate(c0_to_c0d(geom%nc0))
call msi(c0_to_c0d)
ic0d = 0
do ic0=1,geom%nc0
   if (lcheck_c0d(ic0)) then
      ic0d = ic0d+1
      c0d_to_c0(ic0d) = ic0
      c0_to_c0d(ic0) = ic0d
   end if
end do

! Inter-halo conversions
allocate(c0a_to_c0d(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0d = c0_to_c0d(ic0)
   c0a_to_c0d(ic0a) = ic0d
end do

! Global <-> local conversions for data
d_n_s_max_loc = 0
do its=2,nam%nts
   do il0=1,geom%nl0
      d_n_s_max_loc = max(d_n_s_max_loc,ndata%d(il0,its-1)%n_s)
   end do
end do
allocate(interpd_lg(d_n_s_max_loc,geom%nl0,2:nam%nts))
do its=2,nam%nts
   do il0=1,geom%nl0
      i_s_loc = 0
      do i_s=1,dfull(il0,its-1)%n_s
         if (lcheck_d(i_s,il0,its)) then
            i_s_loc = i_s_loc+1
            interpd_lg(i_s_loc,il0,its) = i_s
         end if
      end do
   end do
end do

! Local data
do its=2,nam%nts
   do il0=1,geom%nl0
      ndata%d(il0,its-1)%prefix = 'd'
      ndata%d(il0,its-1)%n_src = ndata%nc0d
      ndata%d(il0,its-1)%n_dst = geom%nc0a
      call linop_alloc(ndata%d(il0,its-1))
      do i_s_loc=1,ndata%d(il0,its-1)%n_s
         i_s = interpd_lg(i_s_loc,il0,its)
         ndata%d(il0,its-1)%row(i_s_loc) = geom%c0_to_c0a(dfull(il0,its-1)%row(i_s))
         ndata%d(il0,its-1)%col(i_s_loc) = c0_to_c0d(dfull(il0,its-1)%col(i_s))
         ndata%d(il0,its-1)%S(i_s_loc) = dfull(il0,its-1)%S(i_s)
      end do
      call linop_reorder(ndata%d(il0,its-1))
   end do
end do

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = geom%nc0a
         nc0d = ndata%nc0d
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
         call mpl_recv(nc0d,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(c0d_to_c0_copy(nc0d))
      allocate(c0a_to_c0d_copy(nc0a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0d_to_c0_copy = c0d_to_c0
         c0a_to_c0d_copy = c0a_to_c0d
      else
         ! Receive data on ioproc
         call mpl_recv(nc0d,c0d_to_c0_copy,iproc,mpl%tag+2)
         call mpl_recv(nc0a,c0a_to_c0d_copy,iproc,mpl%tag+3)
      end if

      ! Allocation
      comAD(iproc)%nred = nc0a
      comAD(iproc)%next = nc0d
      allocate(comAD(iproc)%ext_to_proc(comAD(iproc)%next))
      allocate(comAD(iproc)%ext_to_red(comAD(iproc)%next))
      allocate(comAD(iproc)%red_to_ext(comAD(iproc)%nred))

      ! AD communication
      do ic0d=1,nc0d
         ic0 = c0d_to_c0_copy(ic0d)
         comAD(iproc)%ext_to_proc(ic0d) = geom%c0_to_proc(ic0)
         ic0a = geom%c0_to_c0a(ic0)
         comAD(iproc)%ext_to_red(ic0d) = ic0a
      end do
      comAD(iproc)%red_to_ext = c0a_to_c0d_copy

      ! Release memory
      deallocate(c0d_to_c0_copy)
      deallocate(c0a_to_c0d_copy)
   end do

   ! Communication setup
   call com_setup(comAD)
else
   ! Send dimensions to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)
   call mpl_send(ndata%nc0d,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl_send(ndata%nc0d,c0d_to_c0,mpl%ioproc,mpl%tag+2)
   call mpl_send(geom%nc0a,c0a_to_c0d,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4

! Communication broadcast
ndata%AD%prefix = 'AD'
call com_bcast(comAD,ndata%AD)

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0d =       ',ndata%nc0d
do its=2,nam%nts
   do il0=1,geom%nl0
      write(mpl%unit,'(a10,a,i3,a,i2,a,i8)') '','d(',il0,',',its,')%n_s = ',ndata%d(il0,its-1)%n_s
   end do
end do

! End associate
end associate

end subroutine compute_adv

end module nicas_parameters_adv
