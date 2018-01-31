!----------------------------------------------------------------------
! Module: hdiag_mpi.f90
!> Purpose: compute HDIAG MPI distribution
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_mpi

use tools_display, only: msgerror,prog_init,prog_print
use tools_interp, only: compute_interp
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use type_com, only: comtype,com_dealloc,com_setup,com_bcast
use type_displ, only: displtype
use type_hdata, only: hdatatype
use type_linop, only: linoptype,linop_alloc,linop_copy,linop_reorder
use type_mpl, only: mpl,mpl_send,mpl_recv
use type_nam, only: namtype

implicit none

private
public :: compute_mpi_a,compute_mpi_d,compute_mpi_c

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi_a
!> Purpose: compute HDIAG MPI distribution, halo A
!----------------------------------------------------------------------
subroutine compute_mpi_a(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: ic0a,ic0,ic1a,ic1

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Allocation
allocate(hdata%lcheck_c0a(geom%nc0))
allocate(hdata%lcheck_c1a(nam%nc1))

! Halo definitions

! Halo A
hdata%lcheck_c0a = .false.
hdata%lcheck_c1a = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   if (any(geom%mask(ic0,:).and.(geom%c0_to_proc(ic0)==mpl%myproc))) hdata%lcheck_c0a(ic0) = .true.
end do
do ic1=1,nam%nc1
   ic0 = hdata%c1_to_c0(ic1)
   if (any(geom%mask(ic0,:).and.(geom%c0_to_proc(ic0)==mpl%myproc))) hdata%lcheck_c1a(ic1) = .true.
end do
hdata%nc1a = count(hdata%lcheck_c1a)

! Global <-> local conversions for fields

! Halo A
allocate(hdata%c1a_to_c1(hdata%nc1a))
allocate(hdata%c1_to_c1a(nam%nc1))
ic1a = 0
do ic1=1,nam%nc1
   if (hdata%lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      hdata%c1a_to_c1(ic1a) = ic1
      hdata%c1_to_c1a(ic1) = ic1a
   end if
end do

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0 =       ',geom%nc0
write(mpl%unit,'(a10,a,i8)') '','nc0a =      ',geom%nc0a
write(mpl%unit,'(a10,a,i8)') '','nl0 =       ',geom%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =       ',nam%nc1
write(mpl%unit,'(a10,a,i8)') '','nc1a =      ',hdata%nc1a

! End associate
end associate

end subroutine compute_mpi_a

!----------------------------------------------------------------------
! Subroutine: compute_mpi_d
!> Purpose: compute HDIAG MPI distribution, halo D
!----------------------------------------------------------------------
subroutine compute_mpi_d(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: iproc,ic0,ic0a,ic0d,ic1,ic2,ic2a,jc0,jc1
integer :: nc0a,nc0d
integer,allocatable :: c0_to_c0a(:),c0a_to_c0(:),c0d_to_c0(:),c0a_to_c0d(:)
type(comtype) :: comAD(mpl%nproc)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Allocation
allocate(hdata%lcheck_c0d(geom%nc0))

! Halo definitions

! Halo D
hdata%lcheck_c0d = hdata%lcheck_c0a
do ic2a=1,hdata%nc2a
   ic2 = hdata%c2a_to_c2(ic2a)
   ic1 = hdata%c2_to_c1(ic2)
   if (any(hdata%c1l0_log(ic1,:))) then
      do jc1=1,nam%nc1
         if (any(hdata%displ_mask(jc1,ic2,:))) then
            jc0 = hdata%c1_to_c0(jc1)
            hdata%lcheck_c0d(jc0) = .true.
         end if
      end do
   end if
end do
hdata%nc0d = count(hdata%lcheck_c0d)

! Global <-> local conversions for fields

! Halo D
allocate(hdata%c0d_to_c0(hdata%nc0d))
allocate(hdata%c0_to_c0d(geom%nc0))
call msi(hdata%c0_to_c0d)
ic0d = 0
do ic0=1,geom%nc0
   if (hdata%lcheck_c0d(ic0)) then
      ic0d = ic0d+1
      hdata%c0d_to_c0(ic0d) = ic0
      hdata%c0_to_c0d(ic0) = ic0d
   end if
end do

! Inter-halo conversions
allocate(hdata%c0a_to_c0d(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0d = hdata%c0_to_c0d(ic0)
   hdata%c0a_to_c0d(ic0a) = ic0d
end do

! Get global distribution of the subgrid on ioproc
if (mpl%main) then
   ! Allocation
   allocate(c0_to_c0a(geom%nc0))

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         nc0a = geom%nc0a
      else
         ! Receive dimension on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
      end if

      ! Allocation
      allocate(c0a_to_c0(nc0a))

      if (iproc==mpl%ioproc) then
         ! Copy data
         c0a_to_c0 = geom%c0a_to_c0
      else
         ! Receive data on ioproc
         call mpl_recv(nc0a,c0a_to_c0,iproc,mpl%tag+1)
      end if

      ! Fill c0_to_c0a
      do ic0a=1,nc0a
         c0_to_c0a(c0a_to_c0(ic0a)) = ic0a
      end do

      ! Release memory
      deallocate(c0a_to_c0)
   end do
else
   ! Send dimensions to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)

   ! Send data to ioproc
   call mpl_send(geom%nc0a,geom%c0a_to_c0,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = geom%nc0a
         nc0d = hdata%nc0d
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
         call mpl_recv(nc0d,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(c0d_to_c0(nc0d))
      allocate(c0a_to_c0d(nc0a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0d_to_c0 = hdata%c0d_to_c0
         c0a_to_c0d = hdata%c0a_to_c0d
      else
         ! Receive data on ioproc
         call mpl_recv(nc0d,c0d_to_c0,iproc,mpl%tag+2)
         call mpl_recv(nc0a,c0a_to_c0d,iproc,mpl%tag+3)
      end if

      ! Allocation
      comAD(iproc)%nred = nc0a
      comAD(iproc)%next = nc0d
      allocate(comAD(iproc)%ext_to_proc(comAD(iproc)%next))
      allocate(comAD(iproc)%ext_to_red(comAD(iproc)%next))
      allocate(comAD(iproc)%red_to_ext(comAD(iproc)%nred))

      ! AD communication
      do ic0d=1,nc0d
         ic0 = c0d_to_c0(ic0d)
         comAD(iproc)%ext_to_proc(ic0d) = geom%c0_to_proc(ic0)
         ic0a = c0_to_c0a(ic0)
         comAD(iproc)%ext_to_red(ic0d) = ic0a
      end do
      comAD(iproc)%red_to_ext = c0a_to_c0d

      ! Release memory
      deallocate(c0d_to_c0)
      deallocate(c0a_to_c0d)
   end do

   ! Communication setup
   call com_setup(comAD)

   ! Release memory
   deallocate(c0_to_c0a)
else
   ! Send dimensions to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)
   call mpl_send(hdata%nc0d,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl_send(hdata%nc0d,hdata%c0d_to_c0,mpl%ioproc,mpl%tag+2)
   call mpl_send(geom%nc0a,hdata%c0a_to_c0d,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4

! Communication broadcast
hdata%AD%prefix = 'AD'
call com_bcast(comAD,hdata%AD)

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0d =      ',hdata%nc0d

! End associate
end associate

end subroutine compute_mpi_d

!----------------------------------------------------------------------
! Subroutine: compute_mpi_c
!> Purpose: compute HDIAG MPI distribution, halo C
!----------------------------------------------------------------------
subroutine compute_mpi_c(hdata,displ)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data
type(displtype),intent(in) :: displ    !< Displacement

! Local variables
integer :: iproc,jc3,ic0,ic0a,ic0c,ic1,ic1a,its,il0,d_n_s_max,d_n_s_max_loc,i_s,i_s_loc
integer :: nc0a,nc0c
integer,allocatable :: interpd_lg(:,:,:),c0_to_c0a(:),c0a_to_c0(:),c0c_to_c0(:),c0a_to_c0c(:)
real(kind_real),allocatable :: lon_c1(:),lat_c1(:)
type(comtype) :: comAC(mpl%nproc)
type(linoptype),allocatable :: dfull(:,:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

if (nam%displ_diag) then
   ! Allocation
   allocate(dfull(geom%nl0,nam%nts))
   allocate(lon_c1(nam%nc1))
   allocate(lat_c1(nam%nc1))

   ! Prepare displacement interpolation
   do its=1,nam%nts
      do il0=1,geom%nl0
         ! Copy Sc1 points
         do ic1=1,nam%nc1
            ic0 = hdata%c1_to_c0(ic1)
            lon_c1(ic1) = displ%lon_c0_flt(ic0,il0,its)
            lat_c1(ic1) = displ%lat_c0_flt(ic0,il0,its)
         end do

         ! Compute interpolation
         call compute_interp(geom%mesh,geom%ctree,geom%nc0,geom%mask(:,il0),nam%nc1,lon_c1,lat_c1,hdata%c1l0_log(:,il0), &
       & nam%diag_interp,dfull(il0,its))
      end do
   end do
end if

! Allocation
if (nam%displ_diag) then
   d_n_s_max = 0
   do its=1,nam%nts
      do il0=1,geom%nl0
         d_n_s_max = max(d_n_s_max,dfull(il0,its)%n_s)
      end do
   end do
   allocate(hdata%lcheck_d(d_n_s_max,geom%nl0,nam%nts))
   allocate(hdata%d(geom%nl0,nam%nts))
end if
allocate(hdata%lcheck_c0c(geom%nc0))

! Halo C
hdata%lcheck_c0c = hdata%lcheck_c0a
do jc3=1,nam%nc3
   do ic1a=1,hdata%nc1a
      ic1 = hdata%c1a_to_c1(ic1a)
      if (any(hdata%c1c3l0_log(ic1,jc3,:))) then
         ic0 = hdata%c1c3_to_c0(ic1,jc3)
         if (any(geom%mask(ic0,:))) hdata%lcheck_c0c(ic0) = .true.
      end if
   end do
end do
if (nam%displ_diag) then
   hdata%lcheck_d = .false.
   do its=1,nam%nts
      do il0=1,geom%nl0
         do i_s=1,dfull(il0,its)%n_s
            ic0 = dfull(il0,its)%col(i_s)
            ic1 = dfull(il0,its)%row(i_s)
            if (hdata%lcheck_c1a(ic1)) then
               hdata%lcheck_c0c(ic0) = .true.
               hdata%lcheck_d(i_s,il0,its) = .true.
            end if
         end do
      end do
   end do
end if
hdata%nc0c = count(hdata%lcheck_c0c)
if (nam%displ_diag) then
   do its=1,nam%nts
      do il0=1,geom%nl0
         hdata%d(il0,its)%n_s = count(hdata%lcheck_d(:,il0,its))
      end do
   end do
end if

! Global <-> local conversions for fields

! Halo C
allocate(hdata%c0c_to_c0(hdata%nc0c))
allocate(hdata%c0_to_c0c(geom%nc0))
call msi(hdata%c0_to_c0c)
ic0c = 0
do ic0=1,geom%nc0
   if (hdata%lcheck_c0c(ic0)) then
      ic0c = ic0c+1
      hdata%c0c_to_c0(ic0c) = ic0
      hdata%c0_to_c0c(ic0) = ic0c
   end if
end do

! Inter-halo conversions
allocate(hdata%c0a_to_c0c(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0c = hdata%c0_to_c0c(ic0)
   hdata%c0a_to_c0c(ic0a) = ic0c
end do

if (nam%displ_diag) then
   ! Global <-> local conversions for data
   d_n_s_max_loc = 0
   do its=1,nam%nts
      do il0=1,geom%nl0
         d_n_s_max_loc = max(d_n_s_max_loc,hdata%d(il0,its)%n_s)
      end do
   end do
   allocate(interpd_lg(d_n_s_max_loc,geom%nl0,nam%nts))
   do its=1,nam%nts
      do il0=1,geom%nl0
         i_s_loc = 0
         do i_s=1,dfull(il0,its)%n_s
            if (hdata%lcheck_d(i_s,il0,its)) then
               i_s_loc = i_s_loc+1
               interpd_lg(i_s_loc,il0,its) = i_s
            end if
         end do
      end do
   end do

   ! Local data
   do its=1,nam%nts
      do il0=1,geom%nl0
         hdata%d(il0,its)%prefix = 'd'
         hdata%d(il0,its)%n_src = hdata%nc0c
         hdata%d(il0,its)%n_dst = hdata%nc1a
         call linop_alloc(hdata%d(il0,its))
         do i_s_loc=1,hdata%d(il0,its)%n_s
            i_s = interpd_lg(i_s_loc,il0,its)
            hdata%d(il0,its)%row(i_s_loc) = hdata%c1_to_c1a(dfull(il0,its)%row(i_s))
            hdata%d(il0,its)%col(i_s_loc) = hdata%c0_to_c0c(dfull(il0,its)%col(i_s))
            hdata%d(il0,its)%S(i_s_loc) = dfull(il0,its)%S(i_s)
         end do
         call linop_reorder(hdata%d(il0,its))
      end do
   end do
end if

! Get global distribution of the subgrid on ioproc
if (mpl%main) then
   ! Allocation
   allocate(c0_to_c0a(geom%nc0))

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         nc0a = geom%nc0a
      else
         ! Receive dimension on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
      end if

      ! Allocation
      allocate(c0a_to_c0(nc0a))

      if (iproc==mpl%ioproc) then
         ! Copy data
         c0a_to_c0 = geom%c0a_to_c0
      else
         ! Receive data on ioproc
         call mpl_recv(nc0a,c0a_to_c0,iproc,mpl%tag+1)
      end if

      ! Fill c0_to_c0a
      do ic0a=1,nc0a
         c0_to_c0a(c0a_to_c0(ic0a)) = ic0a
      end do

      ! Release memory
      deallocate(c0a_to_c0)
   end do
else
   ! Send dimensions to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)

   ! Send data to ioproc
   call mpl_send(geom%nc0a,geom%c0a_to_c0,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nc0a = geom%nc0a
         nc0c = hdata%nc0c
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
         call mpl_recv(nc0c,iproc,mpl%tag+1)
      end if

      ! Allocation
      allocate(c0c_to_c0(nc0c))
      allocate(c0a_to_c0c(nc0a))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0c_to_c0 = hdata%c0c_to_c0
         c0a_to_c0c = hdata%c0a_to_c0c
      else
         ! Receive data on ioproc
         call mpl_recv(nc0c,c0c_to_c0,iproc,mpl%tag+2)
         call mpl_recv(nc0a,c0a_to_c0c,iproc,mpl%tag+3)
      end if

      ! Allocation
      comAC(iproc)%nred = nc0a
      comAC(iproc)%next = nc0c
      allocate(comAC(iproc)%ext_to_proc(comAC(iproc)%next))
      allocate(comAC(iproc)%ext_to_red(comAC(iproc)%next))
      allocate(comAC(iproc)%red_to_ext(comAC(iproc)%nred))

      ! AC communication
      do ic0c=1,nc0c
         ic0 = c0c_to_c0(ic0c)
         comAC(iproc)%ext_to_proc(ic0c) = geom%c0_to_proc(ic0)
         ic0a = c0_to_c0a(ic0)
         comAC(iproc)%ext_to_red(ic0c) = ic0a
      end do
      comAC(iproc)%red_to_ext = c0a_to_c0c

      ! Release memory
      deallocate(c0c_to_c0)
      deallocate(c0a_to_c0c)
   end do

   ! Communication setup
   call com_setup(comAC)

   ! Release memory
   deallocate(c0_to_c0a)
else
   ! Send dimensions to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)
   call mpl_send(hdata%nc0c,mpl%ioproc,mpl%tag+1)

   ! Send data to ioproc
   call mpl_send(hdata%nc0c,hdata%c0c_to_c0,mpl%ioproc,mpl%tag+2)
   call mpl_send(geom%nc0a,hdata%c0a_to_c0c,mpl%ioproc,mpl%tag+3)
end if
mpl%tag = mpl%tag+4

! Communication broadcast
hdata%AC%prefix = 'AC'
call com_bcast(comAC,hdata%AC)

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0c =      ',hdata%nc0c

! End associate
end associate

end subroutine compute_mpi_c

end module hdiag_mpi
