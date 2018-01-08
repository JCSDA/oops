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
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use type_com, only: comtype,com_dealloc,com_setup,com_bcast
use type_hdata, only: hdatatype
use type_linop, only: linop_alloc,linop_copy,linop_reorder
use type_mpl, only: mpl,mpl_send,mpl_recv
use type_nam, only: namtype

implicit none

private
public :: compute_mpi

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi
!> Purpose: compute HDIAG MPI distribution
!----------------------------------------------------------------------
subroutine compute_mpi(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: iproc,jc3,ic0,ic0a,ic0c,ic0d,ic1,ic1a,ic2,ic2a,jc1,jc0
integer :: nc0a,nc0c,nc0d
integer,allocatable :: c0_to_c0a(:),c0a_to_c0(:),c0c_to_c0(:),c0a_to_c0c(:),c0d_to_c0(:),c0a_to_c0d(:)
type(comtype) :: comAC(mpl%nproc),comAD(mpl%nproc)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Allocation
allocate(hdata%lcheck_c0a(geom%nc0))
allocate(hdata%lcheck_c1a(nam%nc1))
allocate(hdata%lcheck_c0c(geom%nc0))
if (nam%displ_diag) allocate(hdata%lcheck_c0d(geom%nc0))

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

! Halo C
hdata%lcheck_c0c = hdata%lcheck_c0a
do jc3=1,nam%nc3
   do ic1=1,nam%nc1
      if (any(hdata%c1c3l0_log(ic1,jc3,:))) then
         ic0 = hdata%c1c3_to_c0(ic1,jc3)
         if (any(geom%mask(ic0,:))) hdata%lcheck_c0c(ic0) = .true.
      end if
   end do
end do
hdata%nc0c = count(hdata%lcheck_c0c)

if (nam%displ_diag) then
   ! Halo D
   hdata%lcheck_c0d = hdata%lcheck_c0a
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      if (any(geom%mask(ic0,:).and.(geom%c0_to_proc(ic0)==mpl%myproc))) hdata%lcheck_c0d(ic0) = .true.
   end do
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
end if

! Global <-> local conversions for fields

! Halo A
allocate(hdata%c1a_to_c1(hdata%nc1a))
ic1a = 0
do ic1=1,nam%nc1
   if (hdata%lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      hdata%c1a_to_c1(ic1a) = ic1
   end if
end do

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

if (nam%displ_diag) then
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
end if

! Inter-halo conversions
allocate(hdata%c0a_to_c0c(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0c = hdata%c0_to_c0c(ic0)
   hdata%c0a_to_c0c(ic0a) = ic0c
end do
if (nam%displ_diag) then
   allocate(hdata%c0a_to_c0d(geom%nc0a))
   do ic0a=1,geom%nc0a
      ic0 = geom%c0a_to_c0(ic0a)
      ic0d = hdata%c0_to_c0d(ic0)
      hdata%c0a_to_c0d(ic0a) = ic0d
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
         if (nam%displ_diag) nc0d = hdata%nc0d
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nc0a,iproc,mpl%tag)
         call mpl_recv(nc0c,iproc,mpl%tag+1)
         if (nam%displ_diag) call mpl_recv(nc0d,iproc,mpl%tag+2)
      end if

      ! Allocation
      allocate(c0c_to_c0(nc0c))
      allocate(c0a_to_c0c(nc0a))
      if (nam%displ_diag) then
         allocate(c0d_to_c0(nc0d))
         allocate(c0a_to_c0d(nc0a))
      end if

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         c0c_to_c0 = hdata%c0c_to_c0
         c0a_to_c0c = hdata%c0a_to_c0c
         if (nam%displ_diag) then
            c0d_to_c0 = hdata%c0d_to_c0
            c0a_to_c0d = hdata%c0a_to_c0d
         end if
      else
         ! Receive data on ioproc
         call mpl_recv(nc0c,c0c_to_c0,iproc,mpl%tag+3)
         call mpl_recv(nc0a,c0a_to_c0c,iproc,mpl%tag+4)
         if (nam%displ_diag) then
            call mpl_recv(nc0d,c0d_to_c0,iproc,mpl%tag+5)
            call mpl_recv(nc0a,c0a_to_c0d,iproc,mpl%tag+6)
         end if
      end if

      ! Allocation
      comAC(iproc)%nred = nc0a
      comAC(iproc)%next = nc0c
      allocate(comAC(iproc)%ext_to_proc(comAC(iproc)%next))
      allocate(comAC(iproc)%ext_to_red(comAC(iproc)%next))
      allocate(comAC(iproc)%red_to_ext(comAC(iproc)%nred))
      if (nam%displ_diag) then
         comAD(iproc)%nred = nc0a
         comAD(iproc)%next = nc0d
         allocate(comAD(iproc)%ext_to_proc(comAD(iproc)%next))
         allocate(comAD(iproc)%ext_to_red(comAD(iproc)%next))
         allocate(comAD(iproc)%red_to_ext(comAD(iproc)%nred))
      end if

      ! AC communication
      do ic0c=1,nc0c
         ic0 = c0c_to_c0(ic0c)
         comAC(iproc)%ext_to_proc(ic0c) = geom%c0_to_proc(ic0)
         ic0a = c0_to_c0a(ic0)
         comAC(iproc)%ext_to_red(ic0c) = ic0a
      end do
      comAC(iproc)%red_to_ext = c0a_to_c0c

      if (nam%displ_diag) then
         ! AD communication
         do ic0d=1,nc0d
            ic0 = c0d_to_c0(ic0d)
            comAD(iproc)%ext_to_proc(ic0d) = geom%c0_to_proc(ic0)
            ic0a = c0_to_c0a(ic0)
            comAD(iproc)%ext_to_red(ic0d) = ic0a
         end do
         comAD(iproc)%red_to_ext = c0a_to_c0d
      end if

      ! Release memory
      deallocate(c0c_to_c0)
      deallocate(c0a_to_c0c)
      if (nam%displ_diag) then
         deallocate(c0d_to_c0)
         deallocate(c0a_to_c0d)
      end if
   end do

   ! Communication setup
   call com_setup(comAC)
   if (nam%displ_diag) call com_setup(comAD)
 
   ! Release memory
   deallocate(c0_to_c0a)
   if (nam%displ_diag) deallocate(c0_to_c0a)
else
   ! Send dimensions to ioproc
   call mpl_send(geom%nc0a,mpl%ioproc,mpl%tag)
   call mpl_send(hdata%nc0c,mpl%ioproc,mpl%tag+1)
   if (nam%displ_diag) call mpl_send(hdata%nc0d,mpl%ioproc,mpl%tag+2)

   ! Send data to ioproc
   call mpl_send(hdata%nc0c,hdata%c0c_to_c0,mpl%ioproc,mpl%tag+3)
   call mpl_send(geom%nc0a,hdata%c0a_to_c0c,mpl%ioproc,mpl%tag+4)
   if (nam%displ_diag) then
      call mpl_send(hdata%nc0d,hdata%c0d_to_c0,mpl%ioproc,mpl%tag+5)
      call mpl_send(geom%nc0a,hdata%c0a_to_c0d,mpl%ioproc,mpl%tag+6)
   end if
end if
mpl%tag = mpl%tag+7

! Communication broadcast
hdata%AC%prefix = 'AC'
call com_bcast(comAC,hdata%AC)
if (nam%displ_diag) then
   hdata%AD%prefix = 'AD'
   call com_bcast(comAD,hdata%AD)
end if

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0 =       ',geom%nc0
write(mpl%unit,'(a10,a,i8)') '','nc0a =      ',geom%nc0a
write(mpl%unit,'(a10,a,i8)') '','nc0c =      ',hdata%nc0c
if (nam%displ_diag) write(mpl%unit,'(a10,a,i8)') '','nc0d =      ',hdata%nc0d
write(mpl%unit,'(a10,a,i8)') '','nl0 =       ',geom%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =       ',nam%nc1
write(mpl%unit,'(a10,a,i8)') '','nc1a =      ',hdata%nc1a

! End associate
end associate

end subroutine compute_mpi

end module hdiag_mpi
