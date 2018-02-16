!----------------------------------------------------------------------
! Module: nicas_parameters_mpi.f90
!> Purpose: compute NICAS parameters MPI distribution
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_parameters_mpi

use netcdf
use omp_lib
use tools_const, only: pi,rad2deg,req,sphere_dist
use tools_display, only: msgerror,prog_init,prog_print
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_com, only: comtype,com_setup,com_bcast
use type_linop, only: linop_alloc,linop_reorder
use type_mpl, only: mpl,mpl_send,mpl_recv
use type_nam, only: namtype
use type_ndata, only: ndatatype

implicit none

private
public :: compute_mpi_ab,compute_mpi_c

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi_ab
!> Purpose: compute NICAS MPI distribution, halos A-B
!----------------------------------------------------------------------
subroutine compute_mpi_ab(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: il0i,ic0,iproc,ic1,jc1,ic1a,ic1b,il0,il1,isa,isb,i_s,is,js,h_n_s_max,s_n_s_max
logical :: lcheck_c1b_h(ndata%nc1)
! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Allocation
h_n_s_max = 0
do il0i=1,geom%nl0i
   h_n_s_max = max(h_n_s_max,ndata%hfull(il0i)%n_s)
end do
s_n_s_max = 0
do il1=1,ndata%nl1
   s_n_s_max = max(s_n_s_max,ndata%sfull(il1)%n_s)
end do
allocate(ndata%h(geom%nl0i))
allocate(ndata%s(ndata%nl1))
allocate(ndata%lcheck_c1a(ndata%nc1))
allocate(ndata%lcheck_c1b(ndata%nc1))
allocate(ndata%lcheck_sa(ndata%ns))
allocate(ndata%lcheck_sb(ndata%ns))
allocate(ndata%lcheck_h(h_n_s_max,geom%nl0i))
allocate(ndata%lcheck_s(s_n_s_max,ndata%nl1))

! Halo definitions

! Halo A
ndata%lcheck_c1a = .false.
ndata%lcheck_sa = .false.
do is=1,ndata%ns
   ic1 = ndata%s_to_c1(is)
   ic0 = ndata%c1_to_c0(ic1)
   il1 = ndata%s_to_l1(is)
   il0 = ndata%l1_to_l0(il1)
   if (geom%c0_to_proc(ic0)==mpl%myproc) then
      ndata%lcheck_c1a(ic1) = .true.
      ndata%lcheck_sa(is) = .true.
   end if
end do

! Halo B

! Horizontal interpolation
ndata%lcheck_h = .false.
lcheck_c1b_h = .false.
do il0i=1,geom%nl0i
   do i_s=1,ndata%hfull(il0i)%n_s
      ic0 = ndata%hfull(il0i)%row(i_s)
      iproc = geom%c0_to_proc(ic0)
      if (iproc==mpl%myproc) then
         jc1 = ndata%hfull(il0i)%col(i_s)
         ndata%lcheck_h(i_s,il0i) = .true.
         lcheck_c1b_h(jc1) = .true.
      end if
   end do
end do

! Subsampling horizontal interpolation
ndata%lcheck_sb = .false.
ndata%lcheck_s = .false.
ndata%lcheck_c1b = lcheck_c1b_h
do il1=1,ndata%nl1
   do i_s=1,ndata%sfull(il1)%n_s
      ic1 = ndata%sfull(il1)%row(i_s)
      if (lcheck_c1b_h(ic1)) then
         jc1 = ndata%sfull(il1)%col(i_s)
         js = ndata%c1l1_to_s(jc1,il1)
         ndata%lcheck_c1b(jc1) = .true.
         ndata%lcheck_sb(js) = .true.
         ndata%lcheck_s(i_s,il1) = .true.
      end if
   end do
end do

! Check halos consistency
do is=1,ndata%ns
   if (ndata%lcheck_sa(is).and.(.not.ndata%lcheck_sb(is))) call msgerror('point in halo A but not in halo B')
end do

! Halo sizes
ndata%nc1a = count(ndata%lcheck_c1a)
ndata%nsa = count(ndata%lcheck_sa)
do il0i=1,geom%nl0i
   ndata%h(il0i)%n_s = count(ndata%lcheck_h(:,il0i))
end do
ndata%nc1b = count(ndata%lcheck_c1b)
ndata%nsb = count(ndata%lcheck_sb)
do il1=1,ndata%nl1
   ndata%s(il1)%n_s = count(ndata%lcheck_s(:,il1))
end do

! Global <-> local conversions for fields

! Halo A
allocate(ndata%c1a_to_c1(ndata%nc1a))
allocate(ndata%c1_to_c1a(ndata%nc1))
call msi(ndata%c1_to_c1a)
ic1a = 0
do ic1=1,ndata%nc1
   if (ndata%lcheck_c1a(ic1)) then
      ic1a = ic1a+1
      ndata%c1a_to_c1(ic1a) = ic1
      ndata%c1_to_c1a(ic1) = ic1a
   end if
end do

allocate(ndata%sa_to_s(ndata%nsa))
isa = 0
do is=1,ndata%ns
   if (ndata%lcheck_sa(is)) then
      isa = isa+1
      ndata%sa_to_s(isa) = is
   end if
end do

! Halo B
allocate(ndata%c1b_to_c1(ndata%nc1b))
allocate(ndata%c1_to_c1b(ndata%nc1))
call msi(ndata%c1_to_c1b)
ic1b = 0
do ic1=1,ndata%nc1
   if (ndata%lcheck_c1b(ic1)) then
      ic1b = ic1b+1
      ndata%c1b_to_c1(ic1b) = ic1
      ndata%c1_to_c1b(ic1) = ic1b
   end if
end do

allocate(ndata%sb_to_s(ndata%nsb))
allocate(ndata%s_to_sb(ndata%ns))
call msi(ndata%s_to_sb)
isb = 0
do is=1,ndata%ns
   if (ndata%lcheck_sb(is)) then
      isb = isb+1
      ndata%sb_to_s(isb) = is
      ndata%s_to_sb(is) = isb
   end if
end do

! End associate
end associate

end subroutine compute_mpi_ab

!----------------------------------------------------------------------
! Subroutine: compute_mpi_c
!> Purpose: compute NICAS MPI distribution, halo C
!----------------------------------------------------------------------
subroutine compute_mpi_c(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: iproc,il0i,ic0,ic1,ic1b,il1,isa,isb,isc,i_s,i_s_loc,is,js,h_n_s_max_loc,s_n_s_max_loc
integer :: nsa,nsb,nsc
integer,allocatable :: interph_lg(:,:),interps_lg(:,:)
integer,allocatable :: s_to_sa(:),sa_to_s(:),sb_to_s(:),sc_to_s(:),sa_to_sb(:),sa_to_sc(:)
type(comtype) :: comAB(mpl%nproc),comAC(mpl%nproc)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Allocation
allocate(ndata%lcheck_sc(ndata%ns))
allocate(ndata%lcheck_sc_nor(ndata%ns))

! Halo definitions

! Halo C
if (nam%mpicom==1) then
   ! 1 communication step
   ndata%lcheck_sc = ndata%lcheck_sb
   do i_s=1,ndata%c%n_s
      is = ndata%c%row(i_s)
      js = ndata%c%col(i_s)
      ndata%lcheck_sc(is) = .true.
      ndata%lcheck_sc(js) = .true.
   end do
elseif (nam%mpicom==2) then
   ! 2 communication steps
   ndata%lcheck_sc = ndata%lcheck_sb
   do i_s=1,ndata%c%n_s
      is = ndata%c%row(i_s)
      js = ndata%c%col(i_s)
      ndata%lcheck_sc(is) = .true.
      ndata%lcheck_sc(js) = .true.
   end do
end if
ndata%nsc = count(ndata%lcheck_sc)
ndata%lcheck_sc_nor = ndata%lcheck_sb
do i_s=1,ndata%c_nor%n_s
   is = ndata%c_nor%row(i_s)
   js = ndata%c_nor%col(i_s)
   ndata%lcheck_sc_nor(is) = .true.
   ndata%lcheck_sc_nor(js) = .true.
end do
ndata%nsc_nor = count(ndata%lcheck_sc_nor)

! Check halos consistency
do is=1,ndata%ns
   if (ndata%lcheck_sa(is).and.(.not.ndata%lcheck_sc(is))) call msgerror('point in halo A but not in halo C')
   if (ndata%lcheck_sb(is).and.(.not.ndata%lcheck_sc(is))) call msgerror('point in halo B but not in halo C')
end do

! Global <-> local conversions for fields

! Halo C
allocate(ndata%sc_to_s(ndata%nsc))
allocate(ndata%s_to_sc(ndata%ns))
call msi(ndata%s_to_sc)
isc = 0
do is=1,ndata%ns
   if (ndata%lcheck_sc(is)) then
      isc = isc+1
      ndata%sc_to_s(isc) = is
      ndata%s_to_sc(is) = isc
   end if
end do

allocate(ndata%sc_nor_to_s(ndata%nsc_nor))
allocate(ndata%s_to_sc_nor(ndata%ns))
call msi(ndata%s_to_sc_nor)
isc = 0
do is=1,ndata%ns
   if (ndata%lcheck_sc_nor(is)) then
      isc = isc+1
      ndata%sc_nor_to_s(isc) = is
      ndata%s_to_sc_nor(is) = isc
   end if
end do

! Inter-halo conversions
allocate(ndata%sa_to_sb(ndata%nsa))
allocate(ndata%sa_to_sc(ndata%nsa))
do isa=1,ndata%nsa
   is = ndata%sa_to_s(isa)
   isb = ndata%s_to_sb(is)
   isc = ndata%s_to_sc(is)
   ndata%sa_to_sb(isa) = isb
   ndata%sa_to_sc(isa) = isc
end do
allocate(ndata%sb_to_sc(ndata%nsb))
allocate(ndata%sc_to_sb(ndata%nsc))
allocate(ndata%sb_to_sc_nor(ndata%nsb))
call msi(ndata%sc_to_sb)
do isb=1,ndata%nsb
   is = ndata%sb_to_s(isb)
   isc = ndata%s_to_sc(is)
   ndata%sb_to_sc(isb) = isc
   ndata%sc_to_sb(isc) = isb
   isc = ndata%s_to_sc_nor(is)
   ndata%sb_to_sc_nor(isb) = isc
end do

! Global <-> local conversions for data
h_n_s_max_loc = 0
do il0i=1,geom%nl0i
   h_n_s_max_loc = max(h_n_s_max_loc,ndata%h(il0i)%n_s)
end do
allocate(interph_lg(h_n_s_max_loc,geom%nl0i))
do il0i=1,geom%nl0i
   i_s_loc = 0
   do i_s=1,ndata%hfull(il0i)%n_s
      if (ndata%lcheck_h(i_s,il0i)) then
         i_s_loc = i_s_loc+1
         interph_lg(i_s_loc,il0i) = i_s
      end if
   end do
end do
s_n_s_max_loc = 0
do il1=1,ndata%nl1
   s_n_s_max_loc = max(s_n_s_max_loc,ndata%s(il1)%n_s)
end do
allocate(interps_lg(s_n_s_max_loc,ndata%nl1))
do il1=1,ndata%nl1
   i_s_loc = 0
   do i_s=1,ndata%sfull(il1)%n_s
      if (ndata%lcheck_s(i_s,il1)) then
         i_s_loc = i_s_loc+1
         interps_lg(i_s_loc,il1) = i_s
      end if
   end do
end do

! Local data

! Horizontal interpolation
do il0i=1,geom%nl0i
   ndata%h(il0i)%prefix = 'h'
   ndata%h(il0i)%n_src = ndata%nc1b
   ndata%h(il0i)%n_dst = geom%nc0a
   call linop_alloc(ndata%h(il0i))
   do i_s_loc=1,ndata%h(il0i)%n_s
      i_s = interph_lg(i_s_loc,il0i)
      ndata%h(il0i)%row(i_s_loc) = geom%c0_to_c0a(ndata%hfull(il0i)%row(i_s))
      ndata%h(il0i)%col(i_s_loc) = ndata%c1_to_c1b(ndata%hfull(il0i)%col(i_s))
      ndata%h(il0i)%S(i_s_loc) = ndata%hfull(il0i)%S(i_s)
   end do
   call linop_reorder(ndata%h(il0i))
end do

! Subsampling horizontal interpolation
do il1=1,ndata%nl1
   ndata%s(il1)%prefix = 's'
   ndata%s(il1)%n_src = ndata%nc1b
   ndata%s(il1)%n_dst = ndata%nc1b
   call linop_alloc(ndata%s(il1))
   do i_s_loc=1,ndata%s(il1)%n_s
      i_s = interps_lg(i_s_loc,il1)
      ndata%s(il1)%row(i_s_loc) = ndata%c1_to_c1b(ndata%sfull(il1)%row(i_s))
      ndata%s(il1)%col(i_s_loc) = ndata%c1_to_c1b(ndata%sfull(il1)%col(i_s))
      ndata%s(il1)%S(i_s_loc) = ndata%sfull(il1)%S(i_s)
   end do
   call linop_reorder(ndata%s(il1))
end do

! Convolution
ndata%c%n_src = ndata%nsc
ndata%c%n_dst = ndata%nsc
do i_s=1,ndata%c%n_s
   ndata%c%row(i_s) = ndata%s_to_sc(ndata%c%row(i_s))
   ndata%c%col(i_s) = ndata%s_to_sc(ndata%c%col(i_s))
end do
call linop_reorder(ndata%c)

! Convolution for normalization
ndata%c_nor%n_src = ndata%nsc
ndata%c_nor%n_dst = ndata%nsc
do i_s=1,ndata%c_nor%n_s
   ndata%c_nor%row(i_s) = ndata%s_to_sc_nor(ndata%c_nor%row(i_s))
   ndata%c_nor%col(i_s) = ndata%s_to_sc_nor(ndata%c_nor%col(i_s))
end do
call linop_reorder(ndata%c_nor)

! Conversions
allocate(ndata%sb_to_c1b(ndata%nsb))
allocate(ndata%sb_to_l1(ndata%nsb))
allocate(ndata%c1bl1_to_sb(ndata%nc1b,ndata%nl1))
call msi(ndata%c1bl1_to_sb)
do isb=1,ndata%nsb
   is = ndata%sb_to_s(isb)
   il1 = ndata%s_to_l1(is)
   ic1 = ndata%s_to_c1(is)
   ic1b = ndata%c1_to_c1b(ic1)
   ndata%sb_to_c1b(isb) = ic1b
   ndata%sb_to_l1(isb) = il1
   ndata%c1bl1_to_sb(ic1b,il1) = isb
end do

! Get global distribution of the subgrid on ioproc
if (mpl%main) then
   ! Allocation
   allocate(s_to_sa(ndata%ns))

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         nsa = ndata%nsa
      else
         ! Receive dimension on ioproc
         call mpl_recv(nsa,iproc,mpl%tag)
      end if

      ! Allocation
      allocate(sa_to_s(nsa))

      if (iproc==mpl%ioproc) then
         ! Copy data
         sa_to_s = ndata%sa_to_s
      else
         ! Receive data on ioproc
         call mpl_recv(nsa,sa_to_s,iproc,mpl%tag+1)
      end if

      ! Fill s_to_sa
      do isa=1,nsa
         s_to_sa(sa_to_s(isa)) = isa
      end do

      ! Release memory
      deallocate(sa_to_s)
   end do
else
   ! Send dimensions to ioproc
   call mpl_send(ndata%nsa,mpl%ioproc,mpl%tag)

   ! Send data to ioproc
   call mpl_send(ndata%nsa,ndata%sa_to_s,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nsa = ndata%nsa
         nsb = ndata%nsb
         nsc = ndata%nsc
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nsa,iproc,mpl%tag)
         call mpl_recv(nsb,iproc,mpl%tag+1)
         call mpl_recv(nsc,iproc,mpl%tag+2)
      end if

      ! Allocation
      allocate(sb_to_s(nsb))
      allocate(sc_to_s(nsc))
      allocate(sa_to_sb(nsa))
      allocate(sa_to_sc(nsa))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         sb_to_s = ndata%sb_to_s
         sc_to_s = ndata%sc_to_s
         sa_to_sb = ndata%sa_to_sb
         sa_to_sc = ndata%sa_to_sc
      else
         ! Receive data on ioproc
         call mpl_recv(nsb,sb_to_s,iproc,mpl%tag+3)
         call mpl_recv(nsc,sc_to_s,iproc,mpl%tag+4)
         call mpl_recv(nsa,sa_to_sb,iproc,mpl%tag+5)
         call mpl_recv(nsa,sa_to_sc,iproc,mpl%tag+6)
      end if

      ! Allocation
      comAB(iproc)%nred = nsa
      comAB(iproc)%next = nsb
      allocate(comAB(iproc)%ext_to_proc(comAB(iproc)%next))
      allocate(comAB(iproc)%ext_to_red(comAB(iproc)%next))
      allocate(comAB(iproc)%red_to_ext(comAB(iproc)%nred))
      comAC(iproc)%nred = nsa
      comAC(iproc)%next = nsc
      allocate(comAC(iproc)%ext_to_proc(comAC(iproc)%next))
      allocate(comAC(iproc)%ext_to_red(comAC(iproc)%next))
      allocate(comAC(iproc)%red_to_ext(comAC(iproc)%nred))

      ! AB communication
      do isb=1,nsb
         is = sb_to_s(isb)
         ic1 = ndata%s_to_c1(is)
         ic0 = ndata%c1_to_c0(ic1)
         comAB(iproc)%ext_to_proc(isb) = geom%c0_to_proc(ic0)
         isa = s_to_sa(is)
         comAB(iproc)%ext_to_red(isb) = isa
      end do
      comAB(iproc)%red_to_ext = sa_to_sb

      ! AC communication
      do isc=1,nsc
         is = sc_to_s(isc)
         ic1 = ndata%s_to_c1(is)
         ic0 = ndata%c1_to_c0(ic1)
         comAC(iproc)%ext_to_proc(isc) = geom%c0_to_proc(ic0)
         isa = s_to_sa(is)
         comAC(iproc)%ext_to_red(isc) = isa
      end do
      comAC(iproc)%red_to_ext = sa_to_sc

      ! Release memory
      deallocate(sb_to_s)
      deallocate(sc_to_s)
      deallocate(sa_to_sb)
      deallocate(sa_to_sc)
   end do

   ! Communication setup
   call com_setup(comAB)
   call com_setup(comAC)

   ! Release memory
   deallocate(s_to_sa)
else
   ! Send dimensions to ioproc
   call mpl_send(ndata%nsa,mpl%ioproc,mpl%tag)
   call mpl_send(ndata%nsb,mpl%ioproc,mpl%tag+1)
   call mpl_send(ndata%nsc,mpl%ioproc,mpl%tag+2)

   ! Send data to ioproc
   call mpl_send(ndata%nsb,ndata%sb_to_s,mpl%ioproc,mpl%tag+3)
   call mpl_send(ndata%nsc,ndata%sc_to_s,mpl%ioproc,mpl%tag+4)
   call mpl_send(ndata%nsa,ndata%sa_to_sb,mpl%ioproc,mpl%tag+5)
   call mpl_send(ndata%nsa,ndata%sa_to_sc,mpl%ioproc,mpl%tag+6)
end if
mpl%tag = mpl%tag+7

! Communication broadcast
ndata%AB%prefix = 'AB'
call com_bcast(comAB,ndata%AB)
ndata%AC%prefix = 'AC'
call com_bcast(comAC,ndata%AC)

! Copy mpicom
ndata%mpicom = nam%mpicom

! End associate
end associate

end subroutine compute_mpi_c

end module nicas_parameters_mpi
