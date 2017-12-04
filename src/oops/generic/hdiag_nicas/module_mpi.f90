!----------------------------------------------------------------------
! Module: module_mpi.f90
!> Purpose: compute NICAS parameters MPI distribution
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_mpi

use netcdf
use omp_lib
use tools_const, only: pi,rad2deg,req,sphere_dist
use tools_display, only: msgerror,prog_init,prog_print
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_com, only: comtype,com_dealloc,com_setup,com_bcast
use type_linop, only: linop_alloc,linop_copy,linop_reorder
use type_mpl, only: mpl,mpl_send,mpl_recv
use type_nam, only: namtype
use type_ndata, only: ndatatype,ndataloctype
use type_randgen, only: initialize_sampling

implicit none

private
public :: compute_mpi

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi
!> Purpose: compute NICAS MPI distribution
!----------------------------------------------------------------------
subroutine compute_mpi(ndata,ndataloc)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata       !< NICAS data
type(ndataloctype),intent(inout) :: ndataloc !< NICAS data, local

! Local variables
integer :: iproc,il0i,ic0,ic0a,ic1,ic2,jc2,ic1b,ic2b,il0,il1,isa,isb,isc,i_s,i_s_loc,is,js,s_n_s_max,s_n_s_max_loc
integer :: interph_row_proc(ndata%h(1)%n_s,ndata%geom%nl0i)
integer :: nsa,nsb,nsc
integer,allocatable :: ic1b_to_ic1(:),ic1_to_ic1b(:),ic2il1_to_ic2b(:,:)
integer,allocatable :: interph_i_s_lg(:,:),interps_i_s_lg(:,:),convol_i_s_lg(:)
integer,allocatable :: isa_to_is(:),is_to_isa(:),isb_to_is(:),is_to_isb(:),isc_to_is(:),is_to_isc(:)
integer,allocatable :: isa_to_is_copy(:),isb_to_is_copy(:),isc_to_is_copy(:),isa_to_isb_copy(:),isa_to_isc_copy(:)
logical :: lcheck_nc1b(ndata%nc1),lcheck_nc2b(ndata%nc1,ndata%nl1)
logical :: lcheck_nsa(ndata%ns),lcheck_nsb(ndata%ns),lcheck_nsc(ndata%ns)
logical :: lcheck_h(ndata%h(1)%n_s,ndata%geom%nl0i),lcheck_c(ndata%c%n_s)
logical,allocatable :: lcheck_s(:,:)
type(comtype) :: comAB(mpl%nproc),comAC(mpl%nproc)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Allocation
s_n_s_max = 0
do il1=1,ndata%nl1
   s_n_s_max = max(s_n_s_max,ndata%s(il1)%n_s)
end do
allocate(lcheck_s(s_n_s_max,ndata%nl1))

! Find on which processor are the grid-points and what is their local index for interpolation
do il0i=1,geom%nl0i
   interph_row_proc(1:ndata%h(il0i)%n_s,il0i) = geom%ic0_to_iproc(ndata%h(il0i)%row)
end do

! Copy number of levels
ndataloc%nl1 = ndata%nl1

! Allocation
allocate(ndataloc%nc2b(ndataloc%nl1))
allocate(ndataloc%h(geom%nl0i))
allocate(ndataloc%s(ndataloc%nl1))

! Halo definitions

! Halo A
lcheck_nsa = .false.
do is=1,ndata%ns
   ic1 = ndata%is_to_ic1(is)
   ic0 = ndata%ic1_to_ic0(ic1)
   il1 = ndata%is_to_il1(is)
   il0 = ndata%il1_to_il0(il1)
   if (geom%mask(ic0,il0).and.(geom%ic0_to_iproc(ic0)==mpl%myproc)) lcheck_nsa(is) = .true.
end do
ndataloc%nsa = count(lcheck_nsa)

! Halo B

! Horizontal interpolation
lcheck_h = .false.
lcheck_nc1b = .false.
do il0i=1,geom%nl0i
   do i_s=1,ndata%h(il0i)%n_s
      if (interph_row_proc(i_s,il0i)==mpl%myproc) then
         ic1 = ndata%h(il0i)%col(i_s)
         lcheck_h(i_s,il0i) = .true.
         lcheck_nc1b(ic1) = .true.
      end if
   end do
   ndataloc%h(il0i)%n_s = count(lcheck_h(:,il0i))
end do
ndataloc%nc1b = count(lcheck_nc1b)

! Subsampling horizontal interpolation
lcheck_nc2b = .false.
lcheck_nsb = .false.
lcheck_s = .false.
s_n_s_max_loc = 0
do il1=1,ndataloc%nl1
   do i_s=1,ndata%s(il1)%n_s
      ic1 = ndata%s(il1)%row(i_s)
      if (lcheck_nc1b(ic1)) then
         jc2 = ndata%s(il1)%col(i_s)
         js = ndata%ic2il1_to_is(jc2,il1)
         lcheck_nc2b(jc2,il1) = .true.
         lcheck_nsb(js) = .true.
         lcheck_s(i_s,il1) = .true.
      end if
   end do
   ndataloc%nc2b(il1) = count(lcheck_nc2b(:,il1))
   ndataloc%s(il1)%n_s = count(lcheck_s(:,il1))
   s_n_s_max_loc = max(s_n_s_max_loc,ndataloc%s(il1)%n_s)
end do
ndataloc%nsb = count(lcheck_nsb)

! Halo C
if (nam%mpicom==1) then
   ! 1 communication step
   lcheck_nsc = lcheck_nsb
   lcheck_c = .false.
   do i_s=1,ndata%c%n_s
      is = ndata%c%row(i_s)
      js = ndata%c%col(i_s)
      if (lcheck_nsb(is).or.lcheck_nsb(js)) then
         lcheck_nsc(is) = .true.
         lcheck_nsc(js) = .true.
         lcheck_c(i_s) = .true.
      end if
   end do
elseif (nam%mpicom==2) then
   ! 2 communication steps
   lcheck_nsc = lcheck_nsb
   lcheck_c = .false.
   do i_s=1,ndata%c%n_s
      is = ndata%c%row(i_s)
      js = ndata%c%col(i_s)
      if (lcheck_nsa(is).or.lcheck_nsa(js)) then
         lcheck_nsc(is) = .true.
         lcheck_nsc(js) = .true.
         lcheck_c(i_s) = .true.
      end if
   end do
end if
ndataloc%nsc = count(lcheck_nsc)
ndataloc%c%n_s = count(lcheck_c)

! Check halos consistency
do is=1,ndata%ns
   if (lcheck_nsa(is).and.(.not.lcheck_nsb(is))) then
      call msgerror('point in halo A but not in halo B')
   end if
   if (lcheck_nsa(is).and.(.not.lcheck_nsc(is))) then
      call msgerror('point in halo A but not in halo C')
   end if
   if (lcheck_nsb(is).and.(.not.lcheck_nsc(is))) then
      call msgerror('point in halo B but not in halo C')
   end if
end do

! Global <-> local conversions for fields

! Halo A
if (ndataloc%nsa>0) allocate(isa_to_is(ndataloc%nsa))
isa = 0
do is=1,ndata%ns
   if (lcheck_nsa(is)) then
      isa = isa+1
      if (ndataloc%nsa>0) isa_to_is(isa) = is
   end if
end do

! Halo B
if (ndataloc%nc1b>0) allocate(ic1b_to_ic1(ndataloc%nc1b))
allocate(ic1_to_ic1b(ndata%nc1))
call msi(ic1_to_ic1b)
ic1b = 0
do ic1=1,ndata%nc1
   if (lcheck_nc1b(ic1)) then
      ic1b = ic1b+1
      if (ndataloc%nc1b>0) ic1b_to_ic1(ic1b) = ic1
      ic1_to_ic1b(ic1) = ic1b
   end if
end do

allocate(ic2il1_to_ic2b(ndata%nc1,ndata%nl1))
call msi(ic2il1_to_ic2b)
do il1=1,ndataloc%nl1
   if (ndataloc%nc2b(il1)>0) then
      ic2b = 0
      do ic2=1,ndata%nc2(il1)
         if (lcheck_nc2b(ic2,il1)) then
            ic2b = ic2b+1
            ic2il1_to_ic2b(ic2,il1) = ic2b
         end if
      end do
   end if
end do

if (ndataloc%nsb>0) allocate(isb_to_is(ndataloc%nsb))
allocate(is_to_isb(ndata%ns))
call msi(is_to_isb)
isb = 0
do is=1,ndata%ns
   if (lcheck_nsb(is)) then
      isb = isb+1
      if (ndataloc%nsb>0) isb_to_is(isb) = is
      is_to_isb(is) = isb
   end if
end do

! Halo C
if (ndataloc%nsc>0) allocate(isc_to_is(ndataloc%nsc))
allocate(is_to_isc(ndata%ns))
call msi(is_to_isc)
isc = 0
do is=1,ndata%ns
   if (lcheck_nsc(is)) then
      isc = isc+1
      if (ndataloc%nsc>0) isc_to_is(isc) = is
      is_to_isc(is) = isc
   end if
end do

! Inter-halo conversions
if ((ndataloc%nsa>0).and.(ndataloc%nsb>0).and.(ndataloc%nsc>0)) then
   allocate(ndataloc%isa_to_isb(ndataloc%nsa))
   allocate(ndataloc%isa_to_isc(ndataloc%nsa))
   do isa=1,ndataloc%nsa
      is = isa_to_is(isa)
      isb = is_to_isb(is)
      isc = is_to_isc(is)
      ndataloc%isa_to_isb(isa) = isb
      ndataloc%isa_to_isc(isa) = isc
   end do
   allocate(ndataloc%isb_to_isc(ndataloc%nsb))
   do isb=1,ndataloc%nsb
      is = isb_to_is(isb)
      isc = is_to_isc(is)
      ndataloc%isb_to_isc(isb) = isc
   end do
end if

! Global <-> local conversions for data
allocate(interph_i_s_lg(ndataloc%h(1)%n_s,geom%nl0i))
do il0i=1,geom%nl0i
   i_s_loc = 0
   do i_s=1,ndata%h(il0i)%n_s
      if (lcheck_h(i_s,il0i)) then
         i_s_loc = i_s_loc+1
         interph_i_s_lg(i_s_loc,il0i) = i_s
      end if
   end do
end do
if (s_n_s_max_loc>0) then
   allocate(interps_i_s_lg(s_n_s_max_loc,ndataloc%nl1))
   do il1=1,ndataloc%nl1
      i_s_loc = 0
      do i_s=1,ndata%s(il1)%n_s
         if (lcheck_s(i_s,il1)) then
            i_s_loc = i_s_loc+1
            interps_i_s_lg(i_s_loc,il1) = i_s
         end if
      end do
   end do
end if
if (ndataloc%c%n_s>0) then
   allocate(convol_i_s_lg(ndataloc%c%n_s))
   i_s_loc = 0
   do i_s=1,ndata%c%n_s
      if (lcheck_c(i_s)) then
         i_s_loc = i_s_loc+1
         convol_i_s_lg(i_s_loc) = i_s
      end if
   end do
end if

! Local data

! Horizontal interpolation
do il0i=1,geom%nl0i
   ndataloc%h(il0i)%prefix = 'h'
   ndataloc%h(il0i)%n_src = ndataloc%nc1b
   ndataloc%h(il0i)%n_dst = geom%nc0a
   call linop_alloc(ndataloc%h(il0i))
   do i_s_loc=1,ndataloc%h(il0i)%n_s
      i_s = interph_i_s_lg(i_s_loc,il0i)
      ndataloc%h(il0i)%row(i_s_loc) = geom%ic0_to_ic0a(ndata%h(il0i)%row(i_s))
      ndataloc%h(il0i)%col(i_s_loc) = ic1_to_ic1b(ndata%h(il0i)%col(i_s))
      ndataloc%h(il0i)%S(i_s_loc) = ndata%h(il0i)%S(i_s)
   end do
   call linop_reorder(ndataloc%h(il0i))
end do

! Vertical interpolation
call linop_copy(ndata%v,ndataloc%v)
if (ndataloc%nc1b>0) then
   allocate(ndataloc%vbot(ndataloc%nc1b))
   allocate(ndataloc%vtop(ndataloc%nc1b))
   ndataloc%vbot = ndata%vbot(ic1b_to_ic1)
   ndataloc%vtop = ndata%vtop(ic1b_to_ic1)
end if

! Subsampling horizontal interpolation
do il1=1,ndata%nl1
   ndataloc%s(il1)%prefix = 's'
   ndataloc%s(il1)%n_src = ndataloc%nc2b(il1)
   ndataloc%s(il1)%n_dst = ndataloc%nc1b
   if (ndataloc%s(il1)%n_s>0) then
      call linop_alloc(ndataloc%s(il1))
      do i_s_loc=1,ndataloc%s(il1)%n_s
         i_s = interps_i_s_lg(i_s_loc,il1)
         ndataloc%s(il1)%row(i_s_loc) = ic1_to_ic1b(ndata%s(il1)%row(i_s))
         ndataloc%s(il1)%col(i_s_loc) = ic2il1_to_ic2b(ndata%s(il1)%col(i_s),il1)
         ndataloc%s(il1)%S(i_s_loc) = ndata%s(il1)%S(i_s)
      end do
      call linop_reorder(ndataloc%s(il1))
   end if
end do

! Copy
if (ndataloc%nsb>0) then
   allocate(ndataloc%isb_to_ic2b(ndataloc%nsb))
   allocate(ndataloc%isb_to_il1(ndataloc%nsb))
   call msi(ndataloc%isb_to_ic2b)
   do isb=1,ndataloc%nsb
      is = isb_to_is(isb)
      il1 = ndata%is_to_il1(is)
      ic2 = ndata%is_to_ic2(is)
      ic2b = ic2il1_to_ic2b(ic2,il1)
      ndataloc%isb_to_ic2b(isb) = ic2b
      ndataloc%isb_to_il1(isb) = il1
   end do
end if

! Convolution
if (ndataloc%c%n_s>0) then
   ndataloc%c%prefix = 'c'
   ndataloc%c%n_src = ndataloc%nsc
   ndataloc%c%n_dst = ndataloc%nsc
   call linop_alloc(ndataloc%c)
   do i_s_loc=1,ndataloc%c%n_s
      i_s = convol_i_s_lg(i_s_loc)
      ndataloc%c%row(i_s_loc) = is_to_isc(ndata%c%row(i_s))
      ndataloc%c%col(i_s_loc) = is_to_isc(ndata%c%col(i_s))
      ndataloc%c%S(i_s_loc) = ndata%c%S(i_s)
   end do
   call linop_reorder(ndataloc%c)
end if

! Print local parameters
write(mpl%unit,'(a7,a,i4)') '','Local parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0a =      ',geom%nc0a
write(mpl%unit,'(a10,a,i8)') '','nc1b =      ',ndataloc%nc1b
do il1=1,ndataloc%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2b(',il1,') =  ',ndataloc%nc2b(il1)
end do
write(mpl%unit,'(a10,a,i8)') '','nsa =       ',ndataloc%nsa
write(mpl%unit,'(a10,a,i8)') '','nsb =       ',ndataloc%nsb
write(mpl%unit,'(a10,a,i8)') '','nsc =       ',ndataloc%nsc
do il0i=1,geom%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',ndataloc%h(il0i)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','v%n_s =     ',ndataloc%v%n_s
do il1=1,ndataloc%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',ndataloc%s(il1)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','c%n_s =     ',ndataloc%c%n_s

if (mpl%main) then
   ! Illustration
   allocate(ndata%halo(geom%nc0))
   ndata%halo = 0
   do i_s=1,ndataloc%c%n_s
      ic0 = ndata%ic1_to_ic0(ndata%is_to_ic1(isc_to_is(ndataloc%c%row(i_s))))
      ndata%halo(ic0) = 1
      ic0 = ndata%ic1_to_ic0(ndata%is_to_ic1(isc_to_is(ndataloc%c%col(i_s))))
      ndata%halo(ic0) = 1
   end do
   do i_s=1,ndataloc%h(1)%n_s
      ic0 = ndata%ic1_to_ic0(ic1b_to_ic1(ndataloc%h(1)%col(i_s)))
      ndata%halo(ic0) = 2
   end do
   do isa=1,ndataloc%nsa
      ic0 = ndata%ic1_to_ic0(ndata%is_to_ic1(isa_to_is(isa)))
      ndata%halo(ic0) = 3
   end do
end if

! Release memory
deallocate(ic2il1_to_ic2b)
if (ndataloc%nc1b>0) then
   deallocate(ic1b_to_ic1)
   deallocate(ic1_to_ic1b)
end if
deallocate(interph_i_s_lg)
if (s_n_s_max_loc>0) deallocate(interps_i_s_lg)
if (ndataloc%c%n_s>0) deallocate(convol_i_s_lg)

! Copy norm over processors
allocate(ndataloc%norm(geom%nc0a,geom%nl0))
if (nam%lsqrt) allocate(ndataloc%norm_sqrt(ndataloc%nsb))
do ic0=1,geom%nc0
   if (geom%ic0_to_iproc(ic0)==mpl%myproc) then
      ic0a = geom%ic0_to_ic0a(ic0)
      ndataloc%norm(ic0a,:) = ndata%norm(ic0,:)
   end if
end do
if (nam%lsqrt) then
   do isb=1,ndataloc%nsb
      is = isb_to_is(isb)
      ndataloc%norm_sqrt(isb) = ndata%norm_sqrt(is)
   end do
end if

! Get global distribution of the subgrid on ioproc
if (mpl%main) then
   ! Allocation
   allocate(is_to_isa(ndata%ns))

   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy dimension
         nsa = ndataloc%nsa
      else
         ! Receive dimension on ioproc
         call mpl_recv(nsa,iproc,mpl%tag)
      end if

      ! Allocation
      allocate(isa_to_is_copy(nsa))

      if (iproc==mpl%ioproc) then
         ! Copy data
         isa_to_is_copy = isa_to_is
      else
         ! Receive data on ioproc
         call mpl_recv(nsa,isa_to_is_copy,iproc,mpl%tag+1)
      end if

      ! Fill is_to_isa
      do isa=1,nsa
         is_to_isa(isa_to_is_copy(isa)) = isa
      end do

      ! Release memory
      deallocate(isa_to_is_copy)
   end do
else
   ! Send dimensions to ioproc
   call mpl_send(ndataloc%nsa,mpl%ioproc,mpl%tag)

   ! Send data to ioproc
   call mpl_send(ndataloc%nsa,isa_to_is,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+1

! Setup communications
if (mpl%main) then
   do iproc=1,mpl%nproc
      ! Communicate dimensions
      if (iproc==mpl%ioproc) then
         ! Copy dimensions
         nsa = ndataloc%nsa
         nsb = ndataloc%nsb
         nsc = ndataloc%nsc
      else
         ! Receive dimensions on ioproc
         call mpl_recv(nsa,iproc,mpl%tag)
         call mpl_recv(nsb,iproc,mpl%tag+1)
         call mpl_recv(nsc,iproc,mpl%tag+2)
      end if

      ! Allocation
      allocate(isb_to_is_copy(nsb))
      allocate(isc_to_is_copy(nsc))
      allocate(isa_to_isb_copy(nsa))
      allocate(isa_to_isc_copy(nsa))

      ! Communicate data
      if (iproc==mpl%ioproc) then
         ! Copy data
         isb_to_is_copy = isb_to_is
         isc_to_is_copy = isc_to_is
         isa_to_isb_copy = ndataloc%isa_to_isb
         isa_to_isc_copy = ndataloc%isa_to_isc
      else
         ! Receive data on ioproc
         call mpl_recv(nsb,isb_to_is_copy,iproc,mpl%tag+3)
         call mpl_recv(nsc,isc_to_is_copy,iproc,mpl%tag+4)
         call mpl_recv(nsa,isa_to_isb_copy,iproc,mpl%tag+5)
         call mpl_recv(nsa,isa_to_isc_copy,iproc,mpl%tag+6)
      end if

      ! Allocation
      comAB(iproc)%nred = nsa
      comAB(iproc)%next = nsb
      allocate(comAB(iproc)%iext_to_iproc(comAB(iproc)%next))
      allocate(comAB(iproc)%iext_to_ired(comAB(iproc)%next))
      allocate(comAB(iproc)%ired_to_iext(comAB(iproc)%nred))
      comAC(iproc)%nred = nsa
      comAC(iproc)%next = nsc
      allocate(comAC(iproc)%iext_to_iproc(comAC(iproc)%next))
      allocate(comAC(iproc)%iext_to_ired(comAC(iproc)%next))
      allocate(comAC(iproc)%ired_to_iext(comAC(iproc)%nred))

      ! Initialization
      do isb=1,nsb
         ! Check for points that are in zone B but are not in zone A
         is = isb_to_is_copy(isb)
         ic1 = ndata%is_to_ic1(is)
         ic0 = ndata%ic1_to_ic0(ic1)
         comAB(iproc)%iext_to_iproc(isb) = geom%ic0_to_iproc(ic0)
         isa = is_to_isa(is)
         comAB(iproc)%iext_to_ired(isb) = isa
      end do
      comAB(iproc)%ired_to_iext = isa_to_isb_copy
      do isc=1,nsc
         ! Check for points that are in zone C but are not in zone A
         is = isc_to_is_copy(isc)
         ic1 = ndata%is_to_ic1(is)
         ic0 = ndata%ic1_to_ic0(ic1)
         comAC(iproc)%iext_to_iproc(isc) = geom%ic0_to_iproc(ic0)
         isa = is_to_isa(is)
         comAC(iproc)%iext_to_ired(isc) = isa
      end do
      comAC(iproc)%ired_to_iext = isa_to_isc_copy

      ! Release memory
      deallocate(isb_to_is_copy)
      deallocate(isc_to_is_copy)
      deallocate(isa_to_isb_copy)
      deallocate(isa_to_isc_copy)
   end do

   ! Communication setup
   call com_setup(mpl%nproc,comAB)
   call com_setup(mpl%nproc,comAC)

   ! Release memory
   deallocate(is_to_isa)
else
   ! Send dimensions to ioproc
   call mpl_send(ndataloc%nsa,mpl%ioproc,mpl%tag)
   call mpl_send(ndataloc%nsb,mpl%ioproc,mpl%tag+1)
   call mpl_send(ndataloc%nsc,mpl%ioproc,mpl%tag+2)

   ! Send data to ioproc
   call mpl_send(ndataloc%nsb,isb_to_is,mpl%ioproc,mpl%tag+3)
   call mpl_send(ndataloc%nsc,isc_to_is,mpl%ioproc,mpl%tag+4)
   call mpl_send(ndataloc%nsa,ndataloc%isa_to_isb,mpl%ioproc,mpl%tag+5)
   call mpl_send(ndataloc%nsa,ndataloc%isa_to_isc,mpl%ioproc,mpl%tag+6)
end if
mpl%tag = mpl%tag+8

! Communication broadcast
ndataloc%AB%prefix = 'AB'
call com_bcast(mpl%nproc,comAB,ndataloc%AB)
ndataloc%AC%prefix = 'AC'
call com_bcast(mpl%nproc,comAC,ndataloc%AC)

! End associate
end associate

end subroutine compute_mpi

end module module_mpi
