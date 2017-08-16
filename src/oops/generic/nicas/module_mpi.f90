!----------------------------------------------------------------------
! Module: module_mpi.f90
!> Purpose: compute NICAS parameters
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_mpi

use module_namelist, only: nam
use netcdf
use omp_lib
use tools_const, only: pi,rad2deg,req,sphere_dist
use tools_display, only: msgerror,prog_init,prog_print
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_linop, only: linop_alloc,linop_reorder
use type_mpl, only: mpl
use type_ndata, only: ndatatype,ndataloctype
use type_randgen, only: initialize_sampling

implicit none

! Conversion derived type (private to module_mpi)
type typeconv
   integer,allocatable :: isa_to_is(:)   !< Subgrid, halo A to global
   integer,allocatable :: is_to_isa(:)   !< Subgrid, global to halo A
   integer,allocatable :: isb_to_is(:)   !< Subgrid, halo B to global
   integer,allocatable :: is_to_isb(:)   !< Subgrid, global to halo B
   integer,allocatable :: isc_to_is(:)   !< Subgrid, halo C to global
   integer,allocatable :: is_to_isc(:)   !< Subgrid, global to halo C
end type

private
public :: compute_mpi

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi
!> Purpose: compute NICAS MPI distribution
!----------------------------------------------------------------------
subroutine compute_mpi(ndata,ndataloc_arr)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata                      !< Sampling data
type(ndataloctype),intent(inout) :: ndataloc_arr(nam%nproc) !< Sampling data, local

! Local variables
integer :: il0i,iproc,ic0,ic0a,ic1,ic2,jc2,ic1b,ic2b,il1,isa,isb,isc,i_s,i_s_loc,is,js,jproc,i,s_n_s_max,s_n_s_max_loc
integer :: interph_row_proc(ndata%h(1)%n_s,ndata%nl0i)
integer,allocatable :: ic1b_to_ic1(:),ic1_to_ic1b(:),ic2il1_to_ic2b(:,:)
integer,allocatable :: interph_i_s_lg(:,:),interps_i_s_lg(:,:),convol_i_s_lg(:)
logical :: lcheck_nc1b(ndata%nc1),lcheck_nc2b(ndata%nc1,ndata%nl1)
logical :: lcheck_nsa(ndata%ns),lcheck_nsb(ndata%ns),lcheck_nsc(ndata%ns)
logical :: lcheck_h(ndata%h(1)%n_s,ndata%nl0i),lcheck_c(ndata%c%n_s)
logical,allocatable :: lcheck_s(:,:)
type(typeconv) :: conv(nam%nproc)

! Allocation
do iproc=1,nam%nproc
   ndataloc_arr(iproc)%AB%prefix = 'AB'
   allocate(ndataloc_arr(iproc)%AB%jhalocounts(nam%nproc))
   allocate(ndataloc_arr(iproc)%AB%jexclcounts(nam%nproc))
   allocate(ndataloc_arr(iproc)%AB%jhalodispl(nam%nproc))
   allocate(ndataloc_arr(iproc)%AB%jexcldispl(nam%nproc))
   ndataloc_arr(iproc)%AC%prefix = 'AC'
   allocate(ndataloc_arr(iproc)%AC%jhalocounts(nam%nproc))
   allocate(ndataloc_arr(iproc)%AC%jexclcounts(nam%nproc))
   allocate(ndataloc_arr(iproc)%AC%jhalodispl(nam%nproc))
   allocate(ndataloc_arr(iproc)%AC%jexcldispl(nam%nproc))
end do
s_n_s_max = 0
do il1=1,ndata%nl1
   s_n_s_max = max(s_n_s_max,ndata%s(il1)%n_s)
end do
allocate(lcheck_s(s_n_s_max,ndata%nl1))

! Initialization
ndata%nc0amax = 0
do iproc=1,nam%nproc
   ndataloc_arr(iproc)%AB%jhalocounts = 0
   ndataloc_arr(iproc)%AB%jexclcounts = 0
   ndataloc_arr(iproc)%AB%jhalodispl = 0
   ndataloc_arr(iproc)%AB%jexcldispl = 0
   ndataloc_arr(iproc)%AC%jhalocounts = 0
   ndataloc_arr(iproc)%AC%jexclcounts = 0
   ndataloc_arr(iproc)%AC%jhalodispl = 0
   ndataloc_arr(iproc)%AC%jexcldispl = 0
end do

! Find on which processor are the grid-points and what is their local index for interpolation
do il0i=1,ndata%nl0i
   interph_row_proc(1:ndata%h(il0i)%n_s,il0i) = ndata%ic0_to_iproc(ndata%h(il0i)%row)
end do

do iproc=1,nam%nproc
   ! Copy number of levels
   ndataloc_arr(iproc)%nl0 = ndata%nl0
   ndataloc_arr(iproc)%nl1 = ndata%nl1
   ndataloc_arr(iproc)%nl0i = ndata%nl0i

   ! Allocation
   allocate(ndataloc_arr(iproc)%nc2b(ndataloc_arr(iproc)%nl1))
   allocate(ndataloc_arr(iproc)%h(ndataloc_arr(iproc)%nl0i))
   allocate(ndataloc_arr(iproc)%v(ndataloc_arr(iproc)%nl0i))
   allocate(ndataloc_arr(iproc)%s(ndataloc_arr(iproc)%nl1))

   ! Halo definitions

   ! Halo A
   lcheck_nsa = .false.
   do is=1,ndata%ns
      ic1 = ndata%is_to_ic1(is)
      ic0 = ndata%ic1_to_ic0(ic1)
      if (ndata%ic0_to_iproc(ic0)==iproc) lcheck_nsa(is) = .true.
   end do
   ndataloc_arr(iproc)%nsa = count(lcheck_nsa)

   ! Halo B

   ! Horizontal interpolation
   lcheck_h = .false.
   lcheck_nc1b = .false.
   do il0i=1,ndata%nl0i
      do i_s=1,ndata%h(il0i)%n_s
         if (interph_row_proc(i_s,il0i)==iproc) then
            ic1 = ndata%h(il0i)%col(i_s)
            lcheck_h(i_s,il0i) = .true.
            lcheck_nc1b(ic1) = .true.
         end if
      end do
      ndataloc_arr(iproc)%h(il0i)%n_s = count(lcheck_h(:,il0i))
   end do
   ndataloc_arr(iproc)%nc1b = count(lcheck_nc1b)

   ! Subsampling horizontal interpolation
   lcheck_nc2b = .false.
   lcheck_nsb = .false.
   lcheck_s = .false.
   s_n_s_max_loc = 0
   do il1=1,ndataloc_arr(iproc)%nl1
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
      ndataloc_arr(iproc)%nc2b(il1) = count(lcheck_nc2b(:,il1))
      ndataloc_arr(iproc)%s(il1)%n_s = count(lcheck_s(:,il1))
      s_n_s_max_loc = max(s_n_s_max_loc,ndataloc_arr(iproc)%s(il1)%n_s)
   end do
   ndataloc_arr(iproc)%nsb = count(lcheck_nsb)

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
      lcheck_nsc = lcheck_nsa
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
   ndataloc_arr(iproc)%nsc = count(lcheck_nsc)
   ndataloc_arr(iproc)%c%n_s = count(lcheck_c)

   ! Global <-> local conversions for fields

   ! Halo A
   if (ndataloc_arr(iproc)%nsa>0) allocate(conv(iproc)%isa_to_is(ndataloc_arr(iproc)%nsa))
   allocate(conv(iproc)%is_to_isa(ndata%ns))
   call msi(conv(iproc)%is_to_isa)
   isa = 0
   do is=1,ndata%ns
      if (lcheck_nsa(is)) then
         isa = isa+1
         if (ndataloc_arr(iproc)%nsa>0) conv(iproc)%isa_to_is(isa) = is
         conv(iproc)%is_to_isa(is) = isa
      end if
   end do

   ! Halo B
   if (ndataloc_arr(iproc)%nc1b>0) allocate(ic1b_to_ic1(ndataloc_arr(iproc)%nc1b))
   allocate(ic1_to_ic1b(ndata%nc1))
   call msi(ic1_to_ic1b)
   ic1b = 0
   do ic1=1,ndata%nc1
      if (lcheck_nc1b(ic1)) then
         ic1b = ic1b+1
         if (ndataloc_arr(iproc)%nc1b>0) ic1b_to_ic1(ic1b) = ic1
         ic1_to_ic1b(ic1) = ic1b
      end if
   end do

   allocate(ic2il1_to_ic2b(ndata%nc1,ndata%nl1))
   call msi(ic2il1_to_ic2b)
   do il1=1,ndataloc_arr(iproc)%nl1
      if (ndataloc_arr(iproc)%nc2b(il1)>0) then
         ic2b = 0
         do ic2=1,ndata%nc2(il1)
            if (lcheck_nc2b(ic2,il1)) then
               ic2b = ic2b+1
               ic2il1_to_ic2b(ic2,il1) = ic2b
            end if
         end do
      end if
   end do

   if (ndataloc_arr(iproc)%nsb>0) allocate(conv(iproc)%isb_to_is(ndataloc_arr(iproc)%nsb))
   allocate(conv(iproc)%is_to_isb(ndata%ns))
   call msi(conv(iproc)%is_to_isb)
   isb = 0
   do is=1,ndata%ns
      if (lcheck_nsb(is)) then
         isb = isb+1
         if (ndataloc_arr(iproc)%nsb>0) conv(iproc)%isb_to_is(isb) = is
         conv(iproc)%is_to_isb(is) = isb
      end if
   end do

   ! Halo C
   if (ndataloc_arr(iproc)%nsc>0) allocate(conv(iproc)%isc_to_is(ndataloc_arr(iproc)%nsc))
   allocate(conv(iproc)%is_to_isc(ndata%ns))
   call msi(conv(iproc)%is_to_isc)
   isc = 0
   do is=1,ndata%ns
      if (lcheck_nsc(is)) then
         isc = isc+1
         if (ndataloc_arr(iproc)%nsc>0) conv(iproc)%isc_to_is(isc) = is
         conv(iproc)%is_to_isc(is) = isc
      end if
   end do

   ! Inter-halo conversions
   if ((ndataloc_arr(iproc)%nsa>0).and.(ndataloc_arr(iproc)%nsb>0).and.(ndataloc_arr(iproc)%nsc>0)) then
      allocate(ndataloc_arr(iproc)%isa_to_isb(ndataloc_arr(iproc)%nsa))
      allocate(ndataloc_arr(iproc)%isa_to_isc(ndataloc_arr(iproc)%nsa))
      do isa=1,ndataloc_arr(iproc)%nsa
         is = conv(iproc)%isa_to_is(isa)
         isb = conv(iproc)%is_to_isb(is)
         isc = conv(iproc)%is_to_isc(is)
         ndataloc_arr(iproc)%isa_to_isb(isa) = isb
         ndataloc_arr(iproc)%isa_to_isc(isa) = isc
      end do
      allocate(ndataloc_arr(iproc)%isb_to_isc(ndataloc_arr(iproc)%nsb))
      do isb=1,ndataloc_arr(iproc)%nsb
         is = conv(iproc)%isb_to_is(isb)
         isc = conv(iproc)%is_to_isc(is)
         ndataloc_arr(iproc)%isb_to_isc(isb) = isc
      end do
   end if

   ! Global <-> local conversions for data
   allocate(interph_i_s_lg(ndataloc_arr(iproc)%h(1)%n_s,ndataloc_arr(iproc)%nl0i))
   do il0i=1,ndata%nl0i
      i_s_loc = 0
      do i_s=1,ndata%h(il0i)%n_s
         if (lcheck_h(i_s,il0i)) then
            i_s_loc = i_s_loc+1
            interph_i_s_lg(i_s_loc,il0i) = i_s
         end if
      end do
   end do
   if (s_n_s_max_loc>0) then
      allocate(interps_i_s_lg(s_n_s_max_loc,ndataloc_arr(iproc)%nl1))
      do il1=1,ndataloc_arr(iproc)%nl1
         i_s_loc = 0
         do i_s=1,ndata%s(il1)%n_s
            if (lcheck_s(i_s,il1)) then
               i_s_loc = i_s_loc+1
               interps_i_s_lg(i_s_loc,il1) = i_s
            end if
         end do
      end do
   end if
   if (ndataloc_arr(iproc)%c%n_s>0) then
      allocate(convol_i_s_lg(ndataloc_arr(iproc)%c%n_s))
      i_s_loc = 0
      do i_s=1,ndata%c%n_s
         if (lcheck_c(i_s)) then
            i_s_loc = i_s_loc+1
            convol_i_s_lg(i_s_loc) = i_s
         end if
      end do
   end if

   ! Number of cells
   ndataloc_arr(iproc)%nc0a = count(ndata%ic0_to_iproc==iproc)
   ndata%nc0amax = max(ndata%nc0amax,ndataloc_arr(iproc)%nc0a)

   ! Local data

   ! Horizontal interpolation
   do il0i=1,ndataloc_arr(iproc)%nl0i
      ndataloc_arr(iproc)%h(il0i)%prefix = 'h'
      ndataloc_arr(iproc)%h(il0i)%n_src = ndataloc_arr(iproc)%nc1b
      ndataloc_arr(iproc)%h(il0i)%n_dst = ndataloc_arr(iproc)%nc0a
      call linop_alloc(ndataloc_arr(iproc)%h(il0i))
      do i_s_loc=1,ndataloc_arr(iproc)%h(il0i)%n_s
         i_s = interph_i_s_lg(i_s_loc,il0i)
         ndataloc_arr(iproc)%h(il0i)%row(i_s_loc) = ndata%ic0_to_ic0a(ndata%h(il0i)%row(i_s))
         ndataloc_arr(iproc)%h(il0i)%col(i_s_loc) = ic1_to_ic1b(ndata%h(il0i)%col(i_s))
         ndataloc_arr(iproc)%h(il0i)%S(i_s_loc) = ndata%h(il0i)%S(i_s)
      end do
      call linop_reorder(ndataloc_arr(iproc)%h(il0i))
   end do

   ! Vertical interpolation
   do il0i=1,ndataloc_arr(iproc)%nl0i
      ndataloc_arr(iproc)%v(il0i) = ndata%v(il0i)
   end do
   if (ndataloc_arr(iproc)%nc1b>0) then
      allocate(ndataloc_arr(iproc)%vbot(ndataloc_arr(iproc)%nc1b))
      ndataloc_arr(iproc)%vbot = ndata%vbot(ic1b_to_ic1)
   end if

   ! Subsampling horizontal interpolation
   do il1=1,ndata%nl1
      ndataloc_arr(iproc)%s(il1)%prefix = 's'
      ndataloc_arr(iproc)%s(il1)%n_src = ndataloc_arr(iproc)%nc2b(il1)
      ndataloc_arr(iproc)%s(il1)%n_dst = ndataloc_arr(iproc)%nc1b
      if (ndataloc_arr(iproc)%s(il1)%n_s>0) then
         call linop_alloc(ndataloc_arr(iproc)%s(il1))
         do i_s_loc=1,ndataloc_arr(iproc)%s(il1)%n_s
            i_s = interps_i_s_lg(i_s_loc,il1)
            ndataloc_arr(iproc)%s(il1)%row(i_s_loc) = ic1_to_ic1b(ndata%s(il1)%row(i_s))
            ndataloc_arr(iproc)%s(il1)%col(i_s_loc) = ic2il1_to_ic2b(ndata%s(il1)%col(i_s),il1)
            ndataloc_arr(iproc)%s(il1)%S(i_s_loc) = ndata%s(il1)%S(i_s)
         end do
         call linop_reorder(ndataloc_arr(iproc)%s(il1))
      end if
   end do

   ! Copy
   if (ndataloc_arr(iproc)%nsb>0) then
      allocate(ndataloc_arr(iproc)%isb_to_ic2b(ndataloc_arr(iproc)%nsb))
      allocate(ndataloc_arr(iproc)%isb_to_il1(ndataloc_arr(iproc)%nsb))
      call msi(ndataloc_arr(iproc)%isb_to_ic2b)
      do isb=1,ndataloc_arr(iproc)%nsb
         is = conv(iproc)%isb_to_is(isb)
         il1 = ndata%is_to_il1(is)
         ic2 = ndata%is_to_ic2(is)
         ic2b = ic2il1_to_ic2b(ic2,il1)
         ndataloc_arr(iproc)%isb_to_ic2b(isb) = ic2b
         ndataloc_arr(iproc)%isb_to_il1(isb) = il1
      end do
   end if

   ! Convolution
   if (ndataloc_arr(iproc)%c%n_s>0) then
      ndataloc_arr(iproc)%c%prefix = 'c'
      ndataloc_arr(iproc)%c%n_src = ndataloc_arr(iproc)%nsc
      ndataloc_arr(iproc)%c%n_dst = ndataloc_arr(iproc)%nsc
      call linop_alloc(ndataloc_arr(iproc)%c)
      do i_s_loc=1,ndataloc_arr(iproc)%c%n_s
         i_s = convol_i_s_lg(i_s_loc)
         ndataloc_arr(iproc)%c%row(i_s_loc) = conv(iproc)%is_to_isc(ndata%c%row(i_s))
         ndataloc_arr(iproc)%c%col(i_s_loc) = conv(iproc)%is_to_isc(ndata%c%col(i_s))
         ndataloc_arr(iproc)%c%S(i_s_loc) = ndata%c%S(i_s)
      end do
      call linop_reorder(ndataloc_arr(iproc)%c)
   end if

   ! Print local parameters
   write(mpl%unit,'(a7,a,i4)') '','Local parameters for processor #',iproc
   write(mpl%unit,'(a10,a,i8)') '','nc0a =      ',ndataloc_arr(iproc)%nc0a
   write(mpl%unit,'(a10,a,i8)') '','nc1b =      ',ndataloc_arr(iproc)%nc1b
   do il1=1,ndataloc_arr(iproc)%nl1
      write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2b(',il1,') =  ',ndataloc_arr(iproc)%nc2b(il1)
   end do
   write(mpl%unit,'(a10,a,i8)') '','nsa =       ',ndataloc_arr(iproc)%nsa
   write(mpl%unit,'(a10,a,i8)') '','nsb =       ',ndataloc_arr(iproc)%nsb
   write(mpl%unit,'(a10,a,i8)') '','nsc =       ',ndataloc_arr(iproc)%nsc
   do il0i=1,ndataloc_arr(iproc)%nl0i
      write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',ndataloc_arr(iproc)%h(il0i)%n_s
   end do
   do il0i=1,ndataloc_arr(iproc)%nl0i
      write(mpl%unit,'(a10,a,i3,a,i8)') '','v(',il0i,')%n_s = ',ndataloc_arr(iproc)%v(il0i)%n_s
   end do
   do il1=1,ndataloc_arr(iproc)%nl1
      write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',ndataloc_arr(iproc)%s(il1)%n_s
   end do
   write(mpl%unit,'(a10,a,i8)') '','c%n_s =     ',ndataloc_arr(iproc)%c%n_s

   if (iproc==nam%nproc/2) then
      ! Illustration
      allocate(ndata%halo(ndata%nc0))
      ndata%halo = 0
      do i_s=1,ndataloc_arr(iproc)%c%n_s
         ic0 = ndata%ic1_to_ic0(ndata%is_to_ic1(conv(iproc)%isc_to_is(ndataloc_arr(iproc)%c%row(i_s))))
         ndata%halo(ic0) = 1
         ic0 = ndata%ic1_to_ic0(ndata%is_to_ic1(conv(iproc)%isc_to_is(ndataloc_arr(iproc)%c%col(i_s))))
         ndata%halo(ic0) = 1
      end do
      do i_s=1,ndataloc_arr(iproc)%h(1)%n_s
         ic0 = ndata%ic1_to_ic0(ic1b_to_ic1(ndataloc_arr(iproc)%h(1)%col(i_s)))
         ndata%halo(ic0) = 2
      end do
      do isa=1,ndataloc_arr(iproc)%nsa
         ic0 = ndata%ic1_to_ic0(ndata%is_to_ic1(conv(iproc)%isa_to_is(isa)))
         ndata%halo(ic0) = 3
      end do
   end if

   ! Release memory
   deallocate(ic2il1_to_ic2b)
   if (ndataloc_arr(iproc)%nc1b>0) then
      deallocate(ic1b_to_ic1)
      deallocate(ic1_to_ic1b)
   end if
   deallocate(interph_i_s_lg)
   if (s_n_s_max_loc>0) deallocate(interps_i_s_lg)
   if (ndataloc_arr(iproc)%c%n_s>0) deallocate(convol_i_s_lg)
end do

! Copy norm over processors
do iproc=1,nam%nproc
   allocate(ndataloc_arr(iproc)%norm(ndataloc_arr(iproc)%nc0a,ndataloc_arr(iproc)%nl0))
   if (nam%lsqrt) allocate(ndataloc_arr(iproc)%norm_sqrt(ndataloc_arr(iproc)%nsb))
end do
do ic0=1,ndata%nc0
   iproc = ndata%ic0_to_iproc(ic0)
   ic0a = ndata%ic0_to_ic0a(ic0)
   ndataloc_arr(iproc)%norm(ic0a,1:ndataloc_arr(iproc)%nl0) = ndata%norm(ic0,1:ndata%nl0)
end do
if (nam%lsqrt) then
   do iproc=1,nam%nproc
      do isb=1,ndataloc_arr(iproc)%nsb
         is = conv(iproc)%isb_to_is(isb)
         ndataloc_arr(iproc)%norm_sqrt(isb) = ndata%norm_sqrt(is)
      end do
   end do
end if

! Communications setup
do iproc=1,nam%nproc
   do isb=1,ndataloc_arr(iproc)%nsb
      ! Check for points that are in zone B but are not in zone A
      is = conv(iproc)%isb_to_is(isb)
      ic1 = ndata%is_to_ic1(is)
      ic0 = ndata%ic1_to_ic0(ic1)
      jproc = ndata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         ndataloc_arr(iproc)%AB%jhalocounts(jproc) = ndataloc_arr(iproc)%AB%jhalocounts(jproc)+1

         ! Count of points received on JPROC from IPROC
         ndataloc_arr(jproc)%AB%jexclcounts(iproc) = ndataloc_arr(jproc)%AB%jexclcounts(iproc)+1
      end if
   end do
   do isc=1,ndataloc_arr(iproc)%nsc
      ! Check for points that are in zone C but are not in zone A
      is = conv(iproc)%isc_to_is(isc)
      ic1 = ndata%is_to_ic1(is)
      ic0 = ndata%ic1_to_ic0(ic1)
      jproc = ndata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         ndataloc_arr(iproc)%AC%jhalocounts(jproc) = ndataloc_arr(iproc)%AC%jhalocounts(jproc)+1

         ! Count of points received on JPROC from IPROC
         ndataloc_arr(jproc)%AC%jexclcounts(iproc) = ndataloc_arr(jproc)%AC%jexclcounts(iproc)+1
      end if
   end do
end do

! Compute displacement
do iproc=1,nam%nproc
   ndataloc_arr(iproc)%AB%jhalodispl(1) = 0
   ndataloc_arr(iproc)%AB%jexcldispl(1) = 0
   ndataloc_arr(iproc)%AC%jhalodispl(1) = 0
   ndataloc_arr(iproc)%AC%jexcldispl(1) = 0
   do jproc=2,nam%nproc
      ndataloc_arr(iproc)%AB%jhalodispl(jproc) = & 
    & ndataloc_arr(iproc)%AB%jhalodispl(jproc-1)+ndataloc_arr(iproc)%AB%jhalocounts(jproc-1)
      ndataloc_arr(iproc)%AB%jexcldispl(jproc) = &
    & ndataloc_arr(iproc)%AB%jexcldispl(jproc-1)+ndataloc_arr(iproc)%AB%jexclcounts(jproc-1)
      ndataloc_arr(iproc)%AC%jhalodispl(jproc) = &
    & ndataloc_arr(iproc)%AC%jhalodispl(jproc-1)+ndataloc_arr(iproc)%AC%jhalocounts(jproc-1)
      ndataloc_arr(iproc)%AC%jexcldispl(jproc) = &
    & ndataloc_arr(iproc)%AC%jexcldispl(jproc-1)+ndataloc_arr(iproc)%AC%jexclcounts(jproc-1)
   end do
end do

! Allocation
do iproc=1,nam%nproc
   ndataloc_arr(iproc)%AB%nhalo = sum(ndataloc_arr(iproc)%AB%jhalocounts)
   ndataloc_arr(iproc)%AB%nexcl = sum(ndataloc_arr(iproc)%AB%jexclcounts)
   ndataloc_arr(iproc)%AC%nhalo = sum(ndataloc_arr(iproc)%AC%jhalocounts)
   ndataloc_arr(iproc)%AC%nexcl = sum(ndataloc_arr(iproc)%AC%jexclcounts)
   allocate(ndataloc_arr(iproc)%AB%halo(ndataloc_arr(iproc)%AB%nhalo))
   allocate(ndataloc_arr(iproc)%AB%excl(ndataloc_arr(iproc)%AB%nexcl))
   allocate(ndataloc_arr(iproc)%AC%halo(ndataloc_arr(iproc)%AC%nhalo))
   allocate(ndataloc_arr(iproc)%AC%excl(ndataloc_arr(iproc)%AC%nexcl))
end do

! Fill halo array
do iproc=1,nam%nproc
   ndataloc_arr(iproc)%AB%jhalocounts = 0
   ndataloc_arr(iproc)%AC%jhalocounts = 0
end do
do iproc=1,nam%nproc
   do isb=1,ndataloc_arr(iproc)%nsb
      ! Check for points that ARE in zone B but ARE NOT in zone A
      is = conv(iproc)%isb_to_is(isb)
      ic1 = ndata%is_to_ic1(is)
      ic0 = ndata%ic1_to_ic0(ic1)
      jproc = ndata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         ndataloc_arr(iproc)%AB%jhalocounts(jproc) = ndataloc_arr(iproc)%AB%jhalocounts(jproc)+1

         ! Global index of points sent from IPROC to JPROC
         ndataloc_arr(iproc)%AB%halo(ndataloc_arr(iproc)%AB%jhalodispl(jproc)+ndataloc_arr(iproc)%AB%jhalocounts(jproc)) = is
      end if
   end do
   do isc=1,ndataloc_arr(iproc)%nsc
      ! Check for points that ARE in zone C but ARE NOT in zone A
      is = conv(iproc)%isc_to_is(isc)
      ic1 = ndata%is_to_ic1(is)
      ic0 = ndata%ic1_to_ic0(ic1)
      jproc = ndata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         ndataloc_arr(iproc)%AC%jhalocounts(jproc) = ndataloc_arr(iproc)%AC%jhalocounts(jproc)+1

         ! Global index of points sent from IPROC to JPROC
         ndataloc_arr(iproc)%AC%halo(ndataloc_arr(iproc)%AC%jhalodispl(jproc)+ndataloc_arr(iproc)%AC%jhalocounts(jproc)) = is
      end if
   end do
end do

! Fill excl array
do jproc=1,nam%nproc
   ! Loop over processors sending data to JPROC
   do iproc=1,nam%nproc
      do i=1,ndataloc_arr(iproc)%AB%jhalocounts(jproc)
         ! Global index of points received on JPROC from IPROC
         ndataloc_arr(jproc)%AB%excl(ndataloc_arr(jproc)%AB%jexcldispl(iproc)+i) = &
       & ndataloc_arr(iproc)%AB%halo(ndataloc_arr(iproc)%AB%jhalodispl(jproc)+i)
      end do
      do i=1,ndataloc_arr(iproc)%AC%jhalocounts(jproc)
         ! Global index of points received on JPROC from IPROC
         ndataloc_arr(jproc)%AC%excl(ndataloc_arr(jproc)%AC%jexcldispl(iproc)+i) = &
       & ndataloc_arr(iproc)%AC%halo(ndataloc_arr(iproc)%AC%jhalodispl(jproc)+i)
      end do
   end do
end do

! Transform to local indices
do iproc=1,nam%nproc
   if ((ndataloc_arr(iproc)%nsa>0).and.(ndataloc_arr(iproc)%nsb>0).and.(ndataloc_arr(iproc)%nsc>0)) then
      ndataloc_arr(iproc)%AB%halo = conv(iproc)%is_to_isb(ndataloc_arr(iproc)%AB%halo)
      ndataloc_arr(iproc)%AB%excl = conv(iproc)%is_to_isa(ndataloc_arr(iproc)%AB%excl)
      ndataloc_arr(iproc)%AC%halo = conv(iproc)%is_to_isc(ndataloc_arr(iproc)%AC%halo)
      ndataloc_arr(iproc)%AC%excl = conv(iproc)%is_to_isa(ndataloc_arr(iproc)%AC%excl)
   end if
end do

end subroutine compute_mpi

end module module_mpi
