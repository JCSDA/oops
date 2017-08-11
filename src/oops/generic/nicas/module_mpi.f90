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

use module_apply_com, only: fld_com_gl
use netcdf
use omp_lib
use tools_const, only: pi,rad2deg,req,sphere_dist
use tools_display, only: msgerror,prog_init,prog_print
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr
use tools_nc, only: ncfloat,ncerr
use type_fields, only: fldtype
use type_linop, only: linop_alloc,linop_reorder
use type_mpl, only: mpl
use type_randgen, only: initialize_sampling
use type_sdata, only: sdatatype
implicit none

private
public :: compute_mpi

contains

!----------------------------------------------------------------------
! Subroutine: compute_mpi
!> Purpose: compute NICAS MPI distribution
!----------------------------------------------------------------------
subroutine compute_mpi(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0i,iproc,ic0,ic1,ic2,jc2,ic1b,ic2b,il1,isa,isb,isc,i_s,i_s_loc,is,js,jproc,i,s_n_s_max,s_n_s_max_loc
integer :: interph_row_proc(sdata%h(1)%n_s,sdata%nl0i)
integer,allocatable :: ic1b_to_ic1(:),ic1_to_ic1b(:),ic2il1_to_ic2b(:,:)
integer,allocatable :: interph_i_s_lg(:,:),interps_i_s_lg(:,:),convol_i_s_lg(:)
logical :: lcheck_nc1b(sdata%nc1),lcheck_nc2b(sdata%nc1,sdata%nl1)
logical :: lcheck_nsa(sdata%ns),lcheck_nsb(sdata%ns),lcheck_nsc(sdata%ns)
logical :: lcheck_h(sdata%h(1)%n_s,sdata%nl0i),lcheck_c(sdata%c%n_s)
logical,allocatable :: lcheck_s(:,:)
type(fldtype) :: norm(sdata%nproc)

! Allocation
allocate(sdata%mpi(sdata%nproc))
do iproc=1,sdata%nproc
   sdata%mpi(iproc)%AB%prefix = 'AB'
   sdata%mpi(iproc)%AB%nproc = sdata%nproc
   allocate(sdata%mpi(iproc)%AB%jhalocounts(sdata%nproc))
   allocate(sdata%mpi(iproc)%AB%jexclcounts(sdata%nproc))
   allocate(sdata%mpi(iproc)%AB%jhalodispl(sdata%nproc))
   allocate(sdata%mpi(iproc)%AB%jexcldispl(sdata%nproc))
   sdata%mpi(iproc)%AC%prefix = 'AC'
   sdata%mpi(iproc)%AC%nproc = sdata%nproc
   allocate(sdata%mpi(iproc)%AC%jhalocounts(sdata%nproc))
   allocate(sdata%mpi(iproc)%AC%jexclcounts(sdata%nproc))
   allocate(sdata%mpi(iproc)%AC%jhalodispl(sdata%nproc))
   allocate(sdata%mpi(iproc)%AC%jexcldispl(sdata%nproc))
end do
s_n_s_max = 0
do il1=1,sdata%nl1
   s_n_s_max = max(s_n_s_max,sdata%s(il1)%n_s)
end do
allocate(lcheck_s(s_n_s_max,sdata%nl1))

! Initialization
do iproc=1,sdata%nproc
   sdata%mpi(iproc)%AB%jhalocounts = 0
   sdata%mpi(iproc)%AB%jexclcounts = 0
   sdata%mpi(iproc)%AB%jhalodispl = 0
   sdata%mpi(iproc)%AB%jexcldispl = 0
   sdata%mpi(iproc)%AC%jhalocounts = 0
   sdata%mpi(iproc)%AC%jexclcounts = 0
   sdata%mpi(iproc)%AC%jhalodispl = 0
   sdata%mpi(iproc)%AC%jexcldispl = 0
end do

! Find on which processor are the grid-points and what is their local index for interpolation
do il0i=1,sdata%nl0i
   interph_row_proc(1:sdata%h(il0i)%n_s,il0i) = sdata%ic0_to_iproc(sdata%h(il0i)%row)
end do

do iproc=1,sdata%nproc
   ! Copy number of levels and levels conversion
   sdata%mpi(iproc)%nl0 = sdata%nl0
   sdata%mpi(iproc)%nl1 = sdata%nl1
   sdata%mpi(iproc)%nl0i = sdata%nl0i

   ! Allocation
   allocate(sdata%mpi(iproc)%nc2b(sdata%mpi(iproc)%nl1))
   allocate(sdata%mpi(iproc)%h(sdata%mpi(iproc)%nl0i))
   allocate(sdata%mpi(iproc)%v(sdata%mpi(iproc)%nl0i))
   allocate(sdata%mpi(iproc)%s(sdata%mpi(iproc)%nl1))

   ! Halo definitions

   ! Halo A
   lcheck_nsa = .false.
   do is=1,sdata%ns
      ic1 = sdata%is_to_ic1(is)
      ic0 = sdata%ic1_to_ic0(ic1)
      if (sdata%ic0_to_iproc(ic0)==iproc) lcheck_nsa(is) = .true.
   end do
   sdata%mpi(iproc)%nsa = count(lcheck_nsa)

   ! Halo B

   ! Horizontal interpolation
   lcheck_h = .false.
   lcheck_nc1b = .false.
   do il0i=1,sdata%nl0i
      do i_s=1,sdata%h(il0i)%n_s
         if (interph_row_proc(i_s,il0i)==iproc) then
            ic1 = sdata%h(il0i)%col(i_s)
            lcheck_h(i_s,il0i) = .true.
            lcheck_nc1b(ic1) = .true.
         end if
      end do
      sdata%mpi(iproc)%h(il0i)%n_s = count(lcheck_h(:,il0i))
   end do
   sdata%mpi(iproc)%nc1b = count(lcheck_nc1b)

   ! Subsampling horizontal interpolation
   lcheck_nc2b = .false.
   lcheck_nsb = .false.
   lcheck_s = .false.
   s_n_s_max_loc = 0
   do il1=1,sdata%mpi(iproc)%nl1
      do i_s=1,sdata%s(il1)%n_s
         ic1 = sdata%s(il1)%row(i_s)
         if (lcheck_nc1b(ic1)) then
            jc2 = sdata%s(il1)%col(i_s)
            js = sdata%ic2il1_to_is(jc2,il1)
            lcheck_nc2b(jc2,il1) = .true.
            lcheck_nsb(js) = .true.
            lcheck_s(i_s,il1) = .true.
         end if
      end do
      sdata%mpi(iproc)%nc2b(il1) = count(lcheck_nc2b(:,il1))
      sdata%mpi(iproc)%s(il1)%n_s = count(lcheck_s(:,il1))
      s_n_s_max_loc = max(s_n_s_max_loc,sdata%mpi(iproc)%s(il1)%n_s)
   end do
   sdata%mpi(iproc)%nsb = count(lcheck_nsb)

   ! Halo C
   if (sdata%mpicom==1) then
      ! 1 communication step
      lcheck_nsc = lcheck_nsb
      lcheck_c = .false.
      do i_s=1,sdata%c%n_s
         is = sdata%c%row(i_s)
         js = sdata%c%col(i_s)
         if (lcheck_nsb(is).or.lcheck_nsb(js)) then
            lcheck_nsc(is) = .true.
            lcheck_nsc(js) = .true.
            lcheck_c(i_s) = .true.
         end if
      end do
   elseif (sdata%mpicom==2) then
      ! 2 communication steps
      lcheck_nsc = lcheck_nsa
      lcheck_c = .false.
      do i_s=1,sdata%c%n_s
         is = sdata%c%row(i_s)
         js = sdata%c%col(i_s)
         if (lcheck_nsa(is).or.lcheck_nsa(js)) then
            lcheck_nsc(is) = .true.
            lcheck_nsc(js) = .true.
            lcheck_c(i_s) = .true.
         end if
      end do
   end if
   sdata%mpi(iproc)%nsc = count(lcheck_nsc)
   sdata%mpi(iproc)%c%n_s = count(lcheck_c)

   ! Global <-> local conversions for fields

   ! Halo A
   if (sdata%mpi(iproc)%nsa>0) allocate(sdata%mpi(iproc)%isa_to_is(sdata%mpi(iproc)%nsa))
   allocate(sdata%mpi(iproc)%is_to_isa(sdata%ns))
   call msi(sdata%mpi(iproc)%is_to_isa)
   isa = 0
   do is=1,sdata%ns
      if (lcheck_nsa(is)) then
         isa = isa+1
         if (sdata%mpi(iproc)%nsa>0) sdata%mpi(iproc)%isa_to_is(isa) = is
         sdata%mpi(iproc)%is_to_isa(is) = isa
      end if
   end do

   ! Halo B
   if (sdata%mpi(iproc)%nc1b>0) allocate(ic1b_to_ic1(sdata%mpi(iproc)%nc1b))
   allocate(ic1_to_ic1b(sdata%nc1))
   call msi(ic1_to_ic1b)
   ic1b = 0
   do ic1=1,sdata%nc1
      if (lcheck_nc1b(ic1)) then
         ic1b = ic1b+1
         if (sdata%mpi(iproc)%nc1b>0) ic1b_to_ic1(ic1b) = ic1
         ic1_to_ic1b(ic1) = ic1b
      end if
   end do

   allocate(ic2il1_to_ic2b(sdata%nc1,sdata%nl1))
   call msi(ic2il1_to_ic2b)
   do il1=1,sdata%mpi(iproc)%nl1
      if (sdata%mpi(iproc)%nc2b(il1)>0) then
         ic2b = 0
         do ic2=1,sdata%nc2(il1)
            if (lcheck_nc2b(ic2,il1)) then
               ic2b = ic2b+1
               ic2il1_to_ic2b(ic2,il1) = ic2b
            end if
         end do
      end if
   end do

   if (sdata%mpi(iproc)%nsb>0) allocate(sdata%mpi(iproc)%isb_to_is(sdata%mpi(iproc)%nsb))
   allocate(sdata%mpi(iproc)%is_to_isb(sdata%ns))
   call msi(sdata%mpi(iproc)%is_to_isb)
   isb = 0
   do is=1,sdata%ns
      if (lcheck_nsb(is)) then
         isb = isb+1
         if (sdata%mpi(iproc)%nsb>0) sdata%mpi(iproc)%isb_to_is(isb) = is
         sdata%mpi(iproc)%is_to_isb(is) = isb
      end if
   end do

   ! Halo C
   if (sdata%mpi(iproc)%nsc>0) allocate(sdata%mpi(iproc)%isc_to_is(sdata%mpi(iproc)%nsc))
   allocate(sdata%mpi(iproc)%is_to_isc(sdata%ns))
   call msi(sdata%mpi(iproc)%is_to_isc)
   isc = 0
   do is=1,sdata%ns
      if (lcheck_nsc(is)) then
         isc = isc+1
         if (sdata%mpi(iproc)%nsc>0) sdata%mpi(iproc)%isc_to_is(isc) = is
         sdata%mpi(iproc)%is_to_isc(is) = isc
      end if
   end do

   ! Inter-halo conversions
   if ((sdata%mpi(iproc)%nsa>0).and.(sdata%mpi(iproc)%nsb>0).and.(sdata%mpi(iproc)%nsc>0)) then
      allocate(sdata%mpi(iproc)%isa_to_isb(sdata%mpi(iproc)%nsa))
      allocate(sdata%mpi(iproc)%isa_to_isc(sdata%mpi(iproc)%nsa))
      do isa=1,sdata%mpi(iproc)%nsa
         is = sdata%mpi(iproc)%isa_to_is(isa)
         isb = sdata%mpi(iproc)%is_to_isb(is)
         isc = sdata%mpi(iproc)%is_to_isc(is)
         sdata%mpi(iproc)%isa_to_isb(isa) = isb
         sdata%mpi(iproc)%isa_to_isc(isa) = isc
      end do
      allocate(sdata%mpi(iproc)%isb_to_isc(sdata%mpi(iproc)%nsb))
      do isb=1,sdata%mpi(iproc)%nsb
         is = sdata%mpi(iproc)%isb_to_is(isb)
         isc = sdata%mpi(iproc)%is_to_isc(is)
         sdata%mpi(iproc)%isb_to_isc(isb) = isc
      end do
   end if

   ! Global <-> local conversions for data
   allocate(interph_i_s_lg(sdata%mpi(iproc)%h(1)%n_s,sdata%mpi(iproc)%nl0i))
   do il0i=1,sdata%nl0i
      i_s_loc = 0
      do i_s=1,sdata%h(il0i)%n_s
         if (lcheck_h(i_s,il0i)) then
            i_s_loc = i_s_loc+1
            interph_i_s_lg(i_s_loc,il0i) = i_s
         end if
      end do
   end do
   if (s_n_s_max_loc>0) then
      allocate(interps_i_s_lg(s_n_s_max_loc,sdata%mpi(iproc)%nl1))
      do il1=1,sdata%mpi(iproc)%nl1
         i_s_loc = 0
         do i_s=1,sdata%s(il1)%n_s
            if (lcheck_s(i_s,il1)) then
               i_s_loc = i_s_loc+1
               interps_i_s_lg(i_s_loc,il1) = i_s
            end if
         end do
      end do
   end if
   if (sdata%mpi(iproc)%c%n_s>0) then
      allocate(convol_i_s_lg(sdata%mpi(iproc)%c%n_s))
      i_s_loc = 0
      do i_s=1,sdata%c%n_s
         if (lcheck_c(i_s)) then
            i_s_loc = i_s_loc+1
            convol_i_s_lg(i_s_loc) = i_s
         end if
      end do
   end if

   ! Number of cells
   sdata%mpi(iproc)%nc0a = count(sdata%ic0_to_iproc==iproc)

   ! Local data

   ! Horizontal interpolation
   do il0i=1,sdata%mpi(iproc)%nl0i
      sdata%mpi(iproc)%h(il0i)%prefix = 'h'
      sdata%mpi(iproc)%h(il0i)%n_src = sdata%mpi(iproc)%nc1b
      sdata%mpi(iproc)%h(il0i)%n_dst = sdata%mpi(iproc)%nc0a
      call linop_alloc(sdata%mpi(iproc)%h(il0i))
      do i_s_loc=1,sdata%mpi(iproc)%h(il0i)%n_s
         i_s = interph_i_s_lg(i_s_loc,il0i)
         sdata%mpi(iproc)%h(il0i)%row(i_s_loc) = sdata%ic0_to_ic0a(sdata%h(il0i)%row(i_s))
         sdata%mpi(iproc)%h(il0i)%col(i_s_loc) = ic1_to_ic1b(sdata%h(il0i)%col(i_s))
         sdata%mpi(iproc)%h(il0i)%S(i_s_loc) = sdata%h(il0i)%S(i_s)
      end do
      call linop_reorder(sdata%mpi(iproc)%h(il0i))
   end do

   ! Vertical interpolation
   do il0i=1,sdata%mpi(iproc)%nl0i
      sdata%mpi(iproc)%v(il0i)%prefix = 'v'
      sdata%mpi(iproc)%v(il0i)%n_src = sdata%v(il0i)%n_src
      sdata%mpi(iproc)%v(il0i)%n_dst = sdata%v(il0i)%n_dst
      sdata%mpi(iproc)%v(il0i)%n_s = sdata%v(il0i)%n_s
      if (sdata%mpi(iproc)%v(il0i)%n_s>0) then
         call linop_alloc(sdata%mpi(iproc)%v(il0i))
         sdata%mpi(iproc)%v(il0i)%row = sdata%v(il0i)%row
         sdata%mpi(iproc)%v(il0i)%col = sdata%v(il0i)%col
         sdata%mpi(iproc)%v(il0i)%S = sdata%v(il0i)%S
         call linop_reorder(sdata%mpi(iproc)%v(il0i))
      end if
   end do
   if (sdata%mpi(iproc)%nc1b>0) then
      allocate(sdata%mpi(iproc)%vbot(sdata%mpi(iproc)%nc1b))
      sdata%mpi(iproc)%vbot = sdata%vbot(ic1b_to_ic1)
   end if

   ! Subsampling horizontal interpolation
   do il1=1,sdata%nl1
      sdata%mpi(iproc)%s(il1)%prefix = 's'
      sdata%mpi(iproc)%s(il1)%n_src = sdata%mpi(iproc)%nc2b(il1)
      sdata%mpi(iproc)%s(il1)%n_dst = sdata%mpi(iproc)%nc1b
      if (sdata%mpi(iproc)%s(il1)%n_s>0) then
         call linop_alloc(sdata%mpi(iproc)%s(il1))
         do i_s_loc=1,sdata%mpi(iproc)%s(il1)%n_s
            i_s = interps_i_s_lg(i_s_loc,il1)
            sdata%mpi(iproc)%s(il1)%row(i_s_loc) = ic1_to_ic1b(sdata%s(il1)%row(i_s))
            sdata%mpi(iproc)%s(il1)%col(i_s_loc) = ic2il1_to_ic2b(sdata%s(il1)%col(i_s),il1)
            sdata%mpi(iproc)%s(il1)%S(i_s_loc) = sdata%s(il1)%S(i_s)
         end do
         call linop_reorder(sdata%mpi(iproc)%s(il1))
      end if
   end do

   ! Copy
   if (sdata%mpi(iproc)%nsb>0) then
      allocate(sdata%mpi(iproc)%isb_to_ic2b(sdata%mpi(iproc)%nsb))
      allocate(sdata%mpi(iproc)%isb_to_il1(sdata%mpi(iproc)%nsb))
      call msi(sdata%mpi(iproc)%isb_to_ic2b)
      do isb=1,sdata%mpi(iproc)%nsb
         is = sdata%mpi(iproc)%isb_to_is(isb)
         il1 = sdata%is_to_il1(is)
         ic2 = sdata%is_to_ic2(is)
         ic2b = ic2il1_to_ic2b(ic2,il1)
         sdata%mpi(iproc)%isb_to_ic2b(isb) = ic2b
         sdata%mpi(iproc)%isb_to_il1(isb) = il1
      end do
   end if

   ! Convolution
   if (sdata%mpi(iproc)%c%n_s>0) then
      sdata%mpi(iproc)%c%prefix = 'c'
      sdata%mpi(iproc)%c%n_src = sdata%mpi(iproc)%nsc
      sdata%mpi(iproc)%c%n_dst = sdata%mpi(iproc)%nsc
      call linop_alloc(sdata%mpi(iproc)%c)
      do i_s_loc=1,sdata%mpi(iproc)%c%n_s
         i_s = convol_i_s_lg(i_s_loc)
         sdata%mpi(iproc)%c%row(i_s_loc) = sdata%mpi(iproc)%is_to_isc(sdata%c%row(i_s))
         sdata%mpi(iproc)%c%col(i_s_loc) = sdata%mpi(iproc)%is_to_isc(sdata%c%col(i_s))
         sdata%mpi(iproc)%c%S(i_s_loc) = sdata%c%S(i_s)
      end do
      call linop_reorder(sdata%mpi(iproc)%c)
   end if

   ! Print local parameters
   write(mpl%unit,'(a7,a,i4)') '','Local parameters for processor #',iproc
   write(mpl%unit,'(a10,a,i8)') '','nc0a =      ',sdata%mpi(iproc)%nc0a
   write(mpl%unit,'(a10,a,i8)') '','nc1b =      ',sdata%mpi(iproc)%nc1b
   do il1=1,sdata%mpi(iproc)%nl1
      write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2b(',il1,') =  ',sdata%mpi(iproc)%nc2b(il1)
   end do
   write(mpl%unit,'(a10,a,i8)') '','nsa =       ',sdata%mpi(iproc)%nsa
   write(mpl%unit,'(a10,a,i8)') '','nsb =       ',sdata%mpi(iproc)%nsb
   write(mpl%unit,'(a10,a,i8)') '','nsc =       ',sdata%mpi(iproc)%nsc
   do il0i=1,sdata%mpi(iproc)%nl0i
      write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',sdata%mpi(iproc)%h(il0i)%n_s
   end do
   do il1=1,sdata%mpi(iproc)%nl1
      write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',sdata%mpi(iproc)%s(il1)%n_s
   end do
   write(mpl%unit,'(a10,a,i8)') '','c%n_s =     ',sdata%mpi(iproc)%c%n_s

   if (iproc==sdata%nproc/2) then
      ! Illustration
      allocate(sdata%halo(sdata%nc0))
      sdata%halo = 0
      do i_s=1,sdata%mpi(iproc)%c%n_s
         ic0 = sdata%ic1_to_ic0(sdata%is_to_ic1(sdata%mpi(iproc)%isc_to_is(sdata%mpi(iproc)%c%row(i_s))))
         sdata%halo(ic0) = 1
         ic0 = sdata%ic1_to_ic0(sdata%is_to_ic1(sdata%mpi(iproc)%isc_to_is(sdata%mpi(iproc)%c%col(i_s))))
         sdata%halo(ic0) = 1
      end do
      do i_s=1,sdata%mpi(iproc)%h(1)%n_s
         ic0 = sdata%ic1_to_ic0(ic1b_to_ic1(sdata%mpi(iproc)%h(1)%col(i_s)))
         sdata%halo(ic0) = 2
      end do
      do isa=1,sdata%mpi(iproc)%nsa
         ic0 = sdata%ic1_to_ic0(sdata%is_to_ic1(sdata%mpi(iproc)%isa_to_is(isa)))
         sdata%halo(ic0) = 3
      end do
   end if

   ! Release memory
   deallocate(ic2il1_to_ic2b)
   if (sdata%mpi(iproc)%nc1b>0) then
      deallocate(ic1b_to_ic1)
      deallocate(ic1_to_ic1b)
   end if
   deallocate(interph_i_s_lg)
   if (s_n_s_max_loc>0) deallocate(interps_i_s_lg)
   if (sdata%mpi(iproc)%c%n_s>0) deallocate(convol_i_s_lg)
end do

! Norm distribution
call fld_com_gl(sdata,sdata%norm,norm)
do iproc=1,sdata%nproc
   allocate(sdata%mpi(iproc)%norm%vala(sdata%mpi(iproc)%nc0a,sdata%mpi(iproc)%nl0))
   sdata%mpi(iproc)%norm%vala = norm(iproc)%vala
end do

! Communications setup
do iproc=1,sdata%nproc
   do isb=1,sdata%mpi(iproc)%nsb
      ! Check for points that are in zone B but are not in zone A
      is = sdata%mpi(iproc)%isb_to_is(isb)
      ic1 = sdata%is_to_ic1(is)
      ic0 = sdata%ic1_to_ic0(ic1)
      jproc = sdata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         sdata%mpi(iproc)%AB%jhalocounts(jproc) = sdata%mpi(iproc)%AB%jhalocounts(jproc)+1

         ! Count of points received on JPROC from IPROC
         sdata%mpi(jproc)%AB%jexclcounts(iproc) = sdata%mpi(jproc)%AB%jexclcounts(iproc)+1
      end if
   end do
   do isc=1,sdata%mpi(iproc)%nsc
      ! Check for points that are in zone C but are not in zone A
      is = sdata%mpi(iproc)%isc_to_is(isc)
      ic1 = sdata%is_to_ic1(is)
      ic0 = sdata%ic1_to_ic0(ic1)
      jproc = sdata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         sdata%mpi(iproc)%AC%jhalocounts(jproc) = sdata%mpi(iproc)%AC%jhalocounts(jproc)+1

         ! Count of points received on JPROC from IPROC
         sdata%mpi(jproc)%AC%jexclcounts(iproc) = sdata%mpi(jproc)%AC%jexclcounts(iproc)+1
      end if
   end do
end do

! Compute displacement
do iproc=1,sdata%nproc
   sdata%mpi(iproc)%AB%jhalodispl(1) = 0
   sdata%mpi(iproc)%AB%jexcldispl(1) = 0
   sdata%mpi(iproc)%AC%jhalodispl(1) = 0
   sdata%mpi(iproc)%AC%jexcldispl(1) = 0
   do jproc=2,sdata%nproc
      sdata%mpi(iproc)%AB%jhalodispl(jproc) = sdata%mpi(iproc)%AB%jhalodispl(jproc-1)+sdata%mpi(iproc)%AB%jhalocounts(jproc-1)
      sdata%mpi(iproc)%AB%jexcldispl(jproc) = sdata%mpi(iproc)%AB%jexcldispl(jproc-1)+sdata%mpi(iproc)%AB%jexclcounts(jproc-1)
      sdata%mpi(iproc)%AC%jhalodispl(jproc) = sdata%mpi(iproc)%AC%jhalodispl(jproc-1)+sdata%mpi(iproc)%AC%jhalocounts(jproc-1)
      sdata%mpi(iproc)%AC%jexcldispl(jproc) = sdata%mpi(iproc)%AC%jexcldispl(jproc-1)+sdata%mpi(iproc)%AC%jexclcounts(jproc-1)
   end do
end do

! Allocation
do iproc=1,sdata%nproc
   sdata%mpi(iproc)%AB%nhalo = sum(sdata%mpi(iproc)%AB%jhalocounts)
   sdata%mpi(iproc)%AB%nexcl = sum(sdata%mpi(iproc)%AB%jexclcounts)
   sdata%mpi(iproc)%AC%nhalo = sum(sdata%mpi(iproc)%AC%jhalocounts)
   sdata%mpi(iproc)%AC%nexcl = sum(sdata%mpi(iproc)%AC%jexclcounts)
   allocate(sdata%mpi(iproc)%AB%halo(sdata%mpi(iproc)%AB%nhalo))
   allocate(sdata%mpi(iproc)%AB%excl(sdata%mpi(iproc)%AB%nexcl))
   allocate(sdata%mpi(iproc)%AC%halo(sdata%mpi(iproc)%AC%nhalo))
   allocate(sdata%mpi(iproc)%AC%excl(sdata%mpi(iproc)%AC%nexcl))
end do

! Fill halo array
do iproc=1,sdata%nproc
   sdata%mpi(iproc)%AB%jhalocounts = 0
   sdata%mpi(iproc)%AC%jhalocounts = 0
end do
do iproc=1,sdata%nproc
   do isb=1,sdata%mpi(iproc)%nsb
      ! Check for points that ARE in zone B but ARE NOT in zone A
      is = sdata%mpi(iproc)%isb_to_is(isb)
      ic1 = sdata%is_to_ic1(is)
      ic0 = sdata%ic1_to_ic0(ic1)
      jproc = sdata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         sdata%mpi(iproc)%AB%jhalocounts(jproc) = sdata%mpi(iproc)%AB%jhalocounts(jproc)+1

         ! Global index of points sent from IPROC to JPROC
         sdata%mpi(iproc)%AB%halo(sdata%mpi(iproc)%AB%jhalodispl(jproc)+sdata%mpi(iproc)%AB%jhalocounts(jproc)) = is
      end if
   end do
   do isc=1,sdata%mpi(iproc)%nsc
      ! Check for points that ARE in zone C but ARE NOT in zone A
      is = sdata%mpi(iproc)%isc_to_is(isc)
      ic1 = sdata%is_to_ic1(is)
      ic0 = sdata%ic1_to_ic0(ic1)
      jproc = sdata%ic0_to_iproc(ic0)
      if (jproc/=iproc) then
         ! Count of points sent from IPROC to JPROC
         sdata%mpi(iproc)%AC%jhalocounts(jproc) = sdata%mpi(iproc)%AC%jhalocounts(jproc)+1

         ! Global index of points sent from IPROC to JPROC
         sdata%mpi(iproc)%AC%halo(sdata%mpi(iproc)%AC%jhalodispl(jproc)+sdata%mpi(iproc)%AC%jhalocounts(jproc)) = is
      end if
   end do
end do

! Fill excl array
do jproc=1,sdata%nproc
   ! Loop over processors sending data to JPROC
   do iproc=1,sdata%nproc
      do i=1,sdata%mpi(iproc)%AB%jhalocounts(jproc)
         ! Global index of points received on JPROC from IPROC
         sdata%mpi(jproc)%AB%excl(sdata%mpi(jproc)%AB%jexcldispl(iproc)+i) = &
       & sdata%mpi(iproc)%AB%halo(sdata%mpi(iproc)%AB%jhalodispl(jproc)+i)
      end do
      do i=1,sdata%mpi(iproc)%AC%jhalocounts(jproc)
         ! Global index of points received on JPROC from IPROC
         sdata%mpi(jproc)%AC%excl(sdata%mpi(jproc)%AC%jexcldispl(iproc)+i) = &
       & sdata%mpi(iproc)%AC%halo(sdata%mpi(iproc)%AC%jhalodispl(jproc)+i)
      end do
   end do
end do

! Transform to local indices
do iproc=1,sdata%nproc
   if ((sdata%mpi(iproc)%nsa>0).and.(sdata%mpi(iproc)%nsb>0).and.(sdata%mpi(iproc)%nsc>0)) then
      sdata%mpi(iproc)%AB%halo = sdata%mpi(iproc)%is_to_isb(sdata%mpi(iproc)%AB%halo)
      sdata%mpi(iproc)%AB%excl = sdata%mpi(iproc)%is_to_isa(sdata%mpi(iproc)%AB%excl)
      sdata%mpi(iproc)%AC%halo = sdata%mpi(iproc)%is_to_isc(sdata%mpi(iproc)%AC%halo)
      sdata%mpi(iproc)%AC%excl = sdata%mpi(iproc)%is_to_isa(sdata%mpi(iproc)%AC%excl)
   end if
end do

end subroutine compute_mpi

end module module_mpi
