!----------------------------------------------------------------------
! Module: module_parameters_interp.f90
!> Purpose: compute NICAS parameters (interpolation)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_parameters_interp

use module_namelist, only: nam
use netcdf
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_nc, only: ncfloat,ncerr
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_esmf, only: esmf_create_field,esmf_create_interp
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_reorder
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_randgen, only: initialize_sampling,rand_integer
use type_sdata, only: sdatatype
implicit none

real(kind_real),parameter :: S_inf = 1.0e-3 !< Minimum value for the interpolation coefficients

private
public :: compute_interp_h,compute_interp_v,compute_interp_s

contains

!----------------------------------------------------------------------
! Subroutine: compute_interp_h
!> Purpose: compute basic horizontal interpolation
!----------------------------------------------------------------------
subroutine compute_interp_h(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ic0,ic1,i_s,ibnd,jc0,jc1,progint,il0i
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
real(kind_real) :: renorm(sdata%nc0)
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:),valid(:),validg(:,:),missing(:)
type(linoptype) :: hbase,htmp
type(ctreetype) :: ctree

! Create ESMF field for C1
call esmf_create_field(sdata,sdata%nc1,sdata%lon(sdata%ic1_to_ic0),sdata%lat(sdata%ic1_to_ic0), &
 & any(sdata%mask(sdata%ic1_to_ic0,:),dim=2),sdata%c1field)

! Compute interpolation
call esmf_create_interp(sdata%c1field,sdata%c0field,sdata%regional,hbase)

! Allocation
allocate(sdata%h(sdata%nl0i))

do il0i=1,sdata%nl0i
   ! Allocation and copy
   htmp%n_s = hbase%n_s
   call linop_alloc(htmp)
   htmp%row = hbase%row
   htmp%col = hbase%col
   htmp%S = hbase%S
   allocate(valid(htmp%n_s))

   ! Check interpolation coefficient
   valid = (htmp%S>S_inf)

   ! Check mask
   do i_s=1,htmp%n_s
      if (valid(i_s)) then
         ic0 = htmp%row(i_s)
         jc1 = htmp%col(i_s)
         jc0 = sdata%ic1_to_ic0(jc1)
         valid(i_s) = sdata%mask(ic0,il0i).and.sdata%mask(jc0,il0i)
      end if
   end do

   if (nam%mask_check) then
      ! MPI splitting
      do iproc=1,mpl%nproc
         i_s_s(iproc) = (iproc-1)*(htmp%n_s/mpl%nproc+1)+1
         i_s_e(iproc) = min(iproc*(htmp%n_s/mpl%nproc+1),htmp%n_s)
         n_s_loc(iproc) = i_s_e(iproc)-i_s_s(iproc)+1
      end do

      ! Allocation
      allocate(done(n_s_loc(mpl%myproc)))

      ! Check that interpolations are not crossing mask boundaries
      write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0i),': '
      call prog_init(progint,done)
      !$omp parallel do private(i_s_loc,i_s,x,y,z,v1,v2,va,vp,t,ic0,jc1,jc0)
      do i_s_loc=1,n_s_loc(mpl%myproc)
         ! Indices
         i_s = i_s_s(mpl%myproc)+i_s_loc-1

         if (valid(i_s)) then
            ! Allocation
            allocate(x(2))
            allocate(y(2))
            allocate(z(2))
            allocate(v1(3))
            allocate(v2(3))
            allocate(va(3))
            allocate(vp(3))
            allocate(t(4))

            ! Indices
            ic0 = htmp%row(i_s)
            jc1 = htmp%col(i_s)
            jc0 = sdata%ic1_to_ic0(jc1)

            ! Transform to cartesian coordinates
            call trans(2,sdata%lat((/ic0,jc0/)),sdata%lon((/ic0,jc0/)),x,y,z)

            ! Compute arc orthogonal vector
            v1 = (/x(1),y(1),z(1)/)
            v2 = (/x(2),y(2),z(2)/)
            call vector_product(v1,v2,va)

            ! Check if arc is crossing boundary arcs
            do ibnd=1,sdata%nbnd(il0i)
               call vector_product(va,sdata%vbnd(:,ibnd,il0i),vp)
               v1 = (/x(1),y(1),z(1)/)
               call vector_triple_product(v1,va,vp,t(1))
               v1 = (/x(2),y(2),z(2)/)
               call vector_triple_product(v1,va,vp,t(2))
               v1 = (/sdata%xbnd(1,ibnd,il0i),sdata%ybnd(1,ibnd,il0i),sdata%zbnd(1,ibnd,il0i)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0i),vp,t(3))
               v1 = (/sdata%xbnd(2,ibnd,il0i),sdata%ybnd(2,ibnd,il0i),sdata%zbnd(2,ibnd,il0i)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0i),vp,t(4))
               t(1) = -t(1)
               t(3) = -t(3)
               if (all(t>0).or.(all(t<0))) then
                  valid(i_s) = .false.
                  exit
               end if
            end do

            ! Memory release
            deallocate(x)
            deallocate(y)
            deallocate(z)
            deallocate(v1)
            deallocate(v2)
            deallocate(va)
            deallocate(vp)
            deallocate(t)
         end if

         ! Print progression
         done(i_s_loc) = .true.
         call prog_print(progint,done)
      end do
      !$omp end parallel do
      write(mpl%unit,'(a)') '100%'

      ! Communication
      if (mpl%main) then
         ! Allocation
         allocate(validg(maxval(n_s_loc),mpl%nproc))

         do iproc=1,mpl%nproc
            if (iproc==mpl%ioproc) then
               ! Copy data
               validg(1:n_s_loc(iproc),iproc) = valid(i_s_s(iproc):i_s_e(iproc))
            else
               ! Receive data on ioproc
               call mpl_recv(n_s_loc(iproc),validg(1:n_s_loc(iproc),iproc),iproc,mpl%tag)
            end if
         end do

         ! Format data
         do iproc=1,mpl%nproc
            valid(i_s_s(iproc):i_s_e(iproc)) = validg(1:n_s_loc(iproc),iproc)
         end do

         ! Release memory
         deallocate(validg)
      else
         ! Send data to ioproc
         call mpl_send(n_s_loc(mpl%myproc),valid(i_s_s(mpl%myproc):i_s_e(mpl%myproc)),mpl%ioproc,mpl%tag)
      end if
      mpl%tag = mpl%tag+1

      ! Broadcast
      call mpl_bcast(valid,mpl%ioproc)

      ! Release memory
      deallocate(done)
   else
      write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(il0i)
   end if

   ! Renormalization
   renorm = 0.0
   do i_s=1,htmp%n_s
      if (valid(i_s)) renorm(htmp%row(i_s)) = renorm(htmp%row(i_s))+htmp%S(i_s)
   end do

   ! Initialize object
   sdata%h(il0i)%prefix = 'h'
   sdata%h(il0i)%n_src = sdata%nc1
   sdata%h(il0i)%n_dst = sdata%nc0
   sdata%h(il0i)%n_s = count(valid)
   call linop_alloc(sdata%h(il0i))
   sdata%h(il0i)%n_s = 0
   do i_s=1,htmp%n_s
      if (valid(i_s)) then
         sdata%h(il0i)%n_s = sdata%h(il0i)%n_s+1
         sdata%h(il0i)%row(sdata%h(il0i)%n_s) = htmp%row(i_s)
         sdata%h(il0i)%col(sdata%h(il0i)%n_s) = htmp%col(i_s)
         sdata%h(il0i)%S(sdata%h(il0i)%n_s) = htmp%S(i_s)/renorm(htmp%row(i_s))
      end if
   end do

   ! Release memory
   call linop_dealloc(htmp)
   deallocate(valid)

   ! Allocation
   allocate(missing(sdata%nc0))

   ! Count points that are not interpolated
   missing = .false.
   do ic0=1,sdata%nc0
      if (sdata%mask(ic0,il0i)) missing(ic0) = .true.
   end do
   do i_s=1,sdata%h(il0i)%n_s
      missing(sdata%h(il0i)%row(i_s)) = .false.
   end do
   if (count(missing)>0) then
      ! Allocate temporary arrays
      htmp%n_s = sdata%h(il0i)%n_s
      call linop_alloc(htmp)

      ! Fill temporary arrays
      htmp%row = sdata%h(il0i)%row
      htmp%col = sdata%h(il0i)%col
      htmp%S = sdata%h(il0i)%S

      ! Reallocate interpolation
      call linop_dealloc(sdata%h(il0i))
      sdata%h(il0i)%n_s = sdata%h(il0i)%n_s+count(missing)
      call linop_alloc(sdata%h(il0i))

      ! Fill permanent arrays
      sdata%h(il0i)%row(1:htmp%n_s) = htmp%row
      sdata%h(il0i)%col(1:htmp%n_s) = htmp%col
      sdata%h(il0i)%S(1:htmp%n_s) = htmp%S

      ! Compute cover tree
      allocate(mask_ctree(sdata%nc1))
      mask_ctree = 0
      do ic1=1,sdata%nc1
         if (sdata%mask(sdata%ic1_to_ic0(ic1),il0i)) mask_ctree(ic1) = 1
      end do
      ctree = create_ctree(sdata%nc1,dble(sdata%lon(sdata%ic1_to_ic0)),dble(sdata%lat(sdata%ic1_to_ic0)),mask_ctree)
      deallocate(mask_ctree)

      ! Compute nearest neighbors
      do ic0=1,sdata%nc0
         if (missing(ic0)) then
            htmp%n_s = htmp%n_s+1
            sdata%h(il0i)%row(htmp%n_s) = ic0
            call find_nearest_neighbors(ctree,dble(sdata%lon(ic0)),dble(sdata%lat(ic0)),1,sdata%h(il0i)%col(htmp%n_s:htmp%n_s),dum)
            sdata%h(il0i)%S(htmp%n_s) = 1.0
         end if
      end do

      ! Release memory
      call linop_dealloc(htmp)
      call delete_ctree(ctree)
   end if

   ! Reorder linear operator
   call linop_reorder(sdata%h(il0i))

   ! Release memory
   deallocate(missing)
end do

! Release memory
call linop_dealloc(hbase)

end subroutine compute_interp_h

!----------------------------------------------------------------------
! Subroutine: compute_interp_v
!> Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine compute_interp_v(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: il0,jl0,il1,il0inf,il0sup,il0i,i_s,ic0,ic1
logical,allocatable :: valid(:)
type(linoptype) :: vbase,vtmp

! Linear interpolation
vbase%n_s = sdata%nl1
il0inf = 1
do il0=1,sdata%nl0
   if (sdata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         vbase%n_s = vbase%n_s+2
      end do
      il0inf = il0
   end if
end do
call linop_alloc(vbase)
do il1=1,sdata%nl1
   il0 = sdata%il1_to_il0(il1)
   vbase%row(il1) = il0
   vbase%col(il1) = il0
   vbase%S(il1) = 1.0
end do
vbase%n_s = sdata%nl1
il0inf = 1
do il0=1,sdata%nl0
   if (sdata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         vbase%n_s = vbase%n_s+1
         vbase%row(vbase%n_s) = jl0
         vbase%col(vbase%n_s) = il0inf
         vbase%S(vbase%n_s) = abs(sdata%vunit(il0sup)-sdata%vunit(jl0))/abs(sdata%vunit(il0sup)-sdata%vunit(il0inf))

         vbase%n_s = vbase%n_s+1
         vbase%row(vbase%n_s) = jl0
         vbase%col(vbase%n_s) = il0sup
         vbase%S(vbase%n_s) = abs(sdata%vunit(jl0)-sdata%vunit(il0inf))/abs(sdata%vunit(il0sup)-sdata%vunit(il0inf))
      end do
      il0inf = il0
   end if
end do

! Allocation
allocate(sdata%v(sdata%nl0i))
allocate(valid(vbase%n_s))

do il0i=1,sdata%nl0i
   ! Initialize vertical interpolation
   sdata%v(il0i)%prefix = 'v'
   sdata%v(il0i)%n_src = sdata%nl1
   sdata%v(il0i)%n_dst = sdata%nl0

   if (sdata%nl0i==1) then
      ! Copy basic vertical interpolation
      sdata%v(il0i)%n_s = vbase%n_s
      call linop_alloc(sdata%v(il0i))
      sdata%v(il0i)%row = vbase%row
      sdata%v(il0i)%col = vbase%col
      sdata%v(il0i)%S = vbase%S
   else
      ! Check valid operations
      do i_s=1,vbase%n_s
         valid(i_s) = (vbase%row(i_s)<=il0i).and.(vbase%col(i_s)<=il0i)
      end do

      ! Copy valid operations
      vtmp%n_s = count(valid)
      call linop_alloc(vtmp)
      vtmp%n_s = 0
      do i_s=1,vbase%n_s
         if (valid(i_s)) then
            vtmp%n_s = vtmp%n_s+1
            vtmp%row(vtmp%n_s) = vbase%row(i_s)
            vtmp%col(vtmp%n_s) = vbase%col(i_s)
            vtmp%S(vtmp%n_s) = vbase%S(i_s)
         end if
      end do

      ! Add missing levels to the interpolation
      sdata%v(il0i)%n_s = vtmp%n_s+(il0i-maxval(vtmp%col))
      call linop_alloc(sdata%v(il0i))
      sdata%v(il0i)%row(1:vtmp%n_s) = vtmp%row
      sdata%v(il0i)%col(1:vtmp%n_s) = vtmp%col
      sdata%v(il0i)%S(1:vtmp%n_s) = vtmp%S
      do jl0=maxval(vtmp%col)+1,il0i
         vtmp%n_s = vtmp%n_s+1
         sdata%v(il0i)%row(vtmp%n_s) = jl0
         sdata%v(il0i)%col(vtmp%n_s) = jl0
         sdata%v(il0i)%S(vtmp%n_s) = 1.0
      end do

      ! Release memory
      call linop_dealloc(vtmp)
   end if

   ! Conversion
   sdata%v(il0i)%col = sdata%il0_to_il1(sdata%v(il0i)%col)

   ! Reorder linear operator
   call linop_reorder(sdata%v(il0i))
end do

! Release memory
call linop_dealloc(vbase)
deallocate(valid)

! Find the bottom for each point of S1
allocate(sdata%vbot(sdata%nc1))
!$omp parallel do private(ic1,ic0,il0)
do ic1=1,sdata%nc1
   ic0 = sdata%ic1_to_ic0(ic1)
   il0 = 1
   do while (sdata%mask(ic0,il0).and.(il0<sdata%nl0))
      il0 = il0+1
   end do
   sdata%vbot(ic1) = min(il0,sdata%nl0i)
end do
!$omp end parallel do

end subroutine compute_interp_v

!----------------------------------------------------------------------
! Subroutine: compute_interp_s
!> Purpose: compute horizontal subsampling interpolation
!----------------------------------------------------------------------
subroutine compute_interp_s(sdata)

implicit none

! Passed variables
type(sdatatype),intent(inout) :: sdata !< Sampling data

! Local variables
integer :: ic1,ic2,il1,i_s,ibnd,il0,ic0,jc2,jc1,jc0,progint
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
real(kind_real) :: renorm(sdata%nc1)
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:),valid(:),validg(:,:),missing(:)
type(linoptype) :: stmp
type(ctreetype) :: ctree

! Allocation
allocate(sdata%c2field(sdata%nl1))
allocate(sdata%s(sdata%nl1))

do il1=1,sdata%nl1
   ! Create ESMF field for C2
   call esmf_create_field(sdata,sdata%nc2(il1),sdata%lon(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1)), & 
 & sdata%lat(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1)),sdata%mask(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1), &
 & sdata%il1_to_il0(il1)),sdata%c2field(il1))

   ! Compute interpolation
   call esmf_create_interp(sdata%c2field(il1),sdata%c1field,sdata%regional,interp=stmp)

   ! Allocation
   allocate(valid(stmp%n_s))

   ! Check interpolation coefficient
   valid = .not.(stmp%S<S_inf)

   if (nam%mask_check) then
      ! MPI splitting
      do iproc=1,mpl%nproc
         i_s_s(iproc) = (iproc-1)*(stmp%n_s/mpl%nproc+1)+1
         i_s_e(iproc) = min(iproc*(stmp%n_s/mpl%nproc+1),stmp%n_s)
         n_s_loc(iproc) = i_s_e(iproc)-i_s_s(iproc)+1
      end do

      ! Allocation
      allocate(done(n_s_loc(mpl%myproc)))

      ! Check that interpolations are not crossing mask boundaries
      write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(sdata%il1_to_il0(il1)),': '
      call prog_init(progint,done)
      !$omp parallel do private(i_s_loc,i_s,x,y,z,v1,v2,va,vp,t,ic1,ic0,jc2,jc1,jc0,il0)
      do i_s_loc=1,n_s_loc(mpl%myproc)
         ! Indices
         i_s = i_s_s(mpl%myproc)+i_s_loc-1

         if (valid(i_s)) then
            ! Allocation
            allocate(x(2))
            allocate(y(2))
            allocate(z(2))
            allocate(v1(3))
            allocate(v2(3))
            allocate(va(3))
            allocate(vp(3))
            allocate(t(4))

            ! Indices
            ic1 = stmp%row(i_s)
            ic0 = sdata%ic1_to_ic0(ic1)
            jc2 = stmp%col(i_s)
            jc1 = sdata%ic2il1_to_ic1(jc2,il1)
            jc0 = sdata%ic1_to_ic0(jc1)
            il0 = sdata%il1_to_il0(il1)

            ! Transform to cartesian coordinates
            call trans(2,sdata%lat((/ic0,jc0/)),sdata%lon((/ic0,jc0/)),x,y,z)

            ! Compute arc orthogonal vector
            v1 = (/x(1),y(1),z(1)/)
            v2 = (/x(2),y(2),z(2)/)
            call vector_product(v1,v2,va)

            ! Check if arc is crossing boundary arcs
            do ibnd=1,sdata%nbnd(il0)
               call vector_product(va,sdata%vbnd(:,ibnd,il0),vp)
               v1 = (/x(1),y(1),z(1)/)
               call vector_triple_product(v1,va,vp,t(1))
               v1 = (/x(2),y(2),z(2)/)
               call vector_triple_product(v1,va,vp,t(2))
               v1 = (/sdata%xbnd(1,ibnd,il0),sdata%ybnd(1,ibnd,il0),sdata%zbnd(1,ibnd,il0)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0),vp,t(3))
               v1 = (/sdata%xbnd(2,ibnd,il0),sdata%ybnd(2,ibnd,il0),sdata%zbnd(2,ibnd,il0)/)
               call vector_triple_product(v1,sdata%vbnd(:,ibnd,il0),vp,t(4))
               t(1) = -t(1)
               t(3) = -t(3)
               if (all(t>0).or.(all(t<0))) then
                  valid(i_s) = .false.
                  exit
               end if
            end do

            ! Memory release
            deallocate(x)
            deallocate(y)
            deallocate(z)
            deallocate(v1)
            deallocate(v2)
            deallocate(va)
            deallocate(vp)
            deallocate(t)
         end if

         ! Print progression
         done(i_s_loc) = .true.
         call prog_print(progint,done)
      end do
      !$omp end parallel do
      write(mpl%unit,'(a)') '100%'

      ! Communication
      if (mpl%main) then
         ! Allocation
         allocate(validg(maxval(n_s_loc),mpl%nproc))

         do iproc=1,mpl%nproc
            if (iproc==mpl%ioproc) then
               ! Copy data
               validg(1:n_s_loc(iproc),iproc) = valid(i_s_s(iproc):i_s_e(iproc))
            else
               ! Receive data on ioproc
               call mpl_recv(n_s_loc(iproc),validg(1:n_s_loc(iproc),iproc),iproc,mpl%tag)
            end if
         end do

         ! Format data
         do iproc=1,mpl%nproc
            valid(i_s_s(iproc):i_s_e(iproc)) = validg(1:n_s_loc(iproc),iproc)
         end do

         ! Release memory
         deallocate(validg)
      else
         ! Send data to ioproc
         call mpl_send(n_s_loc(mpl%myproc),valid(i_s_s(mpl%myproc):i_s_e(mpl%myproc)),mpl%ioproc,mpl%tag)
      end if
      mpl%tag = mpl%tag+1

      ! Broadcast
      call mpl_bcast(valid,mpl%ioproc)

      ! Release memory
      deallocate(done)
   else
      write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(sdata%il1_to_il0(il1))
   end if

   ! Renormalization
   renorm = 0.0
   do i_s=1,stmp%n_s
      if (valid(i_s)) renorm(stmp%row(i_s)) = renorm(stmp%row(i_s))+stmp%S(i_s)
   end do

   ! Initialize object
   sdata%s(il1)%prefix = 's'
   sdata%s(il1)%n_src = sdata%nc2(il1)
   sdata%s(il1)%n_dst = sdata%nc1
   sdata%s(il1)%n_s = count(valid)
   call linop_alloc(sdata%s(il1))
   sdata%s(il1)%n_s = 0
   do i_s=1,stmp%n_s
      if (valid(i_s)) then
         sdata%s(il1)%n_s = sdata%s(il1)%n_s+1
         sdata%s(il1)%row(sdata%s(il1)%n_s) = stmp%row(i_s)
         sdata%s(il1)%col(sdata%s(il1)%n_s) = stmp%col(i_s)
         sdata%s(il1)%S(sdata%s(il1)%n_s) = stmp%S(i_s)/renorm(stmp%row(i_s))
      end if
   end do

   ! Release memory
   call linop_dealloc(stmp)
   deallocate(valid)

   ! Allocation
   allocate(missing(sdata%nc1))

   ! Count points that are not interpolated
   missing = .false.
   do ic1=1,sdata%nc1
      if (sdata%mask(sdata%ic1_to_ic0(ic1),sdata%il1_to_il0(il1))) missing(ic1) = .true.
   end do
   do i_s=1,sdata%s(il1)%n_s
      missing(sdata%s(il1)%row(i_s)) = .false.
   end do
   if (count(missing)>0) then
      ! Allocate temporary arrays
      stmp%n_s = sdata%s(il1)%n_s
      call linop_alloc(stmp)

      ! Fill temporary arrays
      stmp%row = sdata%s(il1)%row
      stmp%col = sdata%s(il1)%col
      stmp%S = sdata%s(il1)%S

      ! Reallocate permanent arrays
      call linop_dealloc(sdata%s(il1))
      sdata%s(il1)%n_s = sdata%s(il1)%n_s+count(missing)
      call linop_alloc(sdata%s(il1))

      ! Fill permanent arrays
      sdata%s(il1)%row(1:stmp%n_s) = stmp%row
      sdata%s(il1)%col(1:stmp%n_s) = stmp%col
      sdata%s(il1)%S(1:stmp%n_s) = stmp%S

      ! Compute cover tree
      allocate(mask_ctree(sdata%nc2(il1)))
      mask_ctree = 0
      do ic2=1,sdata%nc2(il1)
         if (sdata%mask(sdata%ic1_to_ic0(sdata%ic2il1_to_ic1(ic2,il1)),sdata%il1_to_il0(il1))) mask_ctree(ic2) = 1
      end do
      ctree = create_ctree(sdata%nc2(il1),dble(sdata%lon(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1))), &
    & dble(sdata%lat(sdata%ic2il1_to_ic0(1:sdata%nc2(il1),il1))),mask_ctree)
      deallocate(mask_ctree)

      ! Compute nearest neighbors
      do ic1=1,sdata%nc1
         if (missing(ic1)) then
            stmp%n_s = stmp%n_s+1
            sdata%s(il1)%row(stmp%n_s) = ic1
            call find_nearest_neighbors(ctree,dble(sdata%lon(sdata%ic1_to_ic0(ic1))), &
          & dble(sdata%lat(sdata%ic1_to_ic0(ic1))),1,sdata%s(il1)%col(stmp%n_s:stmp%n_s),dum)
            sdata%s(il1)%S(stmp%n_s) = 1.0
         end if
      end do

      ! Release memory
      call linop_dealloc(stmp)
      call delete_ctree(ctree)
   end if

   ! Reorder linear operator
   call linop_reorder(sdata%s(il1))

   ! Release memory
   deallocate(missing)
end do

end subroutine compute_interp_s

end module module_parameters_interp
