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

use module_namelist, only: namtype
use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist,vector_product,vector_triple_product
use tools_display, only: msgerror,prog_init,prog_print
use tools_interp, only: interp_horiz
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_stripack, only: trans
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_reorder
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_ndata, only: ndatatype
use type_randgen, only: rng,initialize_sampling

implicit none

private
public :: compute_interp_h,compute_interp_v,compute_interp_s

contains

!----------------------------------------------------------------------
! Subroutine: compute_interp_h
!> Purpose: compute basic horizontal interpolation
!----------------------------------------------------------------------
subroutine compute_interp_h(nam,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ic0,ic1,i_s,il0i
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
real(kind_real) :: renorm(ndata%nc0)
logical,allocatable :: mask_extra(:),valid(:),missing(:)
type(linoptype) :: hbase,htmp
type(ctreetype) :: ctree

! Compute interpolation
call interp_horiz(rng,ndata%nc1,ndata%geom%lon(ndata%ic1_to_ic0),ndata%geom%lat(ndata%ic1_to_ic0), &
 & any(ndata%geom%mask(ndata%ic1_to_ic0,:),dim=2),ndata%nc0,ndata%geom%lon,ndata%geom%lat, &
 & any(ndata%geom%mask,dim=2),hbase) 

! Allocation
allocate(valid(hbase%n_s))
allocate(mask_extra(ndata%nc1))
allocate(ndata%h(ndata%nl0i))

do il0i=1,ndata%nl0i
   ! Copy
   call linop_copy(hbase,htmp)
   valid = .true.

   ! Check mask boundaries
   if (nam%mask_check) then
      call check_mask_bnd(nam,ndata,htmp,valid,il0i,col_to_ic0=ndata%ic1_to_ic0)
   else
      write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(il0i)
   end if

   if (ndata%nl0i>1) then
      ! Check extrapolated points because of masks vertical variations (unsampled levels)
      mask_extra = .false.
      do ic1=1,ndata%nc1
         ic0 = ndata%ic1_to_ic0(ic1)
         if (ndata%geom%mask(ic0,il0i).and.((il0i<ndata%vbot(ic1)).or.(il0i>ndata%vtop(ic1)))) mask_extra(ic1) = .true.
      end do

      ! Remove operations for extrapolated points
      do i_s=1,htmp%n_s
         if (valid(i_s)) then
            ic1 = hbase%col(i_s)
            if (mask_extra(ic1)) valid(i_s) = .false.
         end if
      end do
      if (count(mask_extra)>0) write(mpl%unit,'(a10,a,i5)') '','Extrapolated points: ',count(mask_extra)
   end if

   ! Renormalization
   renorm = 0.0
   do i_s=1,htmp%n_s
      if (valid(i_s)) renorm(htmp%row(i_s)) = renorm(htmp%row(i_s))+htmp%S(i_s)
   end do

   ! Initialize object
   ndata%h(il0i)%prefix = 'h'
   ndata%h(il0i)%n_src = ndata%nc1
   ndata%h(il0i)%n_dst = ndata%nc0
   ndata%h(il0i)%n_s = count(valid)
   call linop_alloc(ndata%h(il0i))
   ndata%h(il0i)%n_s = 0
   do i_s=1,htmp%n_s
      if (valid(i_s)) then
         ndata%h(il0i)%n_s = ndata%h(il0i)%n_s+1
         ndata%h(il0i)%row(ndata%h(il0i)%n_s) = htmp%row(i_s)
         ndata%h(il0i)%col(ndata%h(il0i)%n_s) = htmp%col(i_s)
         ndata%h(il0i)%S(ndata%h(il0i)%n_s) = htmp%S(i_s)/renorm(htmp%row(i_s))
      end if
   end do

   ! Release memory
   call linop_dealloc(htmp)

   ! Allocation
   allocate(missing(ndata%nc0))

   ! Count points that are not interpolated
   missing = .false.
   do ic0=1,ndata%nc0
      if (ndata%geom%mask(ic0,il0i)) missing(ic0) = .true.
   end do
   do i_s=1,ndata%h(il0i)%n_s
      missing(ndata%h(il0i)%row(i_s)) = .false.
   end do
   if (count(missing)>0) then
      ! Copy
      call linop_copy(ndata%h(il0i),htmp)

      ! Reallocate interpolation
      call linop_dealloc(ndata%h(il0i))
      ndata%h(il0i)%n_s = ndata%h(il0i)%n_s+count(missing)
      call linop_alloc(ndata%h(il0i))

      ! Fill permanent arrays
      ndata%h(il0i)%row(1:htmp%n_s) = htmp%row
      ndata%h(il0i)%col(1:htmp%n_s) = htmp%col
      ndata%h(il0i)%S(1:htmp%n_s) = htmp%S

      ! Compute cover tree
      allocate(mask_ctree(ndata%nc1))
      mask_ctree = 0
      do ic1=1,ndata%nc1
         if (ndata%geom%mask(ndata%ic1_to_ic0(ic1),il0i).and.(.not.mask_extra(ic1))) mask_ctree(ic1) = 1
      end do
      ctree = create_ctree(ndata%nc1,dble(ndata%geom%lon(ndata%ic1_to_ic0)),dble(ndata%geom%lat(ndata%ic1_to_ic0)),mask_ctree)
      deallocate(mask_ctree)

      ! Compute nearest neighbors
      do ic0=1,ndata%nc0
         if (missing(ic0)) then
            htmp%n_s = htmp%n_s+1
            ndata%h(il0i)%row(htmp%n_s) = ic0
            call find_nearest_neighbors(ctree,dble(ndata%geom%lon(ic0)),dble(ndata%geom%lat(ic0)),1, &
          & ndata%h(il0i)%col(htmp%n_s:htmp%n_s),dum)
            ndata%h(il0i)%S(htmp%n_s) = 1.0
         end if
      end do

      ! Release memory
      call linop_dealloc(htmp)
      call delete_ctree(ctree)
   end if

   ! Reorder linear operator
   call linop_reorder(ndata%h(il0i))

   ! Release memory
   deallocate(missing)
end do

! Release memory
call linop_dealloc(hbase)
deallocate(valid)
deallocate(mask_extra)

end subroutine compute_interp_h

!----------------------------------------------------------------------
! Subroutine: compute_interp_v
!> Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine compute_interp_v(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: il0,jl0,il1,il0inf,il0sup

! Initialize vertical interpolation
ndata%v%prefix = 'v'
ndata%v%n_src = ndata%nl1
ndata%v%n_dst = ndata%nl0

! Linear interpolation
ndata%v%n_s = ndata%nl1
il0inf = 1
do il0=1,ndata%nl0
   if (ndata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         ndata%v%n_s = ndata%v%n_s+2
      end do
      il0inf = il0
   end if
end do
call linop_alloc(ndata%v)
do il1=1,ndata%nl1
   il0 = ndata%il1_to_il0(il1)
   ndata%v%row(il1) = il0
   ndata%v%col(il1) = il0
   ndata%v%S(il1) = 1.0
end do
ndata%v%n_s = ndata%nl1
il0inf = 1
do il0=1,ndata%nl0
   if (ndata%llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         ndata%v%n_s = ndata%v%n_s+1
         ndata%v%row(ndata%v%n_s) = jl0
         ndata%v%col(ndata%v%n_s) = il0inf
         ndata%v%S(ndata%v%n_s) = abs(ndata%geom%vunit(il0sup)-ndata%geom%vunit(jl0)) & 
                            & /abs(ndata%geom%vunit(il0sup)-ndata%geom%vunit(il0inf))

         ndata%v%n_s = ndata%v%n_s+1
         ndata%v%row(ndata%v%n_s) = jl0
         ndata%v%col(ndata%v%n_s) = il0sup
         ndata%v%S(ndata%v%n_s) = abs(ndata%geom%vunit(jl0)-ndata%geom%vunit(il0inf)) &
                            & /abs(ndata%geom%vunit(il0sup)-ndata%geom%vunit(il0inf))
      end do
      il0inf = il0
   end if
end do

! Conversion
ndata%v%col = ndata%il0_to_il1(ndata%v%col)

! Reorder linear operator
call linop_reorder(ndata%v)

end subroutine compute_interp_v

!----------------------------------------------------------------------
! Subroutine: compute_interp_s
!> Purpose: compute horizontal subsampling interpolation
!----------------------------------------------------------------------
subroutine compute_interp_s(nam,ndata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(inout) :: ndata !< Sampling data

! Local variables
integer :: ic1,ic2,il1,i_s
integer,allocatable :: mask_ctree(:)
real(kind_real) :: dum(1)
real(kind_real) :: renorm(ndata%nc1)
logical,allocatable :: valid(:),missing(:)
type(linoptype) :: stmp
type(ctreetype) :: ctree

! Allocation
allocate(ndata%s(ndata%nl1))

do il1=1,ndata%nl1
   ! Initialize object
   ndata%s(il1)%prefix = 's'
   ndata%s(il1)%n_src = ndata%nc2(il1)
   ndata%s(il1)%n_dst = ndata%nc1

   if (ndata%nc2(il1)==ndata%nc1) then
      ! No interpolation
      ndata%s(il1)%n_s = ndata%nc2(il1)
      call linop_alloc(ndata%s(il1))
      do i_s=1,ndata%s(il1)%n_s
         ndata%s(il1)%row(i_s) = i_s
         ndata%s(il1)%col(i_s) = i_s
         ndata%s(il1)%S(i_s) = 1.0
      end do
   else
      ! Compute interpolation
      call interp_horiz(rng,ndata%nc2(il1),ndata%geom%lon(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1)), & 
    & ndata%geom%lat(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1)), &
    & ndata%geom%mask(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1),ndata%il1_to_il0(il1)), &
    & ndata%nc1,ndata%geom%lon(ndata%ic1_to_ic0),ndata%geom%lat(ndata%ic1_to_ic0), &
    & ndata%geom%mask(ndata%ic1_to_ic0,ndata%il1_to_il0(il1)),stmp)
   
      ! Allocation
      allocate(valid(stmp%n_s))
      valid = .true.
   
      ! Check mask boundaries
      if (nam%mask_check) then
         call check_mask_bnd(nam,ndata,stmp,valid,ndata%il1_to_il0(il1), &
       & row_to_ic0=ndata%ic1_to_ic0,col_to_ic0=ndata%ic1_to_ic0(ndata%ic2il1_to_ic1(1:ndata%nc2(il1),il1)))
      else
         write(mpl%unit,'(a10,a,i3)') '','Level ',nam%levs(ndata%il1_to_il0(il1))
      end if
   
      ! Renormalization
      renorm = 0.0
      do i_s=1,stmp%n_s
         if (valid(i_s)) renorm(stmp%row(i_s)) = renorm(stmp%row(i_s))+stmp%S(i_s)
      end do
   
      ! Copy valid operations
      ndata%s(il1)%n_s = count(valid)
      call linop_alloc(ndata%s(il1))
      ndata%s(il1)%n_s = 0
      do i_s=1,stmp%n_s
         if (valid(i_s)) then
            ndata%s(il1)%n_s = ndata%s(il1)%n_s+1
            ndata%s(il1)%row(ndata%s(il1)%n_s) = stmp%row(i_s)
            ndata%s(il1)%col(ndata%s(il1)%n_s) = stmp%col(i_s)
            ndata%s(il1)%S(ndata%s(il1)%n_s) = stmp%S(i_s)/renorm(stmp%row(i_s))
         end if
      end do
   
      ! Release memory
      call linop_dealloc(stmp)
      deallocate(valid)
   
      ! Allocation
      allocate(missing(ndata%nc1))
   
      ! Count points that are not interpolated
      missing = .false.
      do ic1=1,ndata%nc1
         if (ndata%geom%mask(ndata%ic1_to_ic0(ic1),ndata%il1_to_il0(il1))) missing(ic1) = .true.
      end do
      do i_s=1,ndata%s(il1)%n_s
         missing(ndata%s(il1)%row(i_s)) = .false.
      end do
      if (count(missing)>0) then
         ! Copy
         call linop_copy(ndata%s(il1),stmp)
   
         ! Reallocate permanent arrays
         call linop_dealloc(ndata%s(il1))
         ndata%s(il1)%n_s = ndata%s(il1)%n_s+count(missing)
         call linop_alloc(ndata%s(il1))
   
         ! Fill permanent arrays
         ndata%s(il1)%row(1:stmp%n_s) = stmp%row
         ndata%s(il1)%col(1:stmp%n_s) = stmp%col
         ndata%s(il1)%S(1:stmp%n_s) = stmp%S
   
         ! Compute cover tree
         allocate(mask_ctree(ndata%nc2(il1)))
         mask_ctree = 0
         do ic2=1,ndata%nc2(il1)
            if (ndata%geom%mask(ndata%ic1_to_ic0(ndata%ic2il1_to_ic1(ic2,il1)),ndata%il1_to_il0(il1))) mask_ctree(ic2) = 1
         end do
         ctree = create_ctree(ndata%nc2(il1),dble(ndata%geom%lon(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1))), &
       & dble(ndata%geom%lat(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1))),mask_ctree)
         deallocate(mask_ctree)
   
         ! Compute nearest neighbors
         do ic1=1,ndata%nc1
            if (missing(ic1)) then
               stmp%n_s = stmp%n_s+1
               ndata%s(il1)%row(stmp%n_s) = ic1
               call find_nearest_neighbors(ctree,dble(ndata%geom%lon(ndata%ic1_to_ic0(ic1))), &
             & dble(ndata%geom%lat(ndata%ic1_to_ic0(ic1))),1,ndata%s(il1)%col(stmp%n_s:stmp%n_s),dum)
               ndata%s(il1)%S(stmp%n_s) = 1.0
            end if
         end do
   
         ! Release memory
         call linop_dealloc(stmp)
         call delete_ctree(ctree)
      end if
   
      ! Reorder linear operator
      call linop_reorder(ndata%s(il1))
   
      ! Release memory
      deallocate(missing)
   end if
end do

end subroutine compute_interp_s

!----------------------------------------------------------------------
! Subroutine: check_mask_bnd
!> Purpose: check mask boundaries for interpolations
!----------------------------------------------------------------------
subroutine check_mask_bnd(nam,ndata,interp,valid,il0,row_to_ic0,col_to_ic0)

implicit none

! Passed variables
type(namtype),intent(in) :: nam !< Namelist variables
type(ndatatype),intent(in) :: ndata          !< Sampling data
type(linoptype),intent(inout) :: interp          !< Interpolation
logical,intent(inout) :: valid(interp%n_s)   !< Valid vector
integer,intent(in),optional :: row_to_ic0(interp%n_dst) !< Conversion from row to ic0 (identity if missing)
integer,intent(in),optional :: col_to_ic0(interp%n_src) !< Conversion from col to ic0 (identity if missing)

! Local variables
integer :: ic0,i_s,ibnd,jc0,jc1,progint,il0
integer :: iproc,i_s_s(mpl%nproc),i_s_e(mpl%nproc),n_s_loc(mpl%nproc),i_s_loc
real(kind_real),allocatable :: x(:),y(:),z(:),v1(:),v2(:),va(:),vp(:),t(:)
logical,allocatable :: done(:)

! MPI splitting
do iproc=1,mpl%nproc
   i_s_s(iproc) = (iproc-1)*(interp%n_s/mpl%nproc+1)+1
   i_s_e(iproc) = min(iproc*(interp%n_s/mpl%nproc+1),interp%n_s)
   n_s_loc(iproc) = i_s_e(iproc)-i_s_s(iproc)+1
end do

! Allocation
allocate(done(n_s_loc(mpl%myproc)))

! Check that interpolations are not crossing mask boundaries
write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),': '
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
      if (present(row_to_ic0)) then
         ic0 = row_to_ic0(interp%row(i_s))
      else
         ic0 = interp%row(i_s)
      end if
      if (present(col_to_ic0)) then
         jc0 = col_to_ic0(interp%col(i_s))
      else
         jc0 = interp%col(i_s)
      end if

      ! Transform to cartesian coordinates
      call trans(2,ndata%geom%lat((/ic0,jc0/)),ndata%geom%lon((/ic0,jc0/)),x,y,z)

      ! Compute arc orthogonal vector
      v1 = (/x(1),y(1),z(1)/)
      v2 = (/x(2),y(2),z(2)/)
      call vector_product(v1,v2,va)

      ! Check if arc is crossing boundary arcs
      do ibnd=1,ndata%geom%nbnd(il0)
         call vector_product(va,ndata%geom%vbnd(:,ibnd,il0),vp)
         v1 = (/x(1),y(1),z(1)/)
         call vector_triple_product(v1,va,vp,t(1))
         v1 = (/x(2),y(2),z(2)/)
         call vector_triple_product(v1,va,vp,t(2))
         v1 = (/ndata%geom%xbnd(1,ibnd,il0),ndata%geom%ybnd(1,ibnd,il0),ndata%geom%zbnd(1,ibnd,il0)/)
         call vector_triple_product(v1,ndata%geom%vbnd(:,ibnd,il0),vp,t(3))
         v1 = (/ndata%geom%xbnd(2,ibnd,il0),ndata%geom%ybnd(2,ibnd,il0),ndata%geom%zbnd(2,ibnd,il0)/)
         call vector_triple_product(v1,ndata%geom%vbnd(:,ibnd,il0),vp,t(4))
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

! Broadcast
do iproc=1,mpl%nproc
   call mpl_bcast(valid(i_s_s(iproc):i_s_e(iproc)),iproc)
end do

! Release memory
deallocate(done)

end subroutine check_mask_bnd

end module module_parameters_interp
