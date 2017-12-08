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

use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist
use tools_display, only: msgerror,prog_init,prog_print
use tools_interp, only: compute_grid_interp_bilin,compute_interp_bilin,check_mask_bnd
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_stripack, only: trans
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy,linop_reorder
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_nam, only: namtype
use type_ndata, only: ndatatype
use type_randgen, only: initialize_sampling

implicit none

private
public :: compute_interp_h,compute_interp_v,compute_interp_s

contains

!----------------------------------------------------------------------
! Subroutine: compute_interp_h
!> Purpose: compute basic horizontal interpolation
!----------------------------------------------------------------------
subroutine compute_interp_h(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Allocation
allocate(ndata%h(geom%nl0i))

! Compute grid interpolation
call compute_grid_interp_bilin(geom,ndata%nc1,ndata%ic1_to_ic0,nam%mask_check,ndata%vbot,ndata%vtop,ndata%h)

! End associate
end associate

end subroutine compute_interp_h

!----------------------------------------------------------------------
! Subroutine: compute_interp_v
!> Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine compute_interp_v(ndata,llev)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata     !< NICAS data
logical,intent(in) :: llev(ndata%geom%nl0) !< Levels selection

! Local variables
integer :: il0,jl0,il1,il0inf,il0sup

! Associate
associate(geom=>ndata%geom)

! Initialize vertical interpolation
ndata%v%prefix = 'v'
ndata%v%n_src = ndata%nl1
ndata%v%n_dst = geom%nl0

! Linear interpolation
ndata%v%n_s = ndata%nl1
il0inf = 1
do il0=1,geom%nl0
   if (llev(il0)) then
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
do il0=1,geom%nl0
   if (llev(il0)) then
      il0sup = il0
      do jl0=il0inf+1,il0sup-1
         ndata%v%n_s = ndata%v%n_s+1
         ndata%v%row(ndata%v%n_s) = jl0
         ndata%v%col(ndata%v%n_s) = il0inf
         ndata%v%S(ndata%v%n_s) = abs(geom%vunit(il0sup)-geom%vunit(jl0)) &
                            & /abs(geom%vunit(il0sup)-geom%vunit(il0inf))

         ndata%v%n_s = ndata%v%n_s+1
         ndata%v%row(ndata%v%n_s) = jl0
         ndata%v%col(ndata%v%n_s) = il0sup
         ndata%v%S(ndata%v%n_s) = abs(geom%vunit(jl0)-geom%vunit(il0inf)) &
                            & /abs(geom%vunit(il0sup)-geom%vunit(il0inf))
      end do
      il0inf = il0
   end if
end do

! Conversion
ndata%v%col = ndata%il0_to_il1(ndata%v%col)

! Reorder linear operator
call linop_reorder(ndata%v)

! End associate
end associate

end subroutine compute_interp_v

!----------------------------------------------------------------------
! Subroutine: compute_interp_s
!> Purpose: compute horizontal subsampling interpolation
!----------------------------------------------------------------------
subroutine compute_interp_s(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: ic1,il1,i_s
real(kind_real) :: dum(1)
real(kind_real) :: renorm(ndata%nc1)
logical,allocatable :: valid(:),missing(:)
type(linoptype) :: stmp
type(ctreetype) :: ctree

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

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
      call compute_interp_bilin(ndata%nc2(il1),geom%lon(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1)), &
    & geom%lat(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1)), &
    & geom%mask(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1),ndata%il1_to_il0(il1)), &
    & ndata%nc1,geom%lon(ndata%ic1_to_ic0),geom%lat(ndata%ic1_to_ic0), &
    & geom%mask(ndata%ic1_to_ic0,ndata%il1_to_il0(il1)),stmp)

      ! Allocation
      allocate(valid(stmp%n_s))
      valid = .true.

      ! Check mask boundaries
      if (nam%mask_check) then
         write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Sublevel ',il1,': '
         call check_mask_bnd(geom,stmp,valid,ndata%il1_to_il0(il1), &
       & row_to_ic0=ndata%ic1_to_ic0,col_to_ic0=ndata%ic1_to_ic0(ndata%ic2il1_to_ic1(1:ndata%nc2(il1),il1)))
      else
         write(mpl%unit,'(a10,a,i3)') '','Sublevel ',il1
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
         if (geom%mask(ndata%ic1_to_ic0(ic1),ndata%il1_to_il0(il1))) missing(ic1) = .true.
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
         ctree = create_ctree(ndata%nc2(il1),dble(geom%lon(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1))), &
       & dble(geom%lat(ndata%ic2il1_to_ic0(1:ndata%nc2(il1),il1))), &
       & geom%mask(ndata%ic1_to_ic0(ndata%ic2il1_to_ic1(1:ndata%nc2(il1),il1)),ndata%il1_to_il0(il1)))

         ! Compute nearest neighbors
         do ic1=1,ndata%nc1
            if (missing(ic1)) then
               stmp%n_s = stmp%n_s+1
               ndata%s(il1)%row(stmp%n_s) = ic1
               call find_nearest_neighbors(ctree,dble(geom%lon(ndata%ic1_to_ic0(ic1))), &
             & dble(geom%lat(ndata%ic1_to_ic0(ic1))),1,ndata%s(il1)%col(stmp%n_s:stmp%n_s),dum)
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

! End associate
end associate

end subroutine compute_interp_s

end module module_parameters_interp
