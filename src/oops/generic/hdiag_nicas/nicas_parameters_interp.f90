!----------------------------------------------------------------------
! Module: nicas_parameters_interp.f90
!> Purpose: compute NICAS parameters (interpolation)
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_parameters_interp

use omp_lib
use tools_const, only: pi,req,deg2rad,rad2deg,sphere_dist
use tools_display, only: msgerror,prog_init,prog_print
use tools_interp, only: compute_grid_interp,compute_interp,check_mask_bnd,interp_missing
use tools_kinds,only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr,isnotmsi
use tools_stripack, only: trans
use type_ctree, only: ctreetype,create_ctree,find_nearest_neighbors,delete_ctree
use type_linop, only: linoptype,linop_alloc,linop_dealloc,linop_copy
use type_mpl, only: mpl,mpl_bcast,mpl_recv,mpl_send
use type_nam, only: namtype
use type_ndata, only: ndatatype

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
allocate(ndata%hfull(geom%nl0i))

! Compute grid interpolation
call compute_grid_interp(geom,ndata%nc1,ndata%c1_to_c0,nam%mask_check,ndata%vbot,ndata%vtop,nam%nicas_interp,ndata%hfull)

! End associate
end associate

end subroutine compute_interp_h

!----------------------------------------------------------------------
! Subroutine: compute_interp_v
!> Purpose: compute vertical interpolation
!----------------------------------------------------------------------
subroutine compute_interp_v(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: jl0,il0,jl1,il0inf,il0sup

! Associate
associate(geom=>ndata%geom)

! Initialize vertical interpolation
ndata%v%prefix = 'v'
ndata%v%n_src = ndata%nl1
ndata%v%n_dst = geom%nl0

! Linear interpolation
ndata%v%n_s = ndata%nl1
il0inf = 1
do jl0=1,geom%nl0
   if (ndata%llev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         ndata%v%n_s = ndata%v%n_s+2
      end do
      il0inf = jl0
   end if
end do
call linop_alloc(ndata%v)
do jl1=1,ndata%nl1
   jl0 = ndata%l1_to_l0(jl1)
   ndata%v%row(jl1) = jl0
   ndata%v%col(jl1) = jl0
   ndata%v%S(jl1) = 1.0
end do
ndata%v%n_s = ndata%nl1
il0inf = 1
do jl0=1,geom%nl0
   if (ndata%llev(jl0)) then
      il0sup = jl0
      do il0=il0inf+1,il0sup-1
         ndata%v%n_s = ndata%v%n_s+1
         ndata%v%row(ndata%v%n_s) = il0
         ndata%v%col(ndata%v%n_s) = il0inf
         ndata%v%S(ndata%v%n_s) = abs(geom%vunit(il0sup)-geom%vunit(il0)) &
                            & /abs(geom%vunit(il0sup)-geom%vunit(il0inf))

         ndata%v%n_s = ndata%v%n_s+1
         ndata%v%row(ndata%v%n_s) = il0
         ndata%v%col(ndata%v%n_s) = il0sup
         ndata%v%S(ndata%v%n_s) = abs(geom%vunit(il0)-geom%vunit(il0inf)) &
                            & /abs(geom%vunit(il0sup)-geom%vunit(il0inf))
      end do
      il0inf = jl0
   end if
end do

! Conversion
ndata%v%col = ndata%l0_to_l1(ndata%v%col)

! Deallocate selected levels
deallocate(ndata%llev)

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
integer :: il1,i_s
real(kind_real) :: renorm(ndata%nc1)
logical :: mask_src(ndata%nc1),mask_dst(ndata%nc1)
logical,allocatable :: valid(:)
type(linoptype) :: stmp

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Allocation
allocate(ndata%sfull(ndata%nl1))

do il1=1,ndata%nl1
   ! Initialize object
   ndata%sfull(il1)%prefix = 's'
   ndata%sfull(il1)%n_src = ndata%nc1
   ndata%sfull(il1)%n_dst = ndata%nc1

   if (ndata%nc2(il1)==ndata%nc1) then
      ! No interpolation
      ndata%sfull(il1)%n_s = ndata%nc1
      call linop_alloc(ndata%sfull(il1))
      do i_s=1,ndata%sfull(il1)%n_s
         ndata%sfull(il1)%row(i_s) = i_s
         ndata%sfull(il1)%col(i_s) = i_s
         ndata%sfull(il1)%S(i_s) = 1.0
      end do
   else
      ! Mask
      mask_src = ndata%c2mask(:,il1)
      mask_dst = .true.

      ! Compute interpolation
      call compute_interp(ndata%nc1,geom%lon(ndata%c1_to_c0),geom%lat(ndata%c1_to_c0),mask_src, &
    & ndata%nc1,geom%lon(ndata%c1_to_c0),geom%lat(ndata%c1_to_c0),mask_dst,nam%nicas_interp,stmp)

      ! Allocation
      allocate(valid(stmp%n_s))
      valid = .true.

      ! Check mask boundaries
      if (nam%mask_check) then
         write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Sublevel ',il1,': '
         call check_mask_bnd(geom,stmp,valid,ndata%l1_to_l0(il1), &
       & row_to_ic0=ndata%c1_to_c0,col_to_ic0=ndata%c1_to_c0)
      else
         write(mpl%unit,'(a10,a,i3)') '','Sublevel ',il1
      end if

      ! Renormalization
      renorm = 0.0
      do i_s=1,stmp%n_s
         if (valid(i_s)) renorm(stmp%row(i_s)) = renorm(stmp%row(i_s))+stmp%S(i_s)
      end do

      ! Copy valid operations
      ndata%sfull(il1)%n_s = count(valid)
      call linop_alloc(ndata%sfull(il1))
      ndata%sfull(il1)%n_s = 0
      do i_s=1,stmp%n_s
         if (valid(i_s)) then
            ndata%sfull(il1)%n_s = ndata%sfull(il1)%n_s+1
            ndata%sfull(il1)%row(ndata%sfull(il1)%n_s) = stmp%row(i_s)
            ndata%sfull(il1)%col(ndata%sfull(il1)%n_s) = stmp%col(i_s)
            ndata%sfull(il1)%S(ndata%sfull(il1)%n_s) = stmp%S(i_s)/renorm(stmp%row(i_s))
         end if
      end do

      ! Release memory
      call linop_dealloc(stmp)
      deallocate(valid)

      ! Count points that are not interpolated
      call interp_missing(ndata%nc1,geom%lon(ndata%c1_to_c0),geom%lat(ndata%c1_to_c0),mask_dst,nam%nicas_interp,ndata%sfull(il1))
   end if
end do

! End associate
end associate

end subroutine compute_interp_s

end module nicas_parameters_interp
