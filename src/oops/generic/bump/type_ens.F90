!----------------------------------------------------------------------
! Module: type_ens
! Purpose: ensemble derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_ens

use tools_kinds, only: kind_real
use tools_missing, only: msi,msr
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

! Ensemble derived type
type ens_type
   ! Attributes
   integer :: ne                                  ! Ensemble size
   integer :: nsub                                ! Number of sub-ensembles

   ! Data
   real(kind_real),allocatable :: fld(:,:,:,:,:)  ! Ensemble perturbation
   real(kind_real),allocatable :: mean(:,:,:,:,:) ! Ensemble mean
contains
   procedure :: alloc => ens_alloc
   procedure :: dealloc => ens_dealloc
   procedure :: copy => ens_copy
   procedure :: remove_mean => ens_remove_mean
   procedure :: ens_from
   procedure :: ens_from_nemovar
   generic :: from => ens_from,ens_from_nemovar
end type ens_type

private
public :: ens_type

contains

!----------------------------------------------------------------------
! Subroutine: ens_alloc
! Purpose: ensemble data allocation
!----------------------------------------------------------------------
subroutine ens_alloc(ens,nam,geom,ne,nsub)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble
type(nam_type),intent(in) :: nam     ! Namelist
type(geom_type),intent(in) :: geom   ! Geometry
integer,intent(in) :: ne             ! Ensemble size
integer,intent(in) :: nsub           ! Number of sub-ensembles

! Allocate
if (ne>0) then
   allocate(ens%fld(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne))
   allocate(ens%mean(geom%nc0a,geom%nl0,nam%nv,nam%nts,nsub))
end if

! Initialization
ens%ne = ne
ens%nsub = nsub
if (ens%ne>0) then
   call msr(ens%fld)
   call msr(ens%mean)
end if

end subroutine ens_alloc

!----------------------------------------------------------------------
! Subroutine: ens_dealloc
! Purpose: ensemble data deallocation
!----------------------------------------------------------------------
subroutine ens_dealloc(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble

! Release memory
if (allocated(ens%fld)) deallocate(ens%fld)
if (allocated(ens%mean)) deallocate(ens%mean)

end subroutine ens_dealloc

!----------------------------------------------------------------------
! Subroutine: ens_copy
! Purpose: ensemble data copy
!----------------------------------------------------------------------
subroutine ens_copy(ens_out,ens_in)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens_out ! Ensemble
type(ens_type),intent(in) :: ens_in      ! Ensemble

! Allocate
if (ens_in%ne>0) then
   allocate(ens_out%fld(size(ens_in%fld,1),size(ens_in%fld,2),size(ens_in%fld,3),size(ens_in%fld,4),size(ens_in%fld,5)))
   allocate(ens_out%mean(size(ens_in%mean,1),size(ens_in%mean,2),size(ens_in%mean,3),size(ens_in%mean,4),size(ens_in%mean,5)))
end if

! Initialization
ens_out%ne = ens_in%ne
ens_out%nsub = ens_in%nsub
if (ens_in%ne>0) then
   ens_out%fld = ens_in%fld
   ens_out%mean = ens_in%mean
end if

end subroutine ens_copy

!----------------------------------------------------------------------
! Subroutine: ens_remove_mean
! Purpose: remove ensemble mean
!----------------------------------------------------------------------
subroutine ens_remove_mean(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens ! Ensemble

! Local variables
integer :: isub,ie_sub,ie

if (ens%ne>0) then
   ! Loop over sub-ensembles
   do isub=1,ens%nsub
      ! Compute mean
      ens%mean(:,:,:,:,isub) = 0.0
      do ie_sub=1,ens%ne/ens%nsub
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub
         ens%mean(:,:,:,:,isub) = ens%mean(:,:,:,:,isub)+ens%fld(:,:,:,:,ie)
      end do
      ens%mean(:,:,:,:,isub) = ens%mean(:,:,:,:,isub)/(ens%ne/ens%nsub)

      ! Remove mean
      do ie_sub=1,ens%ne/ens%nsub
         ie = ie_sub+(isub-1)*ens%ne/ens%nsub
         ens%fld(:,:,:,:,ie) = ens%fld(:,:,:,:,ie)-ens%mean(:,:,:,:,isub)
      end do
   end do
end if

end subroutine ens_remove_mean

!----------------------------------------------------------------------
! Subroutine: ens_from
! Purpose: copy ensemble array into ensemble data
!----------------------------------------------------------------------
subroutine ens_from(ens,nam,geom,ne,ens_mga)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                        ! Ensemble
type(nam_type),intent(in) :: nam                                            ! Namelist
type(geom_type),intent(in) :: geom                                          ! Geometry
integer,intent(in) :: ne                                                    ! Ensemble size
real(kind_real),intent(in) :: ens_mga(geom%nmga,geom%nl0,nam%nv,nam%nts,ne) ! Ensemble on model grid, halo A

! Local variables
integer :: ie,its,iv,il0

! Allocation
call ens%alloc(nam,geom,ne,1)

if (ens%ne>0) then
   ! Copy
   do ie=1,ens%ne
      do its=1,nam%nts
         do iv=1,nam%nv
            do il0=1,geom%nl0
               ens%fld(:,il0,iv,its,ie) = ens_mga(geom%c0a_to_mga,il0,iv,its,ie)
            end do
         end do
      end do
   end do

   ! Remove mean
   call ens%remove_mean
end if

end subroutine ens_from

!----------------------------------------------------------------------
! Subroutine: ens_from_nemovar
! Purpose: copy 2d NEMOVAR ensemble into ensemble data
!----------------------------------------------------------------------
subroutine ens_from_nemovar(ens,mpl,nam,geom,nx,ny,nens,ncyc,ens_2d,ens_3d)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                    ! Ensemble
type(mpl_type),intent(in) :: mpl                                        ! MPI data
type(nam_type),intent(in) :: nam                                        ! Namelist
type(geom_type),intent(in) :: geom                                      ! Geometry
integer,intent(in) :: nx                                                ! X-axis size
integer,intent(in) :: ny                                                ! Y-axis size
integer,intent(in) :: nens                                              ! Ensemble size at each cycle
integer,intent(in) :: ncyc                                              ! Number of cycles
real(kind_real),intent(in),optional :: ens_2d(nx,ny,nens,ncyc)          ! Ensemble on model grid, halo A
real(kind_real),intent(in),optional :: ens_3d(nx,ny,geom%nl0,nens,ncyc) ! Ensemble on model grid, halo A

! Local variables
integer :: ie,iens,icyc,its,iv,il0
real(kind_real) :: tmp(geom%nmga)

! Allocation
call ens%alloc(nam,geom,nens*ncyc,ncyc)

! Copy
iv = 1
its = 1
ie = 1
do iens=1,nens
   do icyc=1,ncyc
      ! Pack and copy
      if (present(ens_2d)) then
         il0 = 1
         tmp = pack(ens_2d(:,:,iens,icyc),.true.)
         ens%fld(:,il0,iv,its,ie) = tmp(geom%c0a_to_mga)
      elseif (present(ens_3d)) then
         do il0=1,geom%nl0
            ! Pack and copy
            tmp = pack(ens_3d(:,:,il0,iens,icyc),.true.)
            ens%fld(:,il0,iv,its,ie) = tmp(geom%c0a_to_mga)
         end do
      else
         call mpl%abort('ens_2d or ens_3d should be provided in ens_from_nemovar')
      end if

      ! Update
      ie = ie+1
   end do
end do

! Remove mean
call ens%remove_mean

end subroutine ens_from_nemovar

end module type_ens
