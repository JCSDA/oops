!----------------------------------------------------------------------
! Module: type_ens
!> Purpose: ensemble derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
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
   integer :: ne                                 !< Ensemble size
   integer :: nsub                               !< Number of sub-ensembles

   ! Data
   real(kind_real),allocatable :: fld(:,:,:,:,:) !< Ensemble fields
contains
   procedure :: alloc => ens_alloc
   procedure :: dealloc => ens_dealloc
   procedure :: ens_from
   procedure :: ens_from_oops
   procedure :: ens_from_nemovar
   generic :: from => ens_from,ens_from_oops,ens_from_nemovar
end type ens_type

private
public :: ens_type

contains

!----------------------------------------------------------------------
! Subroutine: ens_alloc
!> Purpose: ensemble data allocation
!----------------------------------------------------------------------
subroutine ens_alloc(ens,nam,geom,ne,nsub)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens !< Ensemble
type(nam_type),intent(in) :: nam     !< Namelist
type(geom_type),intent(in) :: geom   !< Geometry
integer,intent(in) :: ne             !< Ensemble size
integer,intent(in) :: nsub           !< Number of sub-ensembles

! Allocate
allocate(ens%fld(geom%nc0a,geom%nl0,nam%nv,nam%nts,ne))

! Initialization
ens%ne = ne
ens%nsub = nsub
call msr(ens%fld)

end subroutine ens_alloc

!----------------------------------------------------------------------
! Subroutine: ens_dealloc
!> Purpose: ensemble data deallocation
!----------------------------------------------------------------------
subroutine ens_dealloc(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens !< Ensemble

! Release memory
if (allocated(ens%fld)) deallocate(ens%fld)

end subroutine ens_dealloc

!----------------------------------------------------------------------
! Subroutine: ens_from
!> Purpose: copy ensemble array into ensemble data
!----------------------------------------------------------------------
subroutine ens_from(ens,nam,geom,ne,ens_mga)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                        !< Ensemble
type(nam_type),intent(in) :: nam                                            !< Namelist
type(geom_type),intent(in) :: geom                                          !< Geometry
integer,intent(in) :: ne                                                    !< Ensemble size
real(kind_real),intent(in) :: ens_mga(geom%nmga,geom%nl0,nam%nv,nam%nts,ne) !< Ensemble on model grid, halo A

! Local variables
integer :: ie,its,iv,il0

! Allocation
call ens%alloc(nam,geom,ne,1)

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

end subroutine ens_from

!----------------------------------------------------------------------
! Subroutine: ens_from_oops
!> Purpose: copy OOPS ensemble into ensemble data
!----------------------------------------------------------------------
subroutine ens_from_oops(ens,nam,geom,ne,ens_mga)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                        !< Ensemble
type(nam_type),intent(in) :: nam                                            !< Namelist
type(geom_type),intent(in) :: geom                                          !< Geometry
integer,intent(in) :: ne                                                    !< Ensemble size
real(kind_real),intent(in) :: ens_mga(geom%nmga*geom%nl0*nam%nv*nam%nts*ne) !< Ensemble on model grid, halo A

! Local variables
integer :: ie,its,iv,il0,offset,n
real(kind_real) :: fld(geom%nmga,geom%nl0,nam%nv,nam%nts)
logical :: mask_unpack(geom%nmga,geom%nl0,nam%nv,nam%nts)

! Allocation
call ens%alloc(nam,geom,ne,1)

! Initialization
mask_unpack = .true.
offset = 0
n = geom%nmga*geom%nl0*nam%nv*nam%nts

do ie=1,ne
   ! Unpack
   call msr(fld)
   fld = unpack(ens_mga(offset+1:offset+n),mask_unpack,fld)

   ! Copy
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            ens%fld(:,il0,iv,its,ie) = fld(geom%c0a_to_mga,il0,iv,its)
         end do
      end do
   end do

   ! Update
   offset = offset+n
end do

end subroutine ens_from_oops

!----------------------------------------------------------------------
! Subroutine: ens_from_nemovar
!> Purpose: copy 2d NEMOVAR ensemble into ensemble data
!----------------------------------------------------------------------
subroutine ens_from_nemovar(ens,mpl,nam,geom,nx,ny,nens,ncyc,ens_2d,ens_3d)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                    !< Ensemble
type(mpl_type),intent(in) :: mpl                                        !< MPI data
type(nam_type),intent(in) :: nam                                        !< Namelist
type(geom_type),intent(in) :: geom                                      !< Geometry
integer,intent(in) :: nx                                                !< X-axis size
integer,intent(in) :: ny                                                !< Y-axis size
integer,intent(in) :: nens                                              !< Ensemble size at each cycle
integer,intent(in) :: ncyc                                              !< Number of cycles
real(kind_real),intent(in),optional :: ens_2d(nx,ny,nens,ncyc)          !< Ensemble on model grid, halo A
real(kind_real),intent(in),optional :: ens_3d(nx,ny,geom%nl0,nens,ncyc) !< Ensemble on model grid, halo A

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

end subroutine ens_from_nemovar

end module type_ens
