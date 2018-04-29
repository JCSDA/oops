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

use model_interface, only: model_read
use omp_lib
use tools_display, only: msgwarning,msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr
use type_geom, only: geom_type
use type_mpl, only: mpl
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
   procedure :: load => ens_load
   procedure :: from_xyz => ens_from_xyz
   procedure :: from_xz => ens_from_xz
end type ens_type

private
public :: ens_type

contains

!----------------------------------------------------------------------
! Subroutine: ens_alloc
!> Purpose: ensemble allocation
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
!> Purpose: ensemble deallocation
!----------------------------------------------------------------------
subroutine ens_dealloc(ens)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens !< Ensemble

! Allocate
deallocate(ens%fld)

! Initialization
call msi(ens%ne)
call msi(ens%nsub)

end subroutine ens_dealloc

!----------------------------------------------------------------------
! Subroutine: ens_load
!> Purpose: load ensemble
!----------------------------------------------------------------------
subroutine ens_load(ens,nam,geom,filename)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens    !< Ensemble
type(nam_type),intent(in) :: nam        !< Namelist
type(geom_type),intent(in) :: geom      !< Geometry
character(len=*),intent(in) :: filename !< Filename ('ens1' or 'ens2')

! Local variables
integer :: ne,ne_offset,nsub,isub,jsub,ie,ietot

! Setup
select case (trim(filename))
case ('ens1')
   ne = nam%ens1_ne
   ne_offset = nam%ens1_ne_offset
   nsub = nam%ens1_nsub
case ('ens2')
   ne = nam%ens2_ne
   ne_offset = nam%ens2_ne_offset
   nsub = nam%ens2_nsub
case default
   call msi(ne)
   call msi(ne_offset)
   call msi(nsub)
   call msgerror('wrong filename in ens_load')
end select

! Allocation
call ens%alloc(nam,geom,ne,nsub)

! Initialization
call msr(ens%fld)
ietot = 1

! Loop over sub-ensembles
do isub=1,ens%nsub
   if (ens%nsub==1) then
      write(mpl%unit,'(a7,a)',advance='no') '','Full ensemble, member:'
   else
      write(mpl%unit,'(a7,a,i4,a)',advance='no') '','Sub-ensemble ',isub,', member:'
   end if
   call flush(mpl%unit)

   ! Loop over members for a given sub-ensemble
   do ie=1,ens%ne/ens%nsub
      write(mpl%unit,'(i4)',advance='no') ne_offset+ie
      call flush(mpl%unit)

      ! Read member
      if (ens%nsub==1) then
         jsub = 0
      else
         jsub = isub
      end if
      call model_read(nam,geom,trim(filename),ne_offset+ie,jsub,ens%fld(:,:,:,:,ietot))

      ! Update
      ietot = ietot+1
   end do
   write(mpl%unit,'(a)') ''
   call flush(mpl%unit)
end do

end subroutine ens_load

!----------------------------------------------------------------------
! Subroutine: ens_from_xyz
!> Purpose: copy xyz ensemble into ensemble object
!----------------------------------------------------------------------
subroutine ens_from_xyz(ens,nam,geom,nx,ny,ne,ens_ga)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                   !< Ensemble
type(nam_type),intent(in) :: nam                                       !< Namelist
type(geom_type),intent(in) :: geom                                     !< Geometry
integer,intent(in) :: nx                                               !< X-axis size
integer,intent(in) :: ny                                               !< Y-axis size
integer,intent(in) :: ne                                               !< Ensemble size
real(kind_real),intent(in) :: ens_ga(nx,ny,geom%nl0,nam%nv,nam%nts,ne) !< Ensemble on model grid, halo A

! Local variables
integer :: ie,its,iv,il0
real(kind_real),allocatable :: tmp(:)

! Allocation
call ens%alloc(nam,geom,ne,1)

! Copy
!$omp parallel do schedule(static) private(ie,its,iv,il0) firstprivate(tmp)
do ie=1,ens%ne
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            ! Allocation
            allocate(tmp(geom%nmga))

            ! Pack and copy
            tmp = pack(ens_ga(:,:,il0,iv,its,ie),.true.)
            ens%fld(:,il0,iv,its,ie) = tmp(geom%c0a_to_ga)
         end do
      end do
   end do
end do
!$omp end parallel do

end subroutine ens_from_xyz

!----------------------------------------------------------------------
! Subroutine: ens_from_xz
!> Purpose: copy xz ensemble into ensemble object
!----------------------------------------------------------------------
subroutine ens_from_xz(ens,nam,geom,ne,ens_ga)

implicit none

! Passed variables
class(ens_type),intent(inout) :: ens                                       !< Ensemble
type(nam_type),intent(in) :: nam                                           !< Namelist
type(geom_type),intent(in) :: geom                                         !< Geometry
integer,intent(in) :: ne                                                   !< Ensemble size
real(kind_real),intent(in) :: ens_ga(geom%nmga,geom%nl0,nam%nv,nam%nts,ne) !< Ensemble on model grid, halo A

! Local variables
integer :: ie,its,iv,il0

! Allocation
call ens%alloc(nam,geom,ne,1)

! Copy
!$omp parallel do schedule(static) private(ie,its,iv,il0)
do ie=1,ens%ne
   do its=1,nam%nts
      do iv=1,nam%nv
         do il0=1,geom%nl0
            ens%fld(:,il0,iv,its,ie) = ens_ga(geom%c0a_to_ga,il0,iv,its,ie)
         end do
      end do
   end do
end do
!$omp end parallel do

end subroutine ens_from_xz

end module type_ens
