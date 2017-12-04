!----------------------------------------------------------------------
! Module: type_bpar
!> Purpose: block parameters derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_bpar

use tools_kinds,only: kind_real
use tools_missing, only: msi,msr
use type_geom, only: geomtype
use type_nam, only: namtype

implicit none

type bpartype
   ! Block parameters
   integer :: nb                                  !< Number of blocks
   integer,allocatable :: nl0(:)                  !< Effective number of levels
   integer,allocatable :: il0rjl0ib_to_il0(:,:,:) !< Effective level to level
   integer,allocatable :: il0rz(:,:)              !< Effective zero separation level
   integer,allocatable :: icmax(:)                !< Maximum class
   logical,allocatable :: diag_block(:)           !< HDIAG block
   logical,allocatable :: avg_block(:)            !< Averaging block
   logical,allocatable :: fit_block(:)            !< Fit block
   logical,allocatable :: nicas_block(:)          !< NICAS block
   character(len=11),allocatable :: blockname(:)  !< Block name
   integer,allocatable :: ib_to_iv(:)             !< Block to first variable
   integer,allocatable :: ib_to_jv(:)             !< Block to second variable
   integer,allocatable :: ib_to_its(:)            !< Block to first timeslot
   integer,allocatable :: ib_to_jts(:)            !< Block to second timeslot
end type bpartype

private
public :: bpartype
public :: bpar_alloc

contains

!----------------------------------------------------------------------
! Subroutine: bpar_alloc
!> Purpose: allocate general parameters
!----------------------------------------------------------------------
subroutine bpar_alloc(nam,geom,bpar)

implicit none

! Passed variable
type(namtype),intent(in) :: nam      !< Namelist
type(geomtype),intent(in) :: geom    !< Geometry
type(bpartype),intent(inout) :: bpar !< Block parameters

! Local variables
integer :: ib,iv,jv,its,jts,il0r,jl0,il0off

! Number of blocks
if (nam%new_lct) then
   ! Special case for tensors
   bpar%nb = nam%nv*nam%nts
else
   ! General case
   bpar%nb = nam%nv**2*nam%nts**2
end if

! Allocation
allocate(bpar%il0rjl0ib_to_il0(nam%nl0r,geom%nl0,bpar%nb+1))
allocate(bpar%il0rz(geom%nl0,bpar%nb+1))
allocate(bpar%nl0(bpar%nb+1))
allocate(bpar%icmax(bpar%nb+1))
allocate(bpar%diag_block(bpar%nb+1))
allocate(bpar%avg_block(bpar%nb+1))
allocate(bpar%fit_block(bpar%nb+1))
allocate(bpar%nicas_block(bpar%nb+1))
allocate(bpar%blockname(bpar%nb+1))
allocate(bpar%ib_to_iv(bpar%nb))
allocate(bpar%ib_to_jv(bpar%nb))
allocate(bpar%ib_to_its(bpar%nb))
allocate(bpar%ib_to_jts(bpar%nb))

! Initialization
call msi(bpar%il0rjl0ib_to_il0)
call msi(bpar%il0rz)

if (nam%new_lct) then
   ! Individual blocks
   ib = 1
   do iv=1,nam%nv
      do its=1,nam%nts
         ! Classes and levels
         bpar%nl0(ib) = nam%nl0r
         do jl0=1,geom%nl0
            il0off = jl0-(bpar%nl0(ib)-1)/2-1
            if (il0off<1) il0off = 0
            if (il0off+nam%nl0r>geom%nl0) il0off = geom%nl0-nam%nl0r
            do il0r=1,nam%nl0r
               bpar%il0rjl0ib_to_il0(il0r,jl0,ib) = il0off+il0r
               if (bpar%il0rjl0ib_to_il0(il0r,jl0,ib)==jl0) bpar%il0rz(jl0,ib) = il0r
            end do
         end do
         bpar%icmax(ib) = nam%nc
         bpar%icmax(ib) = nam%nc
         bpar%diag_block(ib) = .true.
         bpar%avg_block(ib) = .false.
         bpar%nicas_block(ib) = .false.
         bpar%fit_block(ib) = .false.

         ! Blocks information
         write(bpar%blockname(ib),'(i2.2,a,i2.2)') iv,'_',its
         bpar%ib_to_iv(ib) = iv
         bpar%ib_to_jv(ib) = iv
         bpar%ib_to_its(ib) = its
         bpar%ib_to_jts(ib) = its

         ib = ib+1
      end do
   end do

   ! Common block

   ! Classes and levels
   bpar%il0rjl0ib_to_il0(:,:,bpar%nb+1) = 0
   bpar%il0rz(:,bpar%nb+1) = 0
   bpar%nl0(bpar%nb+1) = 0
   bpar%icmax(bpar%nb+1) = 0
   bpar%diag_block(bpar%nb+1) = .false.
   bpar%avg_block(bpar%nb+1) = .false.
   bpar%nicas_block(bpar%nb+1) = .false.

   ! Blocks information
   bpar%blockname(bpar%nb+1) = 'common'
   bpar%fit_block(bpar%nb+1) = .false.
else
   ! Individual blocks
   ib = 1
   do iv=1,nam%nv
      do jv=1,nam%nv
         do its=1,nam%nts
            do jts=1,nam%nts
               ! Classes and levels
               if ((iv==jv).and.(its==jts)) then
                  bpar%nl0(ib) = nam%nl0r
                  do jl0=1,geom%nl0
                     il0off = jl0-(bpar%nl0(ib)-1)/2-1
                     if (il0off<1) il0off = 0
                     if (il0off+nam%nl0r>geom%nl0) il0off = geom%nl0-nam%nl0r
                     do il0r=1,nam%nl0r
                        bpar%il0rjl0ib_to_il0(il0r,jl0,ib) = il0off+il0r
                        if (bpar%il0rjl0ib_to_il0(il0r,jl0,ib)==jl0) bpar%il0rz(jl0,ib) = il0r
                     end do
                  end do
                  bpar%icmax(ib) = nam%nc
               else
                  do jl0=1,geom%nl0
                     bpar%il0rjl0ib_to_il0(:,jl0,ib) = jl0
                  end do
                  bpar%il0rz(:,ib) = 1
                  bpar%nl0(ib) = 1
                  bpar%icmax(ib) = 1
               end if

               ! Select blocks
               select case (nam%strategy)
               case ('common')
                  bpar%diag_block(ib) = (iv==jv).and.(its==1)
                  bpar%avg_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%nicas_block(ib) = .false.
               case ('specific_univariate','specific_multivariate')
                  bpar%diag_block(ib) = (iv==jv).and.(its==1)
                  bpar%avg_block(ib) = .false.
                  bpar%nicas_block(ib) = (iv==jv).and.(its==1)
               case ('common_weighted')
                  bpar%diag_block(ib) = (its==1)
                  bpar%avg_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%nicas_block(ib) = .false.
               end select
               bpar%fit_block(ib) = bpar%diag_block(ib).and.(iv==jv).and.(its==jts).and.(trim(nam%fit_type)/='none')

               ! Blocks information
               write(bpar%blockname(ib),'(i2.2,a,i2.2,a,i2.2,a,i2.2)') iv,'_',jv,'_',its,'_',jts
               bpar%ib_to_iv(ib) = iv
               bpar%ib_to_jv(ib) = jv
               bpar%ib_to_its(ib) = its
               bpar%ib_to_jts(ib) = jts

               ! Update block index
               ib = ib+1
            end do
         end do
      end do
   end do

   ! Common block
   ib = bpar%nb+1

   ! Classes and levels
   bpar%nl0(ib) = nam%nl0r
   do jl0=1,geom%nl0
      il0off = jl0-(bpar%nl0(ib)-1)/2-1
      if (il0off<1) il0off = 0
      if (il0off+nam%nl0r>geom%nl0) il0off = geom%nl0-nam%nl0r
      do il0r=1,nam%nl0r
         bpar%il0rjl0ib_to_il0(il0r,jl0,ib) = il0off+il0r
         if (bpar%il0rjl0ib_to_il0(il0r,jl0,ib)==jl0) bpar%il0rz(jl0,ib) = il0r
      end do
   end do
   bpar%icmax(ib) = nam%nc

   ! Select blocks
   select case (nam%strategy)
   case ('common','common_weighted')
      bpar%diag_block(ib) = .true.
      bpar%avg_block(ib) = .false.
      bpar%nicas_block(ib) = .true.
   case ('specific_univariate','specific_multivariate')
      bpar%diag_block(ib) = .false.
      bpar%avg_block(ib) = .false.
      bpar%nicas_block(ib) = .false.
   end select

   ! Blocks information
   bpar%blockname(ib) = 'common'
   bpar%fit_block(ib) = bpar%diag_block(ib).and.(trim(nam%fit_type)/='none')
end if

end subroutine bpar_alloc

end module type_bpar
