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
   integer :: nb                                   !< Number of blocks
   integer,allocatable :: nl0r(:)                  !< Effective number of levels
   integer,allocatable :: l0rl0b_to_l0(:,:,:)      !< Effective level to level
   integer,allocatable :: il0rz(:,:)               !< Effective zero separation level
   integer,allocatable :: nc3(:)                   !< Maximum class
   logical,allocatable :: auto_block(:)            !< Autocovariance block
   logical,allocatable :: diag_block(:)            !< HDIAG block
   logical,allocatable :: avg_block(:)             !< Averaging block
   logical,allocatable :: fit_block(:)             !< Fit block
   logical,allocatable :: B_block(:)               !< B-involved block
   logical,allocatable :: nicas_block(:)           !< NICAS block
   logical,allocatable :: cv_block(:)              !< Control variable block
   character(len=11),allocatable :: blockname(:)   !< Block name
   integer,allocatable :: b_to_v1(:)               !< Block to first variable
   integer,allocatable :: b_to_v2(:)               !< Block to second variable
   integer,allocatable :: b_to_ts1(:)              !< Block to first timeslot
   integer,allocatable :: b_to_ts2(:)              !< Block to second timeslot
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
allocate(bpar%l0rl0b_to_l0(nam%nl0r,geom%nl0,bpar%nb+1))
allocate(bpar%il0rz(geom%nl0,bpar%nb+1))
allocate(bpar%nl0r(bpar%nb+1))
allocate(bpar%nc3(bpar%nb+1))
allocate(bpar%auto_block(bpar%nb+1))
allocate(bpar%diag_block(bpar%nb+1))
allocate(bpar%avg_block(bpar%nb+1))
allocate(bpar%fit_block(bpar%nb+1))
allocate(bpar%B_block(bpar%nb+1))
allocate(bpar%nicas_block(bpar%nb+1))
allocate(bpar%cv_block(bpar%nb+1))
allocate(bpar%blockname(bpar%nb+1))
allocate(bpar%b_to_v1(bpar%nb))
allocate(bpar%b_to_v2(bpar%nb))
allocate(bpar%b_to_ts1(bpar%nb))
allocate(bpar%b_to_ts2(bpar%nb))

! Initialization
call msi(bpar%l0rl0b_to_l0)
call msi(bpar%il0rz)

if (nam%new_lct) then
   ! Individual blocks
   ib = 1
   do iv=1,nam%nv
      do its=1,nam%nts
         ! Classes and levels
         bpar%nl0r(ib) = nam%nl0r
         do jl0=1,geom%nl0
            il0off = jl0-(bpar%nl0r(ib)-1)/2-1
            if (il0off<1) il0off = 0
            if (il0off+nam%nl0r>geom%nl0) il0off = geom%nl0-nam%nl0r
            do il0r=1,nam%nl0r
               bpar%l0rl0b_to_l0(il0r,jl0,ib) = il0off+il0r
               if (bpar%l0rl0b_to_l0(il0r,jl0,ib)==jl0) bpar%il0rz(jl0,ib) = il0r
            end do
         end do
         bpar%nc3(ib) = nam%nc3
         bpar%nc3(ib) = nam%nc3
         bpar%auto_block(ib) = .true.
         bpar%diag_block(ib) = .true.
         bpar%avg_block(ib) = .false.
         bpar%fit_block(ib) = .false.
         bpar%B_block(ib) = .false.
         bpar%nicas_block(ib) = .false.
         bpar%cv_block(ib) = .false.

         ! Blocks information
         write(bpar%blockname(ib),'(i2.2,a,i2.2)') iv,'_',its
         bpar%b_to_v1(ib) = iv
         bpar%b_to_v2(ib) = iv
         bpar%b_to_ts1(ib) = its
         bpar%b_to_ts2(ib) = its

         ! Update block index
         ib = ib+1
      end do
   end do

   ! Common block
   ib = bpar%nb+1

   ! Classes and levels
   bpar%l0rl0b_to_l0(:,:,ib) = 0
   bpar%il0rz(:,ib) = 0
   bpar%nl0r(ib) = 0
   bpar%nc3(ib) = 0
   bpar%auto_block(ib) = .false.
   bpar%diag_block(ib) = .false.
   bpar%avg_block(ib) = .false.
   bpar%fit_block(ib) = .false.
   bpar%B_block(ib) = .false.
   bpar%nicas_block(ib) = .false.
   bpar%cv_block(ib) = .false.

   ! Blocks information
   bpar%blockname(ib) = 'common'
else
   ! Individual blocks
   ib = 1
   do iv=1,nam%nv
      do jv=1,nam%nv
         do its=1,nam%nts
            do jts=1,nam%nts
               ! Classes and levels
               if ((iv==jv).and.(its==jts)) then
                  bpar%nl0r(ib) = nam%nl0r
                  do jl0=1,geom%nl0
                     il0off = jl0-(bpar%nl0r(ib)-1)/2-1
                     if (il0off<1) il0off = 0
                     if (il0off+nam%nl0r>geom%nl0) il0off = geom%nl0-nam%nl0r
                     do il0r=1,nam%nl0r
                        bpar%l0rl0b_to_l0(il0r,jl0,ib) = il0off+il0r
                        if (bpar%l0rl0b_to_l0(il0r,jl0,ib)==jl0) bpar%il0rz(jl0,ib) = il0r
                     end do
                  end do
                  bpar%nc3(ib) = nam%nc3
               else
                  do jl0=1,geom%nl0
                     bpar%l0rl0b_to_l0(:,jl0,ib) = jl0
                  end do
                  bpar%il0rz(:,ib) = 1
                  bpar%nl0r(ib) = 1
                  bpar%nc3(ib) = 1
               end if

               ! Select blocks
               bpar%auto_block(ib) = (iv==jv).and.(its==jts)
               select case (nam%strategy)
               case ('common')
                  bpar%diag_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%avg_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%B_block(ib) = .false.
                  bpar%nicas_block(ib) = .false.
                  bpar%cv_block(ib) = .false.
               case ('specific_univariate')
                  bpar%diag_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%avg_block(ib) = .false.
                  bpar%B_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%nicas_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%cv_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
               case ('specific_multivariate')
                  bpar%diag_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%avg_block(ib) = .false.
                  bpar%B_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%nicas_block(ib) = (iv==jv).and.(its==1).and.(jts==1)
                  bpar%cv_block(ib) = .false.
               case ('common_weighted')
                  bpar%diag_block(ib) = .true.
                  bpar%avg_block(ib) = (iv==jv).and.(its==jts)
                  bpar%B_block(ib) = .true.
                  bpar%nicas_block(ib) = .false.
                  bpar%cv_block(ib) = (iv==jv).and.(its==jts)
               end select
               bpar%fit_block(ib) = bpar%diag_block(ib).and.(iv==jv).and.(its==jts).and.(trim(nam%fit_type)/='none')
               if (nam%local_diag) bpar%fit_block(ib) = bpar%fit_block(ib).and.bpar%nicas_block(ib)

               ! Blocks information
               write(bpar%blockname(ib),'(i2.2,a,i2.2,a,i2.2,a,i2.2)') iv,'_',jv,'_',its,'_',jts
               bpar%b_to_v1(ib) = iv
               bpar%b_to_v2(ib) = jv
               bpar%b_to_ts1(ib) = its
               bpar%b_to_ts2(ib) = jts

               ! Update block index
               ib = ib+1
            end do
         end do
      end do
   end do

   ! Common block
   ib = bpar%nb+1

   ! Classes and levels
   bpar%nl0r(ib) = nam%nl0r
   do jl0=1,geom%nl0
      il0off = jl0-(bpar%nl0r(ib)-1)/2-1
      if (il0off<1) il0off = 0
      if (il0off+nam%nl0r>geom%nl0) il0off = geom%nl0-nam%nl0r
      do il0r=1,nam%nl0r
         bpar%l0rl0b_to_l0(il0r,jl0,ib) = il0off+il0r
         if (bpar%l0rl0b_to_l0(il0r,jl0,ib)==jl0) bpar%il0rz(jl0,ib) = il0r
      end do
   end do
   bpar%nc3(ib) = nam%nc3

   ! Select blocks
   bpar%auto_block(ib) = .false.
   select case (nam%strategy)
   case ('common')
      bpar%diag_block(ib) = .true.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .true.
      bpar%nicas_block(ib) = .true.
      bpar%cv_block(ib) = .true.
   case ('specific_univariate')
      bpar%diag_block(ib) = .false.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .false.
      bpar%nicas_block(ib) = .false.
      bpar%cv_block(ib) = .false.
   case ('specific_multivariate')
      bpar%diag_block(ib) = .false.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .false.
      bpar%nicas_block(ib) = .false.
      bpar%cv_block(ib) = .true.
   case ('common_weighted')
      bpar%diag_block(ib) = .true.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .true.
      bpar%nicas_block(ib) = .true.
      bpar%cv_block(ib) = .false.
   end select
   bpar%fit_block(ib) = bpar%diag_block(ib).and.(trim(nam%fit_type)/='none')
   if (nam%local_diag) bpar%fit_block(ib) = bpar%fit_block(ib).and.bpar%nicas_block(ib)

   ! Blocks information
   bpar%blockname(ib) = 'common'
end if

end subroutine bpar_alloc

end module type_bpar
