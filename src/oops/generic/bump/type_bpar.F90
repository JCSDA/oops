!----------------------------------------------------------------------
! Module: type_bpar
! Purpose: block parameters derived type
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module type_bpar

use tools_kinds,only: kind_real
use type_geom, only: geom_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

type bpar_type
   ! Block parameters
   integer :: nb                                 ! Number of blocks
   integer :: nbe                                ! Extended number of blocks
   integer :: nl0rmax                            ! Maximum effective number of levels
   integer,allocatable :: nl0r(:)                ! Effective number of levels
   integer,allocatable :: l0rl0b_to_l0(:,:,:)    ! Effective level to level
   integer,allocatable :: il0rz(:,:)             ! Effective zero separation level
   integer,allocatable :: nc3(:)                 ! Maximum class
   logical,allocatable :: vbal_block(:,:)        ! Vertical balance block
   logical,allocatable :: diag_block(:)          ! HDIAG block
   logical,allocatable :: avg_block(:)           ! Averaging block
   logical,allocatable :: fit_block(:)           ! Fit block
   logical,allocatable :: B_block(:)             ! B-involved block
   logical,allocatable :: nicas_block(:)         ! NICAS block
   integer,allocatable :: cv_block(:)            ! Control variable block index
   character(len=11),allocatable :: blockname(:) ! Block name
   integer,allocatable :: b_to_v1(:)             ! Block to first variable
   integer,allocatable :: b_to_v2(:)             ! Block to second variable
   integer,allocatable :: b_to_ts1(:)            ! Block to first timeslot
   integer,allocatable :: b_to_ts2(:)            ! Block to second timeslot
contains
   procedure :: alloc => bpar_alloc
   procedure :: init => bpar_init
   procedure :: dealloc => bpar_dealloc
end type bpar_type

private
public :: bpar_type

contains

!----------------------------------------------------------------------
! Subroutine: bpar_alloc
! Purpose: allocation
!----------------------------------------------------------------------
subroutine bpar_alloc(bpar,nam,geom)

implicit none

! Passed variable
class(bpar_type),intent(inout) :: bpar ! Block parameters
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Number of blocks
bpar%nb = nam%nv**2*nam%nts**2
if (bpar%nb==1) then
    bpar%nbe = bpar%nb
else
   bpar%nbe = bpar%nb+1
end if

! Allocation
bpar%nl0rmax = min(nam%nl0r,geom%nl0)
allocate(bpar%l0rl0b_to_l0(bpar%nl0rmax,geom%nl0,bpar%nbe))
allocate(bpar%il0rz(geom%nl0,bpar%nbe))
allocate(bpar%nl0r(bpar%nbe))
allocate(bpar%nc3(bpar%nbe))
allocate(bpar%vbal_block(0:nam%nv,0:nam%nv))
allocate(bpar%diag_block(bpar%nbe))
allocate(bpar%avg_block(bpar%nbe))
allocate(bpar%fit_block(bpar%nbe))
allocate(bpar%B_block(bpar%nbe))
allocate(bpar%nicas_block(bpar%nbe))
allocate(bpar%cv_block(bpar%nbe))
allocate(bpar%blockname(bpar%nbe))
allocate(bpar%b_to_v1(bpar%nbe))
allocate(bpar%b_to_v2(bpar%nbe))
allocate(bpar%b_to_ts1(bpar%nbe))
allocate(bpar%b_to_ts2(bpar%nbe))

end subroutine bpar_alloc

!----------------------------------------------------------------------
! Subroutine: bpar_init
! Purpose: initialization
!----------------------------------------------------------------------
subroutine bpar_init(bpar,mpl,nam,geom)

implicit none

! Passed variable
class(bpar_type),intent(inout) :: bpar ! Block parameters
type(mpl_type),intent(inout) :: mpl    ! MPI data
type(nam_type),intent(in) :: nam       ! Namelist
type(geom_type),intent(in) :: geom     ! Geometry

! Local variables
integer :: ib,iv,jv,its,jts,il0,jl0r,jl0off

! Initialization
bpar%vbal_block = .false.

! Loop over variables and timeslots
ib = 1
do iv=1,nam%nv
   do jv=1,nam%nv
      do its=1,nam%nts
         do jts=1,nam%nts
            ! Classes and levels
            if ((trim(nam%strategy)=='diag_all').or.((iv==jv).and.(its==jts))) then
               bpar%nl0r(ib) = bpar%nl0rmax
               do il0=1,geom%nl0
                  jl0off = il0-(bpar%nl0r(ib)-1)/2-1
                  if (jl0off<1) jl0off = 0
                  if (jl0off+bpar%nl0rmax>geom%nl0) jl0off = geom%nl0-bpar%nl0rmax
                  do jl0r=1,bpar%nl0rmax
                     bpar%l0rl0b_to_l0(jl0r,il0,ib) = jl0off+jl0r
                     if (bpar%l0rl0b_to_l0(jl0r,il0,ib)==il0) bpar%il0rz(il0,ib) = jl0r
                  end do
               end do
               bpar%nc3(ib) = nam%nc3
            else
               do il0=1,geom%nl0
                  bpar%l0rl0b_to_l0(:,il0,ib) = il0
               end do
               bpar%il0rz(:,ib) = 1
               bpar%nl0r(ib) = 1
               bpar%nc3(ib) = 1
            end if

            ! Select blocks
            bpar%cv_block(ib) = mpl%msv%vali
            bpar%vbal_block(iv,jv) = (iv>1).and.(jv<iv).and.nam%vbal_block((iv-1)*(iv-2)/2+jv)
            select case (nam%strategy)
            case ('diag_all')
               bpar%diag_block(ib) = .true.
               bpar%avg_block(ib) = .true.
               bpar%B_block(ib) = .false.
               bpar%nicas_block(ib) = .false.
            case ('common')
               bpar%diag_block(ib) = (iv==jv).and.(its==jts)
               bpar%avg_block(ib) = (iv==jv).and.(its==jts)
               bpar%B_block(ib) = (ib==bpar%nbe)
               bpar%nicas_block(ib) = (ib==bpar%nbe)
               if (ib==bpar%nbe) bpar%cv_block(ib) = bpar%nbe
            case ('common_univariate')
               bpar%diag_block(ib) = (iv==jv).and.(its==jts)
               bpar%avg_block(ib) = (iv==jv).and.(its==jts)
               bpar%B_block(ib) = (ib==bpar%nbe)
               bpar%nicas_block(ib) = (ib==bpar%nbe)
               if ((iv==jv).and.(its==jts)) bpar%cv_block(ib) = bpar%nbe
            case ('common_weighted')
               bpar%diag_block(ib) = .true.
               bpar%avg_block(ib) = (iv==jv).and.(its==jts)
               bpar%B_block(ib) = .true.
               bpar%nicas_block(ib) = (bpar%nbe==bpar%nb)
               if ((iv==jv).and.(its==jts)) bpar%cv_block(ib) = bpar%nbe
            case ('specific_univariate')
               bpar%diag_block(ib) = (iv==jv).and.(its==jts)
               bpar%avg_block(ib) = .false.
               bpar%B_block(ib) = (iv==jv).and.(its==jts)
               bpar%nicas_block(ib) = (iv==jv).and.(its==jts)
               if ((iv==jv).and.(its==jts)) bpar%cv_block(ib) = ib
            case ('specific_multivariate')
               bpar%diag_block(ib) = (iv==jv).and.(its==jts)
               bpar%avg_block(ib) = (bpar%nbe==bpar%nb)
               bpar%B_block(ib) = (iv==jv).and.(its==jts)
               bpar%nicas_block(ib) = (iv==jv).and.(its==jts)
               if (ib==1) bpar%cv_block(ib) = 1
            end select
            bpar%fit_block(ib) = bpar%diag_block(ib).and.(iv==jv).and.(its==jts).and.(trim(nam%minim_algo)/='none')
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

if (bpar%nbe>bpar%nb) then
   ! Common block
   ib = bpar%nbe

   ! Classes and levels
   bpar%nl0r(ib) = bpar%nl0rmax
   do il0=1,geom%nl0
      jl0off = il0-(bpar%nl0r(ib)-1)/2-1
      if (jl0off<1) jl0off = 0
      if (jl0off+bpar%nl0rmax>geom%nl0) jl0off = geom%nl0-bpar%nl0rmax
      do jl0r=1,bpar%nl0rmax
         bpar%l0rl0b_to_l0(jl0r,il0,ib) = jl0off+jl0r
         if (bpar%l0rl0b_to_l0(jl0r,il0,ib)==il0) bpar%il0rz(il0,ib) = jl0r
      end do
   end do
   bpar%nc3(ib) = nam%nc3

   ! Select blocks
   bpar%cv_block(ib) = mpl%msv%vali
   select case (nam%strategy)
   case ('diag_all')
      bpar%diag_block(ib) = .true.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .false.
      bpar%nicas_block(ib) = .false.
   case ('common')
      bpar%diag_block(ib) = .true.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .true.
      bpar%nicas_block(ib) = .true.
      bpar%cv_block(ib) = ib
   case ('common_univariate')
      bpar%diag_block(ib) = .true.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .true.
      bpar%nicas_block(ib) = .true.
   case ('common_weighted')
      bpar%diag_block(ib) = .true.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .true.
      bpar%nicas_block(ib) = .true.
   case ('specific_univariate')
      bpar%diag_block(ib) = .false.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .false.
      bpar%nicas_block(ib) = .false.
   case ('specific_multivariate')
      bpar%diag_block(ib) = .false.
      bpar%avg_block(ib) = .false.
      bpar%B_block(ib) = .false.
      bpar%nicas_block(ib) = .false.
   end select
   bpar%fit_block(ib) = bpar%diag_block(ib).and.(trim(nam%minim_algo)/='none')
   if (nam%local_diag) bpar%fit_block(ib) = bpar%fit_block(ib).and.bpar%nicas_block(ib)

   ! Blocks information
   bpar%blockname(ib) = 'common'
   bpar%b_to_v1(ib) = 0
   bpar%b_to_v2(ib) = 0
   bpar%b_to_ts1(ib) = 0
   bpar%b_to_ts2(ib) = 0
end if

! Print summary
do ib=1,bpar%nbe
   iv = bpar%b_to_v1(ib)
   jv = bpar%b_to_v2(ib)
   if (bpar%vbal_block(iv,jv).or.bpar%diag_block(ib).or.bpar%avg_block(ib).or.bpar%fit_block(ib).or.bpar%B_block(ib) &
     & .or.bpar%nicas_block(ib).or.mpl%msv%isnoti(bpar%cv_block(ib))) then
      write(mpl%info,'(a7,a,a,a)') '','Block ',trim(bpar%blockname(ib)),':'
      call mpl%flush
      write(mpl%info,'(a10,a,i3)') '','Effective number of levels:    ',bpar%nl0r(ib)
      call mpl%flush
      write(mpl%info,'(a10,a,i3)') '','Maximum class:                 ',bpar%nc3(ib)
      call mpl%flush
      if ((iv>0).and.(jv>0)) then
         write(mpl%info,'(a10,a,l1)') '','Vertical balance block:          ',bpar%vbal_block(iv,jv)
         call mpl%flush
      end if
      write(mpl%info,'(a10,a,l1)') '','HDIAG block:                     ',bpar%diag_block(ib)
      call mpl%flush
      write(mpl%info,'(a10,a,l1)') '','Averaging block:                 ',bpar%avg_block(ib)
      call mpl%flush
      write(mpl%info,'(a10,a,l1)') '','Fit block:                       ',bpar%fit_block(ib)
      call mpl%flush
      write(mpl%info,'(a10,a,l1)') '','B-involved block:                ',bpar%B_block(ib)
      call mpl%flush
      write(mpl%info,'(a10,a,l1)') '','NICAS block:                     ',bpar%nicas_block(ib)
      call mpl%flush
      if (mpl%msv%isnoti(bpar%cv_block(ib))) then
         write(mpl%info,'(a10,a,i3)') '','Control variable block index:  ',bpar%cv_block(ib)
         call mpl%flush
      end if
   end if
end do

end subroutine bpar_init

!----------------------------------------------------------------------
! Subroutine: bpar_dealloc
! Purpose: release memory
!----------------------------------------------------------------------
subroutine bpar_dealloc(bpar)

implicit none

! Passed variable
class(bpar_type),intent(inout) :: bpar ! Block parameters

! Release memory
if (allocated(bpar%nl0r)) deallocate(bpar%nl0r)
if (allocated(bpar%l0rl0b_to_l0)) deallocate(bpar%l0rl0b_to_l0)
if (allocated(bpar%il0rz)) deallocate(bpar%il0rz)
if (allocated(bpar%nc3)) deallocate(bpar%nc3)
if (allocated(bpar%vbal_block)) deallocate(bpar%vbal_block)
if (allocated(bpar%diag_block)) deallocate(bpar%diag_block)
if (allocated(bpar%avg_block)) deallocate(bpar%avg_block)
if (allocated(bpar%fit_block)) deallocate(bpar%fit_block)
if (allocated(bpar%B_block)) deallocate(bpar%B_block)
if (allocated(bpar%nicas_block)) deallocate(bpar%nicas_block)
if (allocated(bpar%cv_block)) deallocate(bpar%cv_block)
if (allocated(bpar%blockname)) deallocate(bpar%blockname)
if (allocated(bpar%b_to_v1)) deallocate(bpar%b_to_v1)
if (allocated(bpar%b_to_v2)) deallocate(bpar%b_to_v2)
if (allocated(bpar%b_to_ts1)) deallocate(bpar%b_to_ts1)
if (allocated(bpar%b_to_ts2)) deallocate(bpar%b_to_ts2)

end subroutine bpar_dealloc

end module type_bpar
