!----------------------------------------------------------------------
! Module: type_diag
!> Purpose: diagnostic derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_diag

use model_interface, only: model_write
use netcdf
use tools_const, only: reqkm
use tools_display, only: vunitchar,msgerror,msgwarning,prog_init,prog_print,aqua,peach,purple,black
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isnotmsr,isallnotmsr,isnotmsi
use tools_nc, only: ncerr,ncfloat
use type_avg, only: avg_type
use type_bpar, only: bpar_type
use type_diag_blk, only: diag_blk_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_mpl, only: mpl
use type_nam, only: nam_type

implicit none

! Diagnostic derived type
type diag_type
   character(len=1024) :: prefix               !< Prefix
   integer :: nc2a                             !< Number of local points
   type(diag_blk_type),allocatable :: blk(:,:) !< Diagnostic blocks
contains
   procedure :: alloc => diag_alloc
   procedure :: write => diag_write
   procedure :: covariance => diag_covariance
   procedure :: correlation => diag_correlation
   procedure :: localization => diag_localization
   procedure :: hybridization => diag_hybridization
   procedure :: dualens => diag_dualens
end type diag_type

private
public :: diag_type

contains

!----------------------------------------------------------------------
! Subroutine: diag_alloc
!> Purpose: allocation
!----------------------------------------------------------------------
subroutine diag_alloc(diag,nam,geom,bpar,hdata,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
character(len=*),intent(in) :: prefix  !< Block prefix

! Local variables
integer :: ib,ic2a

! Number of local points
if (nam%local_diag) then
   diag%nc2a = hdata%nc2a
else
   diag%nc2a = 0
end if

! Prefix
diag%prefix = trim(prefix)

! Allocation
allocate(diag%blk(0:diag%nc2a,bpar%nb+1))
do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      do ic2a=0,diag%nc2a
         call diag%blk(ic2a,ib)%alloc(nam,geom,bpar,ic2a,ib,prefix)
      end do
   end if
end do

end subroutine diag_alloc

!----------------------------------------------------------------------
! Subroutine: diag_write
!> Purpose: write all diagnostics
!----------------------------------------------------------------------
subroutine diag_write(diag,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data

! Local variables
integer :: ib,i,ic2,il0,il0i,iproc,ic2a,ildw
real(kind_real) :: fld_c2a(hdata%nc2a,geom%nl0),fld_c2b(hdata%nc2b,geom%nl0),fld_c0a(geom%nc0a,geom%nl0)
character(len=7) :: lonchar,latchar
character(len=1024) :: filename,filename_gridded

if (mpl%main) then
   filename = trim(nam%prefix)//'_diag.nc'
   do ib=1,bpar%nb+1
      if (bpar%diag_block(ib)) then
        call diag%blk(0,ib)%write(nam,geom,bpar,filename)
      end if
   end do
end if

if (nam%local_diag) then
   do ib=1,bpar%nb+1
      if (bpar%fit_block(ib)) then
         filename = trim(nam%prefix)//'_local_diag_'//trim(diag%prefix)//'.nc'
         filename_gridded = trim(nam%prefix)//'_local_diag_'//trim(diag%prefix)//'_gridded.nc'
         do i=1,2
            ! Copy data
            do ic2a=1,hdata%nc2a
               if (i==1) then
                  fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rh*reqkm
               elseif (i==2) then
                  fld_c2a(ic2a,:) = diag%blk(ic2a,ib)%fit_rv
               end if
            end do

            ! Interpolation
            call hdata%com_AB%ext(geom%nl0,fld_c2a,fld_c2b)
            do il0=1,geom%nl0
               il0i = min(il0,geom%nl0i)
               call hdata%h(il0i)%apply(fld_c2b(:,il0),fld_c0a(:,il0))
            end do

            ! Write fields
            if (i==1) then
               call geom%fld_write(nam,filename,trim(bpar%blockname(ib))//'_fit_rh',fld_c0a)
               call model_write(nam,geom,filename_gridded,trim(bpar%blockname(ib))//'_fit_rh',fld_c0a)
            elseif (i==2) then
               call geom%fld_write(nam,filename,trim(bpar%blockname(ib))//'_fit_rv',fld_c0a)
               call model_write(nam,geom,filename_gridded,trim(bpar%blockname(ib))//'_fit_rv',fld_c0a)
            end if
         end do
      end if
   end do
end if

do ildw=1,nam%nldwv
   if (isnotmsi(hdata%nn_ldwv_index(ildw))) then
      ic2 = hdata%nn_ldwv_index(ildw)
      iproc = hdata%c2_to_proc(ic2)
      if (mpl%myproc==iproc) then
         ! Build file name
         write(lonchar,'(f7.2)') nam%lon_ldwv(ildw)
         write(latchar,'(f7.2)') nam%lat_ldwv(ildw)
         filename = trim(nam%prefix)//'_diag_'//trim(adjustl(lonchar))//'-'//trim(adjustl(latchar))//'.nc'

         ! Find diagnostic point task
         ic2a = hdata%c2_to_c2a(ic2)
         do ib=1,bpar%nb+1
            if (bpar%diag_block(ib)) call diag%blk(ic2a,ib)%write(nam,geom,bpar,filename)
         end do
      end if
   else
      call msgwarning('missing local profile')
   end if
end do

end subroutine diag_write

!----------------------------------------------------------------------
! Subroutine: diag_covariance
!> Purpose: compute covariance
!----------------------------------------------------------------------
subroutine diag_covariance(diag,nam,geom,bpar,hdata,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(avg_type),intent(in) :: avg       !< Averaged statistics
character(len=*),intent(in) :: prefix  !< Diagnostic prefix

! Local variables
integer :: ib,ic2a,il0
logical,allocatable :: done(:)

! Allocation
call diag%alloc(nam,geom,bpar,hdata,prefix)
allocate(done(0:diag%nc2a))

do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib))

      do ic2a=0,diag%nc2a
         ! Copy
         diag%blk(ic2a,ib)%raw = avg%blk(ic2a,ib)%m11
      end do

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw(1,bpar%il0rz(il0,ib),il0))) write(mpl%unit,'(a13,a,i3,a,a,e9.2,a)') '','Level: ', &
       & nam%levs(il0),' ~> cov. at class zero: ',trim(peach),diag%blk(0,ib)%raw(1,bpar%il0rz(il0,ib),il0),trim(black)
      end do
   end if
end do

! Write
call diag%write(nam,geom,bpar,hdata)

end subroutine diag_covariance

!----------------------------------------------------------------------
! Subroutine: diag_correlation
!> Purpose: compute correlation
!----------------------------------------------------------------------
subroutine diag_correlation(diag,nam,geom,bpar,hdata,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(avg_type),intent(in) :: avg       !< Averaged statistics
character(len=*),intent(in) :: prefix  !< Diagnostic prefix

! Local variables
integer :: ib,ic2a,progint,il0
logical,allocatable :: done(:)

! Allocation
call diag%alloc(nam,geom,bpar,hdata,prefix)
allocate(done(0:diag%nc2a))

do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'

      ! Initialization
      call prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Copy
         diag%blk(ic2a,ib)%raw = avg%blk(ic2a,ib)%cor

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(nam,geom,bpar)

         ! Update
         done(ic2a) = .true.
         call prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'


      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(avg%blk(0,ib)%cor(1,bpar%il0rz(il0,ib),il0))) write(mpl%unit,'(a13,a,i3,a4,a20,a,f10.2,a)') '', &
       & 'Level: ',nam%levs(il0),' ~> ','cor. at class zero: ',trim(peach),avg%blk(0,ib)%cor(1,bpar%il0rz(il0,ib),il0),trim(black)
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) write(mpl%unit,'(a47,a,f10.2,a,f10.2,a)') 'cor. support radii: ', &
          & trim(aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(black)//' km  / '//trim(aqua),diag%blk(0,ib)%fit_rv(il0), &
          & trim(black)//' '//trim(vunitchar)
         end if
      end do
   end if
end do

! Write
call diag%write(nam,geom,bpar,hdata)

end subroutine diag_correlation

!----------------------------------------------------------------------
! Subroutine: diag_localization
!> Purpose: compute diagnostic _localization
!----------------------------------------------------------------------
subroutine diag_localization(diag,nam,geom,bpar,hdata,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(avg_type),intent(in) :: avg       !< Averaged statistics
character(len=*),intent(in) :: prefix  !< Block prefix

! Local variables
integer :: ib,ic2a,progint,il0
logical,allocatable :: done(:)

! Allocation
call diag%alloc(nam,geom,bpar,hdata,prefix)
allocate(done(0:diag%nc2a))

do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'

      ! Initialization
      call prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Compute localization
         call diag%blk(ic2a,ib)%localization(geom,bpar,avg%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(nam,geom,bpar)

         ! Update
         done(ic2a) = .true.
         call prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) write(mpl%unit,'(a13,a,i3,a4,a20,a,f10.2,a)') '','Level: ', &
       & nam%levs(il0),' ~> ','loc. at class zero: ',trim(peach),diag%blk(0,ib)%raw_coef_ens(il0),trim(black)
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) write(mpl%unit,'(a47,a,f10.2,a,f10.2,a)') 'loc. support radii: ', &
          & trim(aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(black)//' km  / '//trim(aqua),diag%blk(0,ib)%fit_rv(il0), &
          & trim(black)//' '//trim(vunitchar)
         end if
      end do
   end if
end do

! Write
call diag%write(nam,geom,bpar,hdata)

end subroutine diag_localization

!----------------------------------------------------------------------
! Subroutine: diag_hybridization
!> Purpose: compute diagnostic _hybridization
!----------------------------------------------------------------------
subroutine diag_hybridization(diag,nam,geom,bpar,hdata,avg,avg_sta,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic (localization)
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(avg_type),intent(in) :: avg       !< Averaged statistics
type(avg_type),intent(in) :: avg_sta   !< Static averaged statistics
character(len=*),intent(in) :: prefix  !< Diagnostic prefix

! Local variables
integer :: ib,ic2a,progint,il0
logical,allocatable :: done(:)

! Allocation
call diag%alloc(nam,geom,bpar,hdata,prefix)
allocate(done(0:diag%nc2a))

do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'

      ! Initialization
      call prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Compute hybridization
         call diag%blk(ic2a,ib)%hybridization(geom,bpar,avg%blk(ic2a,ib),avg_sta%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(nam,geom,bpar)

         ! Update
         done(ic2a) = .true.
         call prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) write(mpl%unit,'(a13,a,i3,a4,a21,a,f10.2,a)') '','Level: ', &
       & nam%levs(il0),' ~> ','loc. at class zero: ',trim(peach),diag%blk(0,ib)%raw_coef_ens(il0),trim(black)
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) write(mpl%unit,'(a48,a,f10.2,a,f10.2,a)') 'loc. support radii: ', &
          & trim(aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(black)//' km  / '//trim(aqua),diag%blk(0,ib)%fit_rv(il0), &
          & trim(black)//' '//trim(vunitchar)
         end if
      end do
      write(mpl%unit,'(a13,a,f10.2,a)') '','Raw static coeff.: ',trim(purple),diag%blk(0,ib)%raw_coef_sta,trim(black)
   end if
end do

! Write
call diag%write(nam,geom,bpar,hdata)

end subroutine diag_hybridization

!----------------------------------------------------------------------
! Subroutine: diag_dualens
!> Purpose: compute diagnostic _dualens
!----------------------------------------------------------------------
subroutine diag_dualens(diag,nam,geom,bpar,hdata,avg,avg_lr,diag_lr,prefix,prefix_lr)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag   !< Diagnostic (localization)
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
type(bpar_type),intent(in) :: bpar       !< Block parameters
type(hdata_type),intent(in) :: hdata     !< HDIAG data
type(avg_type),intent(in) :: avg         !< Averaged statistics
type(avg_type),intent(in) :: avg_lr      !< LR averaged statistics
type(diag_type),intent(inout) :: diag_lr !< Diagnostic (LR localization)
character(len=*),intent(in) :: prefix    !< Diagnostic prefix
character(len=*),intent(in) :: prefix_lr !< LR diagnostic prefix

! Local variables
integer :: ib,ic2a,progint,il0
logical,allocatable :: done(:)

! Allocation
call diag%alloc(nam,geom,bpar,hdata,prefix)
call diag_lr%alloc(nam,geom,bpar,hdata,prefix_lr)
allocate(done(0:diag%nc2a))

do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'

      ! Initialization
      call prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Compute dualens
         call diag%blk(ic2a,ib)%dualens(geom,bpar,avg%blk(ic2a,ib),avg_lr%blk(ic2a,ib),diag_lr%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)
         call diag_lr%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) then
            call diag%blk(ic2a,ib)%fitting(nam,geom,bpar)
            call diag_lr%blk(ic2a,ib)%fitting(nam,geom,bpar)
         end if

         ! Update
         done(ic2a) = .true.
         call prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) write(mpl%unit,'(a10,a,i3,a4,a21,a,f10.2,a)') '','Level: ', &
       & nam%levs(il0),' ~> ','loc. at class zero (HR): ',trim(peach),diag%blk(0,ib)%raw_coef_ens(il0),trim(black)
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) write(mpl%unit,'(a45,a,f10.2,a)') 'loc. at class zero (LR): ', &
       & trim(peach),diag_lr%blk(0,ib)%raw_coef_ens(il0),trim(black)
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) write(mpl%unit,'(a45,a,f10.2,a,f10.2,a)') 'loc. support radii (HR): ', &
          & trim(aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(black)//' km  / '//trim(aqua),diag_lr%blk(0,ib)%fit_rv(il0), &
          & trim(black)//' '//trim(vunitchar)
            if (isnotmsr(diag_lr%blk(0,ib)%fit_rh(il0))) write(mpl%unit,'(a45,a,f10.2,a,f10.2,a)') 'loc. support radii (LR): ', &
          & trim(aqua),diag_lr%blk(0,ib)%fit_rh(il0)*reqkm,trim(black)//' km  / '//trim(aqua),diag_lr%blk(0,ib)%fit_rv(il0), &
          & trim(black)//' '//trim(vunitchar)
         end if
      end do
   end if
end do

! Write
call diag%write(nam,geom,bpar,hdata)

end subroutine diag_dualens

end module type_diag
