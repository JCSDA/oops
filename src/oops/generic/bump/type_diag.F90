!----------------------------------------------------------------------
! Module: type_diag
!> Purpose: diagnostic derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_diag

use netcdf
use tools_const, only: reqkm,rad2deg,pi
use tools_fit, only: ver_smooth
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr,isnotmsi
use tools_nc, only: ncfloat
use type_avg, only: avg_type
use type_bpar, only: bpar_type
use type_diag_blk, only: diag_blk_type
use type_geom, only: geom_type
use type_hdata, only: hdata_type
use type_io, only: io_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type

implicit none

real(kind_real),parameter :: bound = 5.0_kind_real !< Restriction bound

! Diagnostic derived type
type diag_type
   character(len=1024) :: prefix               !< Prefix
   integer :: nc2a                             !< Number of local points
   type(diag_blk_type),allocatable :: blk(:,:) !< Diagnostic blocks
contains
   procedure :: alloc => diag_alloc
   procedure :: fit_filter => diag_fit_filter
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
! Subroutine: diag_fit_filter
!> Purpose: filter fit diagnostics
!----------------------------------------------------------------------
subroutine diag_fit_filter(diag,mpl,nam,geom,bpar,hdata)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(mpl_type),intent(in) :: mpl       !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(hdata_type),intent(in) :: hdata   !< HDIAG data

! Local variables
integer :: ib,il0,ic2a
real(kind_real) :: rh_c2a(hdata%nc2a,geom%nl0),rv_c2a(hdata%nc2a,geom%nl0)

do ib=1,bpar%nb+1
   if (bpar%fit_block(ib)) then
      ! Initialization
      call msr(rh_c2a)
      call msr(rv_c2a)

      do il0=1,geom%nl0
         do ic2a=1,hdata%nc2a
            ! Copy data
            rh_c2a(ic2a,il0) = diag%blk(ic2a,ib)%fit_rh(il0)
            rv_c2a(ic2a,il0) = diag%blk(ic2a,ib)%fit_rv(il0)

            ! Apply bounds
            if (isnotmsr(rh_c2a(ic2a,il0)).and.isnotmsr(diag%blk(0,ib)%fit_rh(il0))) then
               if ((rh_c2a(ic2a,il0)<diag%blk(0,ib)%fit_rh(il0)/bound) &
             & .or.(rh_c2a(ic2a,il0)>diag%blk(0,ib)%fit_rh(il0)*bound)) call msr(rh_c2a(ic2a,il0))
            end if
            if (isnotmsr(rv_c2a(ic2a,il0)).and.isnotmsr(diag%blk(0,ib)%fit_rv(il0))) then
               if ((rv_c2a(ic2a,il0)<diag%blk(0,ib)%fit_rv(il0)/bound) &
             & .or.(rv_c2a(ic2a,il0)>diag%blk(0,ib)%fit_rv(il0)*bound)) call msr(rv_c2a(ic2a,il0))
            end if
         end do

         ! Median filter to remove extreme values
         call hdata%diag_filter(mpl,nam,geom,il0,'median',nam%diag_rhflt,rh_c2a(:,il0))
         call hdata%diag_filter(mpl,nam,geom,il0,'median',nam%diag_rhflt,rv_c2a(:,il0))

         ! Average filter to smooth support radii
         call hdata%diag_filter(mpl,nam,geom,il0,'average',nam%diag_rhflt,rh_c2a(:,il0))
         call hdata%diag_filter(mpl,nam,geom,il0,'average',nam%diag_rhflt,rv_c2a(:,il0))

         ! Fill missing values
         call hdata%diag_fill(mpl,geom,il0,rh_c2a(:,il0))
         call hdata%diag_fill(mpl,geom,il0,rv_c2a(:,il0))

         ! Copy data
         do ic2a=1,hdata%nc2a
            diag%blk(ic2a,ib)%fit_rh(il0) = rh_c2a(ic2a,il0)
            diag%blk(ic2a,ib)%fit_rv(il0) = rv_c2a(ic2a,il0)
         end do
      end do

      ! Smooth vertically
      do ic2a=0,diag%nc2a
         call ver_smooth(mpl,geom%nl0,geom%vunitavg,nam%rvflt,diag%blk(ic2a,ib)%fit_rh)
         call ver_smooth(mpl,geom%nl0,geom%vunitavg,nam%rvflt,diag%blk(ic2a,ib)%fit_rv)
      end do
   end if
end do

end subroutine diag_fit_filter

!----------------------------------------------------------------------
! Subroutine: diag_write
!> Purpose: write all diagnostics
!----------------------------------------------------------------------
subroutine diag_write(diag,mpl,nam,geom,bpar,io,hdata)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(mpl_type),intent(inout) :: mpl    !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(io_type),intent(in) :: io         !< I/O
type(hdata_type),intent(in) :: hdata   !< HDIAG data

! Local variables
integer :: ib,i,ic2,il0,il0i,iproc,ic2a,ildw
real(kind_real) :: fld_c2a(hdata%nc2a,geom%nl0),fld_c2b(hdata%nc2b,geom%nl0),fld_c0a(geom%nc0a,geom%nl0)
character(len=7) :: lonchar,latchar
character(len=1024) :: filename

if (mpl%main) then
   filename = trim(nam%prefix)//'_diag.nc'
   do ib=1,bpar%nb+1
      if (bpar%diag_block(ib)) then
        call diag%blk(0,ib)%write(mpl,nam,geom,bpar,filename)
      end if
   end do
end if

if (nam%local_diag) then
   do ib=1,bpar%nb+1
      if (bpar%fit_block(ib)) then
         filename = trim(nam%prefix)//'_local_diag_'//trim(diag%prefix)
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
            call hdata%com_AB%ext(mpl,geom%nl0,fld_c2a,fld_c2b)
            do il0=1,geom%nl0
               il0i = min(il0,geom%nl0i)
               call hdata%h(il0i)%apply(mpl,fld_c2b(:,il0),fld_c0a(:,il0))
            end do

            ! Write fields
            if (i==1) then
               call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rh',fld_c0a)
            elseif (i==2) then
               call io%fld_write(mpl,nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rv',fld_c0a)
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
         write(lonchar,'(f7.2)') nam%lon_ldwv(ildw)*rad2deg
         write(latchar,'(f7.2)') nam%lat_ldwv(ildw)*rad2deg
         filename = trim(nam%prefix)//'_diag_'//trim(adjustl(lonchar))//'-'//trim(adjustl(latchar))//'.nc'

         ! Find diagnostic point task
         ic2a = hdata%c2_to_c2a(ic2)
         do ib=1,bpar%nb+1
            if (bpar%diag_block(ib)) call diag%blk(ic2a,ib)%write(mpl,nam,geom,bpar,filename)
         end do
      end if
   else
      call mpl%warning('missing local profile')
   end if
end do

end subroutine diag_write

!----------------------------------------------------------------------
! Subroutine: diag_covariance
!> Purpose: compute covariance
!----------------------------------------------------------------------
subroutine diag_covariance(diag,mpl,nam,geom,bpar,io,hdata,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(mpl_type),intent(inout) :: mpl    !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(io_type),intent(in) :: io         !< I/O
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(avg_type),intent(in) :: avg       !< Averaged statistics
character(len=*),intent(in) :: prefix  !< Diagnostic prefix

! Local variables
integer :: ib,ic2a,il0

! Allocation
call diag%alloc(nam,geom,bpar,hdata,prefix)

do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)') '','Block ',trim(bpar%blockname(ib))
      call flush(mpl%unit)

      do ic2a=0,diag%nc2a
         ! Copy
         diag%blk(ic2a,ib)%raw = avg%blk(ic2a,ib)%m11
      end do

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw(1,bpar%il0rz(il0,ib),il0))) then
            write(mpl%unit,'(a13,a,i3,a,a,e9.2,a)') '','Level: ',nam%levs(il0),' ~> cov. at class zero: ',trim(mpl%peach), &
          & diag%blk(0,ib)%raw(1,bpar%il0rz(il0,ib),il0),trim(mpl%black)
            call flush(mpl%unit)
         end if
      end do
   end if
end do

! Write
call diag%write(mpl,nam,geom,bpar,io,hdata)

end subroutine diag_covariance

!----------------------------------------------------------------------
! Subroutine: diag_correlation
!> Purpose: compute correlation
!----------------------------------------------------------------------
subroutine diag_correlation(diag,mpl,nam,geom,bpar,io,hdata,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(mpl_type),intent(inout) :: mpl    !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(io_type),intent(in) :: io         !< I/O
type(hdata_type),intent(in) :: hdata   !< HDIAG data
type(avg_type),intent(in) :: avg       !< Averaged statistics
character(len=*),intent(in) :: prefix  !< Diagnostic prefix

! Local variables
integer :: ib,ic2a,progint,il0
logical,allocatable :: done(:)
type(diag_type) :: ndiag

! Allocation
call diag%alloc(nam,geom,bpar,hdata,prefix)
call ndiag%alloc(nam,geom,bpar,hdata,'n'//trim(prefix))
allocate(done(0:diag%nc2a))

do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      write(mpl%unit,'(a10,a,a,a)',advance='no') '','Block ',trim(bpar%blockname(ib)),':'
      call flush(mpl%unit)

      ! Initialization
      call mpl%prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Copy
         diag%blk(ic2a,ib)%raw = avg%blk(ic2a,ib)%cor

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(mpl,nam,geom,bpar,hdata)

         ! Update
         done(ic2a) = .true.
         call mpl%prog_print(progint,done)
      end do
      ndiag%blk(0,ib)%raw = avg%blk(0,ib)%nc1a_cor
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(avg%blk(0,ib)%cor(1,bpar%il0rz(il0,ib),il0))) then
            write(mpl%unit,'(a13,a,i3,a4,a20,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> ','cor. at class zero: ', &
          & trim(mpl%peach),avg%blk(0,ib)%cor(1,bpar%il0rz(il0,ib),il0),trim(mpl%black)
            call flush(mpl%unit)
         end if
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%unit,'(a47,a,f10.2,a,f10.2,a)') 'cor. support radii: ',trim(mpl%aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm, &
             & trim(mpl%black)//' km  / '//trim(mpl%aqua),diag%blk(0,ib)%fit_rv(il0),trim(mpl%black)//' '//trim(mpl%vunitchar)
               call flush(mpl%unit)
            end if
         end if
      end do
   end if
end do

! Filtering
if (nam%local_diag) call diag%fit_filter(mpl,nam,geom,bpar,hdata)

! Write
call diag%write(mpl,nam,geom,bpar,io,hdata)
call ndiag%write(mpl,nam,geom,bpar,io,hdata)

end subroutine diag_correlation

!----------------------------------------------------------------------
! Subroutine: diag_localization
!> Purpose: compute diagnostic _localization
!----------------------------------------------------------------------
subroutine diag_localization(diag,mpl,nam,geom,bpar,io,hdata,avg,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic
type(mpl_type),intent(inout) :: mpl    !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(io_type),intent(in) :: io         !< I/O
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
      call flush(mpl%unit)

      ! Initialization
      call mpl%prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Compute localization
         call diag%blk(ic2a,ib)%localization(geom,bpar,avg%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(mpl,nam,geom,bpar,hdata)

         ! Update
         done(ic2a) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) then
            write(mpl%unit,'(a13,a,i3,a4,a20,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> ','loc. at class zero: ', &
          & trim(mpl%peach),diag%blk(0,ib)%raw_coef_ens(il0),trim(mpl%black)
            call flush(mpl%unit)
         end if
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%unit,'(a47,a,f10.2,a,f10.2,a)') 'loc. support radii: ',trim(mpl%aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm, &
             & trim(mpl%black)//' km  / '//trim(mpl%aqua),diag%blk(0,ib)%fit_rv(il0),trim(mpl%black)//' '//trim(mpl%vunitchar)
               call flush(mpl%unit)
            end if
         end if
      end do
   end if
end do

! Filtering
if (nam%local_diag) call diag%fit_filter(mpl,nam,geom,bpar,hdata)

! Write
call diag%write(mpl,nam,geom,bpar,io,hdata)

end subroutine diag_localization

!----------------------------------------------------------------------
! Subroutine: diag_hybridization
!> Purpose: compute diagnostic hybridization
!----------------------------------------------------------------------
subroutine diag_hybridization(diag,mpl,nam,geom,bpar,io,hdata,avg,avg_sta,prefix)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag !< Diagnostic (localization)
type(mpl_type),intent(inout) :: mpl    !< MPI data
type(nam_type),intent(in) :: nam       !< Namelist
type(geom_type),intent(in) :: geom     !< Geometry
type(bpar_type),intent(in) :: bpar     !< Block parameters
type(io_type),intent(in) :: io         !< I/O
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
      call flush(mpl%unit)

      ! Initialization
      call mpl%prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Compute hybridization
         call diag%blk(ic2a,ib)%hybridization(geom,bpar,avg%blk(ic2a,ib),avg_sta%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) call diag%blk(ic2a,ib)%fitting(mpl,nam,geom,bpar,hdata)

         ! Update
         done(ic2a) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) then
            write(mpl%unit,'(a13,a,i3,a4,a21,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> ','loc. at class zero: ', &
          & trim(mpl%peach),diag%blk(0,ib)%raw_coef_ens(il0),trim(mpl%black)
            call flush(mpl%unit)
         end if
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%unit,'(a48,a,f10.2,a,f10.2,a)') 'loc. support radii: ',trim(mpl%aqua),diag%blk(0,ib)%fit_rh(il0)*reqkm, &
             & trim(mpl%black)//' km  / '//trim(mpl%aqua),diag%blk(0,ib)%fit_rv(il0),trim(mpl%black)//' '//trim(mpl%vunitchar)
               call flush(mpl%unit)
            end if
         end if
      end do
      write(mpl%unit,'(a13,a,f10.2,a)') '','Raw static coeff.: ',trim(mpl%purple),diag%blk(0,ib)%raw_coef_sta,trim(mpl%black)
      call flush(mpl%unit)
   end if
end do

! Filtering
if (nam%local_diag) call diag%fit_filter(mpl,nam,geom,bpar,hdata)

! Write
call diag%write(mpl,nam,geom,bpar,io,hdata)

end subroutine diag_hybridization

!----------------------------------------------------------------------
! Subroutine: diag_dualens
!> Purpose: compute diagnostic dualens
!----------------------------------------------------------------------
subroutine diag_dualens(diag,mpl,nam,geom,bpar,io,hdata,avg,avg_lr,diag_lr,prefix,prefix_lr)

implicit none

! Passed variables
class(diag_type),intent(inout) :: diag   !< Diagnostic (localization)
type(mpl_type),intent(inout) :: mpl    !< MPI data
type(nam_type),intent(in) :: nam         !< Namelist
type(geom_type),intent(in) :: geom       !< Geometry
type(bpar_type),intent(in) :: bpar       !< Block parameters
type(io_type),intent(in) :: io           !< I/O
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
      call flush(mpl%unit)

      ! Initialization
      call mpl%prog_init(progint,done)

      do ic2a=0,diag%nc2a
         ! Compute dualens
         call diag%blk(ic2a,ib)%dualens(geom,bpar,avg%blk(ic2a,ib),avg_lr%blk(ic2a,ib),diag_lr%blk(ic2a,ib))

         ! Normalization
         call diag%blk(ic2a,ib)%normalization(geom,bpar)
         call diag_lr%blk(ic2a,ib)%normalization(geom,bpar)

         ! Fitting
         if (bpar%fit_block(ib)) then
            call diag%blk(ic2a,ib)%fitting(mpl,nam,geom,bpar,hdata)
            call diag_lr%blk(ic2a,ib)%fitting(mpl,nam,geom,bpar,hdata)
         end if

         ! Update
         done(ic2a) = .true.
         call mpl%prog_print(progint,done)
      end do
      write(mpl%unit,'(a)') '100%'
      call flush(mpl%unit)

      ! Print results
      do il0=1,geom%nl0
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) then
            write(mpl%unit,'(a10,a,i3,a4,a21,a,f10.2,a)') '','Level: ',nam%levs(il0),' ~> ','loc. at class zero (HR): ', &
          & trim(mpl%peach),diag%blk(0,ib)%raw_coef_ens(il0),trim(mpl%black)
            call flush(mpl%unit)
         end if
         if (isnotmsr(diag%blk(0,ib)%raw_coef_ens(il0))) then
            write(mpl%unit,'(a45,a,f10.2,a)') 'loc. at class zero (LR): ',trim(mpl%peach),diag_lr%blk(0,ib)%raw_coef_ens(il0), &
          & trim(mpl%black)
            call flush(mpl%unit)
         end if
         if (bpar%fit_block(ib)) then
            if (isnotmsr(diag%blk(0,ib)%fit_rh(il0))) then
               write(mpl%unit,'(a45,a,f10.2,a,f10.2,a)') 'loc. support radii (HR): ',trim(mpl%aqua), &
             & diag%blk(0,ib)%fit_rh(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),diag_lr%blk(0,ib)%fit_rv(il0), &
             & trim(mpl%black)//' '//trim(mpl%vunitchar)
               call flush(mpl%unit)
            end if
            if (isnotmsr(diag_lr%blk(0,ib)%fit_rh(il0))) then
               write(mpl%unit,'(a45,a,f10.2,a,f10.2,a)') 'loc. support radii (LR): ',trim(mpl%aqua), &
             & diag_lr%blk(0,ib)%fit_rh(il0)*reqkm,trim(mpl%black)//' km  / '//trim(mpl%aqua),diag_lr%blk(0,ib)%fit_rv(il0), &
             & trim(mpl%black)//' '//trim(mpl%vunitchar)
               call flush(mpl%unit)
            end if
         end if
      end do
   end if
end do

! Filtering
if (nam%local_diag) then
   call diag%fit_filter(mpl,nam,geom,bpar,hdata)
   call diag_lr%fit_filter(mpl,nam,geom,bpar,hdata)
end if

! Write
call diag%write(mpl,nam,geom,bpar,io,hdata)
call diag_lr%write(mpl,nam,geom,bpar,io,hdata)

end subroutine diag_dualens

end module type_diag
