!----------------------------------------------------------------------
! Module: module_diag_tools.f90
!> Purpose: diagnostics tools routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_diag_tools

use netcdf
use omp_lib
use tools_const, only: req,sphere_dist,rad2deg,gc99,median
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isnotmsi,isnotmsr,isanynotmsr,isallnotmsr
use tools_nc, only: ncfloat,ncerr
use type_hdata, only: hdatatype
use type_linop, only: apply_linop
use type_mpl, only: mpl
implicit none

interface diag_filter
  module procedure diag_filter_2d
  module procedure diag_filter_3d
end interface
interface diag_interpolation
  module procedure diag_interpolation_2d
  module procedure diag_interpolation_3d
end interface
interface diag_write
  module procedure diag_write_2d
  module procedure diag_write_3d
end interface

private
public :: diag_filter,diag_interpolation,diag_write

contains

!----------------------------------------------------------------------
! Subroutine: diag_filter_2d
!> Purpose: filter diagnostics, on a single level
!----------------------------------------------------------------------
subroutine diag_filter_2d(hdata,il0,filter_type,r,diag)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata              !< HDIAG data
integer,intent(in) :: il0                        !< Level
character(len=*),intent(in) :: filter_type       !< Filter type
real(kind_real),intent(in) :: r                  !< Filter support radius
real(kind_real),intent(inout) :: diag(hdata%nc2) !< Filtered diagnostics

! Local variables
integer :: ic2,jc2,nc2eff
real(kind_real) :: diag_tmp(hdata%nc2),distnorm,norm,wgt
real(kind_real),allocatable :: list(:),list_dist(:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Copy diagnostics
diag_tmp = diag

!$omp parallel do schedule(static) private(ic2,list,list_dist,nc2eff,jc2,distnorm,norm,wgt)
do ic2=1,hdata%nc2
   if (hdata%ic1il0_log(hdata%ic2_to_ic1(ic2),il0)) then
      ! Allocation
      allocate(list(hdata%nc2))
      allocate(list_dist(hdata%nc2))

      ! Build list of valid points
      nc2eff = 0
      jc2 = 1
      do while (hdata%nn_nc2_dist(jc2,ic2,min(il0,geom%nl0i))<r)
         ! Check the point validity
         if (isnotmsr(diag_tmp(hdata%nn_nc2_index(jc2,ic2,min(il0,geom%nl0i))))) then
            nc2eff = nc2eff+1
            list(nc2eff) = diag_tmp(hdata%nn_nc2_index(jc2,ic2,min(il0,geom%nl0i)))
            list_dist(nc2eff) = hdata%nn_nc2_dist(jc2,ic2,min(il0,geom%nl0i))
         end if
         jc2 = jc2+1
         if (jc2>hdata%nc2) exit
      end do

      ! Apply filter
      if (nc2eff>0) then
         select case (trim(filter_type))
         case ('average')
            ! Compute average
            diag(ic2) = sum(list(1:nc2eff))/float(nc2eff)
         case ('gc99')
            ! Gaspari-Cohn (1999) kernel
            diag(ic2) = 0.0
            norm = 0.0
            do jc2=1,nc2eff
               distnorm = list_dist(jc2)/r
               wgt = gc99(distnorm)
               diag(ic2) = diag(ic2)+wgt*list(jc2)
               norm = norm+wgt
            end do
            diag(ic2) = diag(ic2)/norm
         case ('median')
            ! Compute median
            diag(ic2) = median(nc2eff,list(1:nc2eff))
         case default
            ! Wrong filter
            call msgerror('wrong filter type')
         end select
      else
         call msr(diag(ic2))
      end if

      ! Release memory
      deallocate(list)
      deallocate(list_dist)
   end if
end do
!$omp end parallel do

! End associate
end associate

end subroutine diag_filter_2d

!----------------------------------------------------------------------
! Subroutine: diag_filter_3d
!> Purpose: filter diagnostics, on multiple levels
!----------------------------------------------------------------------
subroutine diag_filter_3d(hdata,filter_type,r,diag)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                             !< HDIAG data
character(len=*),intent(in) :: filter_type                      !< Filter type
real(kind_real),intent(in) :: r                                 !< Filter support radius
real(kind_real),intent(inout) :: diag(hdata%nc2,hdata%geom%nl0) !< Filtered diagnostics

! Local variables
integer :: il0

! Associate
associate(geom=>hdata%geom)

do il0=1,geom%nl0
   ! Filter diagnostics on a single level
   call diag_filter_2d(hdata,il0,filter_type,r,diag(:,il0))
end do

! End associate
end associate

end subroutine diag_filter_3d

!----------------------------------------------------------------------
! Subroutine: diag_interpolation_2d
!> Purpose: compute interpolation from diagnostic points to grid-points, on a single level
!----------------------------------------------------------------------
subroutine diag_interpolation_2d(hdata,il0,fld_nc2,fld)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                !< HDIAG data
integer,intent(in) :: il0                          !< Level
real(kind_real),intent(in) :: fld_nc2(hdata%nc2)   !< Subgrid field
real(kind_real),intent(out) :: fld(hdata%geom%nc0) !< Field

! Associate
associate(geom=>hdata%geom)

! Apply interpolation
call apply_linop(hdata%h(min(il0,geom%nl0i)),fld_nc2,fld)

! End associate
end associate

end subroutine diag_interpolation_2d

!----------------------------------------------------------------------
! Subroutine: diag_interpolation_3d
!> Purpose: compute interpolation from diagnostic points to grid-points, on multiple levels
!----------------------------------------------------------------------
subroutine diag_interpolation_3d(hdata,fld_nc2,fld)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                               !< HDIAG data
real(kind_real),intent(in) :: fld_nc2(hdata%nc2,hdata%geom%nl0)   !< Subgrid field
real(kind_real),intent(out) :: fld(hdata%geom%nc0,hdata%geom%nl0) !< Field

! Local variables
integer :: il0

! Associate
associate(geom=>hdata%geom)

do il0=1,geom%nl0
   ! Compute interpolation on a single level
   call diag_interpolation_2d(hdata,il0,fld_nc2(:,il0),fld(:,il0))
end do

! End associate
end associate

end subroutine diag_interpolation_3d

!----------------------------------------------------------------------
! Subroutine: diag_write_2d
!> Purpose: write model field on sample grid, on a single level
!----------------------------------------------------------------------
subroutine diag_write_2d(filename,varname,hdata,fld)

implicit none

! Passed variables
character(len=*),intent(in) :: filename      !< File name
character(len=*),intent(in) :: varname       !< Variable name
type(hdatatype),intent(in) :: hdata          !< HDIAG data
real(kind_real),intent(in) :: fld(hdata%nc2) !< Written field

! Local variables
integer :: ierr
integer :: ncid,nc2_id,fld_id,lon_id,lat_id
character(len=1024) :: subr = 'diag_write_2d'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Check if the file exists
ierr = nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_noclobber,nf90_64bit_offset),ncid)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_write,ncid))
   call ncerr(subr,nf90_redef(ncid))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))
end if
call ncerr(subr,nf90_enddef(ncid))

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Create dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'nc2',nc2_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nc2',hdata%nc2,nc2_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nc2_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
call ncerr(subr,nf90_put_var(ncid,fld_id,fld))

! Add coordinates
ierr = nf90_inq_varid(ncid,'longitude',lon_id)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   call ncerr(subr,nf90_def_var(ncid,'longitude',ncfloat,(/nc2_id/),lon_id))
   call ncerr(subr,nf90_def_var(ncid,'latitude',ncfloat,(/nc2_id/),lat_id))
   call ncerr(subr,nf90_enddef(ncid))
   call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon(hdata%ic2_to_ic0)*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat(hdata%ic2_to_ic0)*rad2deg))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine diag_write_2d

!----------------------------------------------------------------------
! Subroutine: diag_write_3d
!> Purpose: write model field on sample grid, on multiple levels
!----------------------------------------------------------------------
subroutine diag_write_3d(filename,varname,hdata,fld)

implicit none

! Passed variables
character(len=*),intent(in) :: filename                     !< File name
character(len=*),intent(in) :: varname                      !< Variable name
type(hdatatype),intent(in) :: hdata                         !< HDIAG data
real(kind_real),intent(in) :: fld(hdata%nc2,hdata%geom%nl0) !< Written field

! Local variables
integer :: ierr
integer :: ncid,nl0_id,nc2_id,fld_id,lon_id,lat_id
character(len=1024) :: subr = 'diag_write_3d'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Check if the file exists
ierr = nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_noclobber,nf90_64bit_offset),ncid)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_write,ncid))
   call ncerr(subr,nf90_redef(ncid))
   call ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))
end if
call ncerr(subr,nf90_enddef(ncid))

! Get variable id
ierr = nf90_inq_varid(ncid,trim(varname),fld_id)

! Create dimensions and variable if necessary
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   ierr = nf90_inq_dimid(ncid,'nl0',nl0_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
   ierr = nf90_inq_dimid(ncid,'nc2',nc2_id)
   if (ierr/=nf90_noerr) call ncerr(subr,nf90_def_dim(ncid,'nc2',hdata%nc2,nc2_id))
   call ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nc2_id,nl0_id/),fld_id))
   call ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_enddef(ncid))
end if

! Write data
call ncerr(subr,nf90_put_var(ncid,fld_id,fld))

! Add coordinates
ierr = nf90_inq_varid(ncid,'longitude',lon_id)
if (ierr/=nf90_noerr) then
   call ncerr(subr,nf90_redef(ncid))
   call ncerr(subr,nf90_def_var(ncid,'longitude',ncfloat,(/nc2_id/),lon_id))
   call ncerr(subr,nf90_def_var(ncid,'latitude',ncfloat,(/nc2_id/),lat_id))
   call ncerr(subr,nf90_enddef(ncid))
   call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon(hdata%ic2_to_ic0)*rad2deg))
   call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat(hdata%ic2_to_ic0)*rad2deg))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine diag_write_3d

end module module_diag_tools
