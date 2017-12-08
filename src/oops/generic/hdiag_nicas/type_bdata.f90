!----------------------------------------------------------------------
! Module: type_bdata
!> Purpose: sample data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_bdata

use module_diag_tools, only: diag_filter,diag_interpolation
use netcdf
use tools_display, only: msgwarning,msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isnotmsr,isallnotmsr
use tools_nc, only: ncerr,ncfloat
use type_curve, only: curvetype
use type_geom, only: geomtype
use type_hdata, only: hdatatype
use type_mpl, only: mpl
use type_nam, only: namtype

implicit none

! B data derived type
type bdatatype
   ! Block name
   character(len=1024) :: cname                 !< Block name

   ! Namelist
   type(namtype),pointer :: nam                 !< Namelist

   ! Geometry
   type(geomtype),pointer :: geom               !< Geometry

   ! Data
   real(kind_real),allocatable :: coef_ens(:,:) !< Ensemble coefficient
   real(kind_real),allocatable :: rh0(:,:)      !< Fit support radius
   real(kind_real),allocatable :: rv0(:,:)      !< Fit support radius
   real(kind_real),allocatable :: coef_sta(:,:) !< Static coefficient
   real(kind_real) :: wgt                       !< Block weight
end type bdatatype

interface diag_to_bdata
  module procedure diag_to_bdata
  module procedure diag_nc2_to_bdata
end interface

private
public :: bdatatype
public :: bdata_alloc,bdata_dealloc,diag_to_bdata,bdata_read,bdata_write

contains

!----------------------------------------------------------------------
! Subroutine: bdata_alloc
!> Purpose: bdata object allocation
!----------------------------------------------------------------------
subroutine bdata_alloc(bdata)

implicit none

! Passed variables
type(bdatatype),intent(inout) :: bdata !< Sampling data

! Associate
associate(nam=>bdata%nam,geom=>bdata%geom)

! Allocation
allocate(bdata%coef_ens(geom%nc0,geom%nl0))
allocate(bdata%rh0(geom%nc0,geom%nl0))
allocate(bdata%rv0(geom%nc0,geom%nl0))
allocate(bdata%coef_sta(geom%nc0,geom%nl0))

! Initialization
call msr(bdata%coef_ens)
call msr(bdata%rh0)
call msr(bdata%rv0)
call msr(bdata%coef_sta)
call msr(bdata%wgt)

! End associate
end associate

end subroutine bdata_alloc

!----------------------------------------------------------------------
! Subroutine: bdata_dealloc
!> Purpose: bdata object deallocation
!----------------------------------------------------------------------
subroutine bdata_dealloc(bdata)

implicit none

! Passed variables
type(bdatatype),intent(inout) :: bdata !< Sampling data

! Release memory
deallocate(bdata%coef_ens)
deallocate(bdata%rh0)
deallocate(bdata%rv0)
deallocate(bdata%coef_sta)

end subroutine bdata_dealloc

!----------------------------------------------------------------------
! Subroutine: diag_to_bdata
!> Purpose: copy diagnostics into bdata object
!----------------------------------------------------------------------
subroutine diag_to_bdata(hdata,ib,diag,bdata)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata
integer,intent(in) :: ib
type(curvetype),intent(in) :: diag
type(bdatatype),intent(inout) :: bdata !< Sampling data

! Local variables
integer :: il0

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

if (bpar%nicas_block(ib)) then
   do il0=1,geom%nl0
      bdata%coef_ens(:,il0) = diag%raw_coef_ens(il0)
      bdata%rh0(:,il0) = diag%fit_rh(il0)
      bdata%rv0(:,il0) = diag%fit_rv(il0)
      bdata%coef_sta(:,il0) = diag%raw_coef_sta
   end do
   bdata%wgt = sum(diag%raw_coef_ens)/float(geom%nl0)
else
   bdata%wgt = sum(diag%raw_coef_ens)/float(geom%nl0)
end if

! End associate
end associate

end subroutine diag_to_bdata

!----------------------------------------------------------------------
! Subroutine: diag_nc2_to_bdata
!> Purpose: copy local diagnostics into bdata object
!----------------------------------------------------------------------
subroutine diag_nc2_to_bdata(hdata,ib,diag,bdata)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata           !< HDIAG data
integer,intent(in) :: ib                      !< Block index
type(curvetype),intent(in) :: diag(hdata%nc2) !< Diagnostic curves
type(bdatatype),intent(inout) :: bdata        !< B data

! Local variables
integer :: ic2
real(kind_real) :: fld_nc2(hdata%nc2,bdata%geom%nl0)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

if (bpar%nicas_block(ib)) then
   do ic2=1,hdata%nc2
      fld_nc2(ic2,:) = diag(ic2)%raw_coef_ens
   end do
   call diag_filter(hdata,nam%flt_type,nam%diag_rhflt,fld_nc2)
   call diag_interpolation(hdata,fld_nc2,bdata%coef_ens)
   do ic2=1,hdata%nc2
      fld_nc2(ic2,:) = diag(ic2)%fit_rh
   end do
   call diag_filter(hdata,nam%flt_type,nam%diag_rhflt,fld_nc2)
   call diag_interpolation(hdata,fld_nc2,bdata%rh0)
   do ic2=1,hdata%nc2
      fld_nc2(ic2,:) = diag(ic2)%fit_rv
   end do
   call diag_filter(hdata,nam%flt_type,nam%diag_rhflt,fld_nc2)
   call diag_interpolation(hdata,fld_nc2,bdata%rv0)
   do ic2=1,hdata%nc2
      fld_nc2(ic2,:) = diag(ic2)%raw_coef_sta
   end do
   call diag_filter(hdata,nam%flt_type,nam%diag_rhflt,fld_nc2)
   call diag_interpolation(hdata,fld_nc2,bdata%coef_sta)
   do ic2=1,hdata%nc2
     fld_nc2(ic2,:) = diag(ic2)%raw_coef_ens
   end do
   bdata%wgt = sum(fld_nc2)/float(hdata%nc2*geom%nl0)
else
   do ic2=1,hdata%nc2
      fld_nc2(ic2,:) = diag(ic2)%raw_coef_ens
   end do
   bdata%wgt = sum(fld_nc2)/float(hdata%nc2*geom%nl0)
end if

! End associate
end associate

end subroutine diag_nc2_to_bdata

!----------------------------------------------------------------------
! Subroutine: bdata_read
!> Purpose: read bdata object
!----------------------------------------------------------------------
subroutine bdata_read(bdata)

implicit none

! Passed variables
type(bdatatype),intent(inout) :: bdata !< B data

! Local variables
integer :: nc0_test,nl0_test,il0
integer :: info,ncid,nc0_id,nl0_id
integer :: coef_ens_id,rh0_id,rv0_id,coef_sta_id
character(len=1024) :: subr = 'bdata_read'

! Associate
associate(nam=>bdata%nam,geom=>bdata%geom)

! Open file
info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_'//trim(bdata%cname)//'.nc',nf90_nowrite,ncid)
if (info==nf90_noerr) then
   ! Check dimensions
   call ncerr(subr,nf90_inq_dimid(ncid,'nc0',nc0_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,nc0_id,len=nc0_test))
   call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
   call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=nl0_test))
   if ((geom%nc0/=nc0_test).or.(geom%nl0/=nl0_test)) call msgerror('wrong dimension when reading B')

   ! Get arrays ID
   call ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))
   call ncerr(subr,nf90_inq_varid(ncid,'rh0',rh0_id))
   call ncerr(subr,nf90_inq_varid(ncid,'rv0',rv0_id))
   call ncerr(subr,nf90_inq_varid(ncid,'coef_sta',coef_sta_id))

   ! Read arrays
   call ncerr(subr,nf90_get_var(ncid,coef_ens_id,bdata%coef_ens))
   call ncerr(subr,nf90_get_var(ncid,rh0_id,bdata%rh0))
   call ncerr(subr,nf90_get_var(ncid,rv0_id,bdata%rv0))
   call ncerr(subr,nf90_get_var(ncid,coef_sta_id,bdata%coef_sta))

   ! Get main weight
   call ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',bdata%wgt))

   ! Close file
   call ncerr(subr,nf90_close(ncid))
else
   call msgwarning('cannot find B data to read, use namelist values')
   bdata%coef_ens = 1.0
   do il0=1,geom%nl0
      bdata%rh0(:,il0) = nam%rh(il0)
      bdata%rv0(:,il0) = nam%rv(il0)
   end do
   bdata%coef_sta = 0.0
   bdata%wgt = 1.0
end if

! Check
if (any((bdata%rh0<0.0).and.isnotmsr(bdata%rh0))) call msgerror('rh0 should be positive')
if (any((bdata%rv0<0.0).and.isnotmsr(bdata%rv0))) call msgerror('rv0 should be positive')

! End associate
end associate

end subroutine bdata_read

!----------------------------------------------------------------------
! Subroutine: bdata_write
!> Purpose: write bdata object
!----------------------------------------------------------------------
subroutine bdata_write(bdata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata !< B data

! Local variables
integer :: ncid,nc0_id,nl0_id
integer :: coef_ens_id,rh0_id,rv0_id,coef_sta_id
character(len=1024) :: subr = 'bdata_write'

! Associate
associate(nam=>bdata%nam,geom=>bdata%geom)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Create file
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_'//trim(bdata%cname)//'.nc', &
 & or(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))

! Define arrays
call ncerr(subr,nf90_def_var(ncid,'coef_ens',ncfloat,(/nc0_id,nl0_id/),coef_ens_id))
call ncerr(subr,nf90_put_att(ncid,coef_ens_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'rh0',ncfloat,(/nc0_id,nl0_id/),rh0_id))
call ncerr(subr,nf90_put_att(ncid,rh0_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'rv0',ncfloat,(/nc0_id,nl0_id/),rv0_id))
call ncerr(subr,nf90_put_att(ncid,rv0_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'coef_sta',ncfloat,(/nc0_id,nl0_id/),coef_sta_id))
call ncerr(subr,nf90_put_att(ncid,coef_sta_id,'_FillValue',msvalr))

! Write main weight
call ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',bdata%wgt))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write arrays
call ncerr(subr,nf90_put_var(ncid,coef_ens_id,bdata%coef_ens))
call ncerr(subr,nf90_put_var(ncid,rh0_id,bdata%rh0))
call ncerr(subr,nf90_put_var(ncid,rv0_id,bdata%rv0))
call ncerr(subr,nf90_put_var(ncid,coef_sta_id,bdata%coef_sta))

! Close file
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine bdata_write

end module type_bdata
