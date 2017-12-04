!----------------------------------------------------------------------
! Module: type_curve
!> Purpose: curve derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_curve

use model_interface, only: model_write
use module_diag_tools, only: diag_write,diag_filter,diag_interpolation
use netcdf
use tools_display, only: vunitchar,msgerror
use tools_jacobi_eigenvalue, only: jacobi_eigenvalue
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isnotmsr,isallnotmsr
use tools_nc, only: ncerr,ncfloat
use type_hdata, only: hdatatype
use type_mpl, only: mpl
use type_nam, only: namtype,namncwrite

implicit none

! Curve derived type
type curvetype
   character(len=1024) :: cname                   !< Curve name
   integer :: npack                               !< Pack buffer size
   real(kind_real),allocatable :: raw(:,:,:)      !< Raw curve
   real(kind_real),allocatable :: raw_coef_ens(:) !< Raw ensemble coefficient
   real(kind_real) :: raw_coef_sta                !< Raw static coefficient
   real(kind_real),allocatable :: fit_wgt(:,:,:)  !< Fit weight
   real(kind_real),allocatable :: fit(:,:,:)      !< Fit
   real(kind_real),allocatable :: fit_rh(:)       !< Fit support radius
   real(kind_real),allocatable :: fit_rv(:)       !< Fit support radius
end type curvetype

private
public :: curvetype
public :: curve_alloc,curve_dealloc,curve_normalization
public :: curve_pack,curve_unpack,curve_write,curve_write_all,curve_write_local

contains

!----------------------------------------------------------------------
! Subroutine: curve_alloc
!> Purpose: curve object allocation
!----------------------------------------------------------------------
subroutine curve_alloc(hdata,ib,cname,curve)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata    !< HDIAG data
integer,intent(in) :: ib               !< Block index
character(len=*),intent(in) :: cname   !< Curve name
type(curvetype),intent(inout) :: curve !< Curve

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Set name
curve%cname = cname

! Allocation
allocate(curve%raw(nam%nc,bpar%nl0(ib),geom%nl0))
allocate(curve%raw_coef_ens(geom%nl0))

! Initialization
curve%npack = nam%nc*geom%nl0*bpar%nl0(ib)+2*geom%nl0
call msr(curve%raw)
call msr(curve%raw_coef_ens)

if (trim(nam%fit_type)/='none') then
   ! Allocation
   allocate(curve%fit_wgt(nam%nc,bpar%nl0(ib),geom%nl0))
   allocate(curve%fit(nam%nc,bpar%nl0(ib),geom%nl0))
   allocate(curve%fit_rh(geom%nl0))
   allocate(curve%fit_rv(geom%nl0))

   ! Initialization
   curve%fit_wgt = 1.0
   call msr(curve%fit)
   call msr(curve%fit_rh)
   call msr(curve%fit_rv)
end if

! End ssociate
end associate

end subroutine curve_alloc

!----------------------------------------------------------------------
! Subroutine: curve_dealloc
!> Purpose: curve object deallocation
!----------------------------------------------------------------------
subroutine curve_dealloc(hdata,curve)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata    !< HDIAG data
type(curvetype),intent(inout) :: curve !< Curve

! Associate
associate(nam=>hdata%nam)

! Deallocation
deallocate(curve%raw)
deallocate(curve%raw_coef_ens)
if (trim(nam%fit_type)/='none') then
   deallocate(curve%fit_wgt)
   deallocate(curve%fit)
   deallocate(curve%fit_rh)
   deallocate(curve%fit_rv)
end if

! End associate
end associate

end subroutine curve_dealloc

!----------------------------------------------------------------------
! Subroutine: curve_normalization
!> Purpose: compute localization normalization
!----------------------------------------------------------------------
subroutine curve_normalization(hdata,ib,curve)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata    !< HDIAG data
integer,intent(in) :: ib               !< Block index
type(curvetype),intent(inout) :: curve !< Curve

! Local variables
integer :: il0r,il0,jl0,ic

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Get diagonal values
do jl0=1,geom%nl0
   il0r = bpar%il0rz(jl0,ib)
   if (isnotmsr(curve%raw(1,il0r,jl0))) curve%raw_coef_ens(jl0) = curve%raw(1,il0r,jl0)
end do

! Normalize
do jl0=1,geom%nl0
   do il0r=1,bpar%nl0(ib)
      il0 = bpar%il0rjl0ib_to_il0(il0r,jl0,ib)
      do ic=1,bpar%icmax(ib)
         if (isnotmsr(curve%raw(ic,il0r,jl0)).and.isnotmsr(curve%raw_coef_ens(il0)).and.isnotmsr(curve%raw_coef_ens(jl0))) &
       & curve%raw(ic,il0r,jl0) = curve%raw(ic,il0r,jl0)/sqrt(curve%raw_coef_ens(il0)*curve%raw_coef_ens(jl0))
      end do
   end do
end do

! End associate
end associate

end subroutine curve_normalization

!----------------------------------------------------------------------
! Subroutine: curve_pack
!> Purpose: curve packing
!----------------------------------------------------------------------
subroutine curve_pack(hdata,curve,buf)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata             !< HDIAG data
type(curvetype),intent(in) :: curve             !< Curve
real(kind_real),intent(out) :: buf(curve%npack) !< Buffer

! Local variables
integer :: offset

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Pack
offset = 0
buf(offset+1:offset+nam%nc*nam%nl0r*geom%nl0) = pack(curve%fit,.true.)
offset = offset+nam%nc*nam%nl0r*geom%nl0
buf(offset+1:offset+geom%nl0) = curve%fit_rh
offset = offset+geom%nl0
buf(offset+1:offset+geom%nl0) = curve%fit_rv

! End associate
end associate

end subroutine curve_pack

!----------------------------------------------------------------------
! Subroutine: curve_unpack
!> Purpose: curve unpacking
!----------------------------------------------------------------------
subroutine curve_unpack(hdata,curve,buf)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata            !< HDIAG data
type(curvetype),intent(inout) :: curve         !< Curve
real(kind_real),intent(in) :: buf(curve%npack) !< Buffer

! Local variables
integer :: offset
logical,allocatable :: mask_unpack(:,:,:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Allocation
allocate(mask_unpack(nam%nc,nam%nl0r,geom%nl0))
mask_unpack = .true.

! Unpack
offset = 0
curve%fit = unpack(buf(offset+1:offset+nam%nc*nam%nl0r*geom%nl0),mask_unpack,curve%fit)
offset = offset+nam%nc*nam%nl0r*geom%nl0
curve%fit_rh = buf(offset+1:offset+geom%nl0)
offset = offset+geom%nl0
curve%fit_rv = buf(offset+1:offset+geom%nl0)

! End associate
end associate

end subroutine curve_unpack

!----------------------------------------------------------------------
! Subroutine: curve_write
!> Purpose: write a curve
!----------------------------------------------------------------------
subroutine curve_write(hdata,ncid,curve)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata    !< HDIAG data
integer,intent(in) :: ncid             !< NetCDF file id
type(curvetype),intent(inout) :: curve !< Curve

! Local variables
integer :: one_id,nc_id,nl0r_id,nl0_id
integer :: raw_id,raw_coef_ens_id,raw_coef_sta_id
integer :: fit_id,fit_rh_id,fit_rv_id
character(len=1024) :: subr = 'curve_write'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Get dimensions ID
call ncerr(subr,nf90_inq_dimid(ncid,'one',one_id))
call ncerr(subr,nf90_inq_dimid(ncid,'nc',nc_id))
call ncerr(subr,nf90_inq_dimid(ncid,'nl0r',nl0r_id))
call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))

! Define variables
call ncerr(subr,nf90_redef(ncid))

! Raw curve
call ncerr(subr,nf90_def_var(ncid,trim(curve%cname)//'_raw',ncfloat,(/nc_id,nl0r_id,nl0_id/),raw_id))
call ncerr(subr,nf90_put_att(ncid,raw_id,'_FillValue',msvalr))
! Raw curve ensemble coefficient
call ncerr(subr,nf90_def_var(ncid,trim(curve%cname)//'_raw_coef_ens',ncfloat,(/nl0_id/),raw_coef_ens_id))
call ncerr(subr,nf90_put_att(ncid,raw_coef_ens_id,'_FillValue',msvalr))
if (isnotmsr(curve%raw_coef_sta)) then
   ! Raw curve static coefficient
   call ncerr(subr,nf90_def_var(ncid,trim(curve%cname)//'_raw_coef_sta',ncfloat,(/one_id/),raw_coef_sta_id))
   call ncerr(subr,nf90_put_att(ncid,raw_coef_sta_id,'_FillValue',msvalr))
end if
if (trim(nam%fit_type)/='none') then
   ! Fitted curve
   call ncerr(subr,nf90_def_var(ncid,trim(curve%cname)//'_fit',ncfloat,(/nc_id,nl0r_id,nl0_id/),fit_id))
   call ncerr(subr,nf90_put_att(ncid,fit_id,'_FillValue',msvalr))
   ! Fitted curve support radius
   call ncerr(subr,nf90_def_var(ncid,trim(curve%cname)//'_fit_rh',ncfloat,(/nl0_id/),fit_rh_id))
   call ncerr(subr,nf90_put_att(ncid,fit_rh_id,'_FillValue',msvalr))
   ! Fitted curve support radius
   call ncerr(subr,nf90_def_var(ncid,trim(curve%cname)//'_fit_rv',ncfloat,(/nl0_id/),fit_rv_id))
   call ncerr(subr,nf90_put_att(ncid,fit_rv_id,'_FillValue',msvalr))
end if

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write variables

! Raw curve
call ncerr(subr,nf90_put_var(ncid,raw_id,curve%raw))
! Raw curve ensemble coefficient
call ncerr(subr,nf90_put_var(ncid,raw_coef_ens_id,curve%raw_coef_ens))
! Raw curve static coefficient
if (isnotmsr(curve%raw_coef_sta)) call ncerr(subr,nf90_put_var(ncid,raw_coef_sta_id,curve%raw_coef_sta))
if (trim(nam%fit_type)/='none') then
   ! Fitted curve
   call ncerr(subr,nf90_put_var(ncid,fit_id,curve%fit))
   ! Fitted curve support radius
   call ncerr(subr,nf90_put_var(ncid,fit_rh_id,curve%fit_rh))
   ! Fitted curve support radius
   call ncerr(subr,nf90_put_var(ncid,fit_rv_id,curve%fit_rv))
end if

! End associate
end associate

end subroutine curve_write

!----------------------------------------------------------------------
! Subroutine: curve_write_all
!> Purpose: write all curves
!----------------------------------------------------------------------
subroutine curve_write_all(hdata,filename,cor_1,cor_2,loc_1,loc_2,loc_3,loc_4)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                     !< HDIAG data
character(len=*),intent(in) :: filename                 !< File name
type(curvetype),intent(inout) :: cor_1(hdata%bpar%nb+1) !< Correlation 1
type(curvetype),intent(inout) :: cor_2(hdata%bpar%nb+1) !< Correlation 2
type(curvetype),intent(inout) :: loc_1(hdata%bpar%nb+1) !< Localization 1
type(curvetype),intent(inout) :: loc_2(hdata%bpar%nb+1) !< Localization 2
type(curvetype),intent(inout) :: loc_3(hdata%bpar%nb+1) !< Localization 3
type(curvetype),intent(inout) :: loc_4(hdata%bpar%nb+1) !< Localization 4

! Local variables
integer :: ncid,one_id,nc_id,nl0r_id,nl0_id,disth_id,vunit_id
integer :: ib
character(len=1024) :: subr = 'curve_write_all'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

call system('rm -f '//trim(nam%datadir)//'/'//trim(filename))
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))
call namncwrite(nam,ncid)
call ncerr(subr,nf90_put_att(ncid,nf90_global,'vunitchar',trim(vunitchar)))
call ncerr(subr,nf90_def_dim(ncid,'one',1,one_id))
call ncerr(subr,nf90_def_dim(ncid,'nc',nam%nc,nc_id))
call ncerr(subr,nf90_def_dim(ncid,'nl0r',nam%nl0r,nl0r_id))
call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
call ncerr(subr,nf90_def_var(ncid,'disth',ncfloat,(/nc_id/),disth_id))
call ncerr(subr,nf90_def_var(ncid,'vunit',ncfloat,(/nl0_id/),vunit_id))
call ncerr(subr,nf90_enddef(ncid))
call ncerr(subr,nf90_put_var(ncid,disth_id,geom%disth(1:nam%nc)))
call ncerr(subr,nf90_put_var(ncid,vunit_id,geom%vunit))
do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      call curve_write(hdata,ncid,cor_1(ib))
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd','dual-ens')
         call curve_write(hdata,ncid,cor_2(ib))
      end select
      select case (trim(nam%method))
      case ('loc','hyb-avg','hyb-rnd','dual-ens')
         call curve_write(hdata,ncid,loc_1(ib))
      end select
      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd','dual-ens')
         call curve_write(hdata,ncid,loc_2(ib))
      end select
      if (trim(nam%method)=='dual-ens') then
         call curve_write(hdata,ncid,loc_3(ib))
         call curve_write(hdata,ncid,loc_4(ib))
      end if
   end if
end do
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine curve_write_all

!----------------------------------------------------------------------
! Subroutine: curve_write_local
!> Purpose: write all curves
!----------------------------------------------------------------------
subroutine curve_write_local(hdata,filename,curve_nc2)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                                !< HDIAG data
character(len=*),intent(in) :: filename                            !< File name
type(curvetype),intent(in) :: curve_nc2(hdata%nc2,hdata%bpar%nb+1) !< Curves array

! Local variables
integer :: ncid
integer :: ic2,ib
real(kind_real) :: fld_nc2(hdata%nc2,hdata%geom%nl0),fld(hdata%geom%nc0,hdata%geom%nl0)
character(len=1024) :: subr = 'curve_write_all'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))
call namncwrite(nam,ncid)
call ncerr(subr,nf90_close(ncid))
do ib=1,bpar%nb+1
   if (bpar%fit_block(ib)) then
      call msr(fld_nc2)
      do ic2=1,hdata%nc2
         fld_nc2(ic2,:) = curve_nc2(ic2,ib)%fit_rh
      end do
      call diag_interpolation(hdata,fld_nc2,fld)
      call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rh',fld)
      if (trim(nam%flt_type)/='none') then
         call diag_filter(hdata,nam%flt_type,nam%diag_rhflt,fld_nc2)
         call diag_interpolation(hdata,fld_nc2,fld)
         call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rh_flt',fld)
      end if
      call msr(fld_nc2)
      do ic2=1,hdata%nc2
         fld_nc2(ic2,:) = curve_nc2(ic2,ib)%fit_rv
      end do
      call diag_interpolation(hdata,fld_nc2,fld)
      call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rv',fld)
      if (trim(nam%flt_type)/='none') then
         call diag_filter(hdata,nam%flt_type,nam%diag_rhflt,fld_nc2)
         call diag_interpolation(hdata,fld_nc2,fld)
         call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_fit_rv_flt',fld)
      end if
   end if
end do

! End associate
end associate

end subroutine curve_write_local

end module type_curve
