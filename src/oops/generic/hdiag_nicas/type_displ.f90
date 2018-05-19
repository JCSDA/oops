!----------------------------------------------------------------------
! Module: type_displ
!> Purpose: displacement data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_displ

use netcdf
use tools_const, only: rad2deg,reqkm
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr,isnotmsr
use tools_nc, only: ncfloat,ncerr
use type_hdata, only: hdatatype
use type_linop, only: linoptype,linop_dealloc
use type_mpl, only: mpl
use type_nam, only: namncwrite

implicit none

! Displacement data derived type
type displtype
   integer :: niter                                  !< Number of stored iterations
   real(kind_real),allocatable :: dist(:,:,:)        !< Displacement distance
   real(kind_real),allocatable :: valid(:,:,:)       !< Displacement validity
   real(kind_real),allocatable :: rhflt(:,:,:)       !< Displacement filtering support radius

   real(kind_real),allocatable :: lon_c2(:,:)        !< Longitude origin
   real(kind_real),allocatable :: lat_c2(:,:)        !< Latitude origin
   real(kind_real),allocatable :: lon_c2_raw(:,:,:)  !< Raw displaced longitude
   real(kind_real),allocatable :: lat_c2_raw(:,:,:)  !< Raw displaced latitude
   real(kind_real),allocatable :: dist_c2_raw(:,:,:) !< Raw displacement distance
   real(kind_real),allocatable :: lon_c2_flt(:,:,:)  !< Filtered displaced longitude
   real(kind_real),allocatable :: lat_c2_flt(:,:,:)  !< Filtered displaced latitude
   real(kind_real),allocatable :: dist_c2_flt(:,:,:) !< Displacement distance, filtered
   real(kind_real),allocatable :: lon_c0_flt(:,:,:)  !< Interpolated displaced longitude
   real(kind_real),allocatable :: lat_c0_flt(:,:,:)  !< Interpolated displaced latitude
end type displtype

private
public :: displtype
public :: displ_alloc,displ_dealloc,displ_write

contains

!----------------------------------------------------------------------
! Subroutine: displ_alloc
!> Purpose: displacement data allocation
!----------------------------------------------------------------------
subroutine displ_alloc(hdata,displ)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata    !< HDIAG data
type(displtype),intent(inout) :: displ !< Displacement data

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Allocation
allocate(displ%dist(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%valid(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%rhflt(0:nam%displ_niter,geom%nl0,2:nam%nts))
allocate(displ%lon_c2(hdata%nc2,geom%nl0))
allocate(displ%lat_c2(hdata%nc2,geom%nl0))
allocate(displ%lon_c2_raw(hdata%nc2,geom%nl0,2:nam%nts))
allocate(displ%lat_c2_raw(hdata%nc2,geom%nl0,2:nam%nts))
allocate(displ%dist_c2_raw(hdata%nc2,geom%nl0,2:nam%nts))
allocate(displ%lon_c2_flt(hdata%nc2,geom%nl0,2:nam%nts))
allocate(displ%lat_c2_flt(hdata%nc2,geom%nl0,2:nam%nts))
allocate(displ%dist_c2_flt(hdata%nc2,geom%nl0,2:nam%nts))
allocate(displ%lon_c0_flt(geom%nc0,geom%nl0,nam%nts))
allocate(displ%lat_c0_flt(geom%nc0,geom%nl0,nam%nts))

! Initialization
call msr(displ%dist)
call msr(displ%valid)
call msr(displ%rhflt)
call msr(displ%lon_c2)
call msr(displ%lat_c2)
call msr(displ%lon_c2_raw)
call msr(displ%lat_c2_raw)
call msr(displ%dist_c2_raw)
call msr(displ%lon_c2_flt)
call msr(displ%lat_c2_flt)
call msr(displ%dist_c2_flt)
call msr(displ%lon_c0_flt)
call msr(displ%lat_c0_flt)

! End associate
end associate

end subroutine displ_alloc

!----------------------------------------------------------------------
! Subroutine: displ_dealloc
!> Purpose: displacement data deallocation
!----------------------------------------------------------------------
subroutine displ_dealloc(displ)

implicit none

! Passed variables
type(displtype),intent(inout) :: displ !< Displacement data

! Deallocation
if (allocated(displ%dist)) deallocate(displ%dist)
if (allocated(displ%valid)) deallocate(displ%valid)
if (allocated(displ%rhflt)) deallocate(displ%rhflt)
if (allocated(displ%lon_c2)) deallocate(displ%lon_c2)
if (allocated(displ%lat_c2)) deallocate(displ%lat_c2)
if (allocated(displ%lon_c2_raw)) deallocate(displ%lon_c2_raw)
if (allocated(displ%lat_c2_raw)) deallocate(displ%lat_c2_raw)
if (allocated(displ%dist_c2_raw)) deallocate(displ%dist_c2_raw)
if (allocated(displ%lon_c2_flt)) deallocate(displ%lon_c2_flt)
if (allocated(displ%lat_c2_flt)) deallocate(displ%lat_c2_flt)
if (allocated(displ%dist_c2_flt)) deallocate(displ%dist_c2_flt)
if (allocated(displ%lon_c0_flt)) deallocate(displ%lon_c0_flt)
if (allocated(displ%lat_c0_flt)) deallocate(displ%lat_c0_flt)

end subroutine displ_dealloc

!----------------------------------------------------------------------
! Subroutine: displ_write
!> Purpose: write displacement data
!----------------------------------------------------------------------
subroutine displ_write(hdata,filename,displ)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata     !< HDIAG data
character(len=*),intent(in) :: filename !< File name
type(displtype),intent(in) :: displ     !< Displacement data

! Local variables
integer :: ncid,nc0_id,nc2_id,nl0_id,nts_id,displ_niter_id,vunit_id,valid_id,dist_id,rhflt_id
integer :: lon_c2_id,lat_c2_id,lon_c2_raw_id,lat_c2_raw_id,dist_c2_raw_id,lon_c2_flt_id,lat_c2_flt_id,dist_c2_flt_id
integer :: lon_c0_id,lat_c0_id,lon_c0_flt_id,lat_c0_flt_id
character(len=1024) :: subr = 'displ_write'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Create file
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(filename),or(nf90_clobber,nf90_64bit_offset),ncid))
call namncwrite(nam,ncid)

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
call ncerr(subr,nf90_def_dim(ncid,'nc2',hdata%nc2,nc2_id))
call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
call ncerr(subr,nf90_def_dim(ncid,'nts',nam%nts-1,nts_id))
call ncerr(subr,nf90_def_dim(ncid,'niter',nam%displ_niter+1,displ_niter_id))

! Define arrays
call ncerr(subr,nf90_def_var(ncid,'vunit',ncfloat,(/nl0_id/),vunit_id))
call ncerr(subr,nf90_def_var(ncid,'valid',ncfloat,(/displ_niter_id,nl0_id,nts_id/),valid_id))
call ncerr(subr,nf90_put_att(ncid,valid_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'dist',ncfloat,(/displ_niter_id,nl0_id,nts_id/),dist_id))
call ncerr(subr,nf90_put_att(ncid,dist_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'rhflt',ncfloat,(/displ_niter_id,nl0_id,nts_id/),rhflt_id))
call ncerr(subr,nf90_put_att(ncid,rhflt_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon_c2',ncfloat,(/nc2_id,nl0_id/),lon_c2_id))
call ncerr(subr,nf90_put_att(ncid,lon_c2_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lat_c2',ncfloat,(/nc2_id,nl0_id/),lat_c2_id))
call ncerr(subr,nf90_put_att(ncid,lat_c2_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),lon_c2_raw_id))
call ncerr(subr,nf90_put_att(ncid,lon_c2_raw_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lat_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),lat_c2_raw_id))
call ncerr(subr,nf90_put_att(ncid,lat_c2_raw_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'dist_c2_raw',ncfloat,(/nc2_id,nl0_id,nts_id/),dist_c2_raw_id))
call ncerr(subr,nf90_put_att(ncid,dist_c2_raw_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),lon_c2_flt_id))
call ncerr(subr,nf90_put_att(ncid,lon_c2_flt_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lat_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),lat_c2_flt_id))
call ncerr(subr,nf90_put_att(ncid,lat_c2_flt_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'dist_c2_flt',ncfloat,(/nc2_id,nl0_id,nts_id/),dist_c2_flt_id))
call ncerr(subr,nf90_put_att(ncid,dist_c2_flt_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon_c0',ncfloat,(/nc0_id/),lon_c0_id))
call ncerr(subr,nf90_put_att(ncid,lon_c0_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lat_c0',ncfloat,(/nc0_id/),lat_c0_id))
call ncerr(subr,nf90_put_att(ncid,lat_c0_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon_c0_flt',ncfloat,(/nc0_id,nl0_id,nts_id/),lon_c0_flt_id))
call ncerr(subr,nf90_put_att(ncid,lon_c0_flt_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lat_c0_flt',ncfloat,(/nc0_id,nl0_id,nts_id/),lat_c0_flt_id))
call ncerr(subr,nf90_put_att(ncid,lat_c0_flt_id,'_FillValue',msvalr))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write arrays
call ncerr(subr,nf90_put_var(ncid,vunit_id,geom%vunit))
call ncerr(subr,nf90_put_var(ncid,valid_id,displ%valid))
call ncerr(subr,nf90_put_var(ncid,dist_id,displ%dist*reqkm))
call ncerr(subr,nf90_put_var(ncid,rhflt_id,displ%rhflt*reqkm))
call ncerr(subr,nf90_put_var(ncid,lon_c2_id,displ%lon_c2*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_c2_id,displ%lat_c2*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lon_c2_raw_id,displ%lon_c2_raw(:,:,2:nam%nts)*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_c2_raw_id,displ%lat_c2_raw(:,:,2:nam%nts)*rad2deg))
call ncerr(subr,nf90_put_var(ncid,dist_c2_raw_id,displ%dist_c2_raw(:,:,2:nam%nts)*reqkm))
call ncerr(subr,nf90_put_var(ncid,lon_c2_flt_id,displ%lon_c2_flt(:,:,2:nam%nts)*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_c2_flt_id,displ%lat_c2_flt(:,:,2:nam%nts)*rad2deg))
call ncerr(subr,nf90_put_var(ncid,dist_c2_flt_id,displ%dist_c2_flt(:,:,2:nam%nts)*reqkm))
call ncerr(subr,nf90_put_var(ncid,lon_c0_id,geom%lon*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_c0_id,geom%lat*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lon_c0_flt_id,displ%lon_c0_flt(:,:,2:nam%nts)*rad2deg))
call ncerr(subr,nf90_put_var(ncid,lat_c0_flt_id,displ%lat_c0_flt(:,:,2:nam%nts)*rad2deg))

! Close file
call ncerr(subr,nf90_close(ncid))

! End associate
end associate

end subroutine displ_write

end module type_displ
