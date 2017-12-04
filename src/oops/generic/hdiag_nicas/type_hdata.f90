!----------------------------------------------------------------------
! Module: type_hdata
!> Purpose: sample data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_hdata

use netcdf
use tools_const, only: rad2deg
use tools_display, only: msgerror,msgwarning
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msi,msr
use tools_nc, only: ncerr,ncfloat
use type_bpar, only: bpartype
use type_ctree, only: ctreetype,delete_ctree
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_dealloc,linop_read,linop_write
use type_mpl, only: mpl
use type_nam, only: namtype

implicit none

! Sampling data derived type
type hdatatype
   ! Namelist
   type(namtype),pointer :: nam                      !< Namelist

   ! Geometry
   type(geomtype),pointer :: geom                    !< Geometry

   ! Block parameters
   type(bpartype),pointer :: bpar                    !< Block parameters

   ! Sampling
   integer,allocatable :: ic1_to_ic0(:)              !< First sampling index
   logical,allocatable :: ic1il0_log(:,:)            !< Log for the first sampling index
   integer,allocatable :: ic1icil0_to_ic0(:,:,:)     !< Second horizontal sampling index
   logical,allocatable :: ic1icil0_log(:,:,:)        !< Log for the second horizontal sampling index
   real(kind_real),allocatable :: swgt(:,:,:,:)      !< Sampling weights
   real(kind_real),allocatable :: bwgtsq(:,:,:,:)    !< Squared block weights

   integer :: nc2                                    !< Subgrid size
   integer,allocatable :: ic2_to_ic1(:)              !< Subgrid to diagnostic points
   integer,allocatable :: ic2_to_ic0(:)              !< Subgrid to grid

   ! Cover tree and nearest neighbors
   logical,allocatable ::  local_mask(:,:,:)         !< Local mask
   logical,allocatable ::  displ_mask(:,:,:)         !< Displacement mask
   integer,allocatable :: nn_nc2_index(:,:,:)        !< Nearest diagnostic neighbors from diagnostic points
   real(kind_real),allocatable :: nn_nc2_dist(:,:,:) !< Nearest diagnostic neighbors distance from diagnostic points
   integer,allocatable :: nn_ldwv_index(:)           !< Nearest diagnostic neighbors for local diagnostics profiles

   ! Sampling mesh
   integer :: nt                                     !< Number of triangles
   integer,allocatable :: ltri(:,:)                  !< Triangles indices
   integer,allocatable :: larc(:,:)                  !< Arcs indices
   real(kind_real),allocatable :: bdist(:)           !< Distance to the closest boundary arc

   ! Interpolations
   type(linoptype),allocatable :: h(:)               !< Horizontal interpolation from Sc2 to Sc0
   type(linoptype),allocatable :: s(:)               !< Horizontal interpolation from Sc2 to Sc1
end type hdatatype

private
public :: hdatatype
public :: hdata_alloc,hdata_dealloc,hdata_read,hdata_write

contains

!----------------------------------------------------------------------
! Subroutine: hdata_alloc
!> Purpose: hdata object allocation
!----------------------------------------------------------------------
subroutine hdata_alloc(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Allocation
allocate(hdata%ic1_to_ic0(nam%nc1))
allocate(hdata%ic1il0_log(nam%nc1,geom%nl0))
allocate(hdata%ic1icil0_to_ic0(nam%nc1,nam%nc,geom%nl0))
allocate(hdata%ic1icil0_log(nam%nc1,nam%nc,geom%nl0))
allocate(hdata%swgt(nam%nc1,nam%nc,geom%nl0,geom%nl0))
allocate(hdata%bwgtsq(nam%nc,geom%nl0,geom%nl0,bpar%nb))
if (nam%local_diag.or.nam%displ_diag) then
   allocate(hdata%ic2_to_ic1(hdata%nc2))
   allocate(hdata%ic2_to_ic0(hdata%nc2))
end if
if (nam%local_diag.or.nam%displ_diag) then
   allocate(hdata%local_mask(nam%nc1,hdata%nc2,geom%nl0i))
   allocate(hdata%displ_mask(nam%nc1,hdata%nc2,geom%nl0i))
   if (trim(nam%flt_type)/='none') then
      allocate(hdata%nn_nc2_index(hdata%nc2,hdata%nc2,geom%nl0i))
      allocate(hdata%nn_nc2_dist(hdata%nc2,hdata%nc2,geom%nl0i))
   end if
   allocate(hdata%h(geom%nl0i))
   allocate(hdata%s(geom%nl0i))
end if

! Initialization
call msi(hdata%ic1_to_ic0)
hdata%ic1il0_log = .false.
call msi(hdata%ic1icil0_to_ic0)
hdata%ic1icil0_log = .false.
call msr(hdata%swgt)
call msr(hdata%bwgtsq)
if (nam%local_diag.or.nam%displ_diag) then
   call msi(hdata%ic2_to_ic1)
   call msi(hdata%ic2_to_ic0)
end if
if (nam%local_diag.or.nam%displ_diag) then
   hdata%local_mask = .false.
   hdata%displ_mask = .false.
   if (trim(nam%flt_type)/='none') then
      call msi(hdata%nn_nc2_index)
      call msr(hdata%nn_nc2_dist)
   end if
end if

! End associate
end associate

end subroutine hdata_alloc

!----------------------------------------------------------------------
! Subroutine: hdata_dealloc
!> Purpose: hdata object deallocation
!----------------------------------------------------------------------
subroutine hdata_dealloc(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: il0

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Release memory
deallocate(hdata%ic1_to_ic0)
deallocate(hdata%ic1il0_log)
deallocate(hdata%ic1icil0_to_ic0)
deallocate(hdata%ic1icil0_log)
deallocate(hdata%swgt)
deallocate(hdata%bwgtsq)
if (nam%local_diag.or.nam%displ_diag) then
   deallocate(hdata%ic2_to_ic1)
   deallocate(hdata%ic2_to_ic0)
   deallocate(hdata%local_mask)
   deallocate(hdata%displ_mask)
   if (trim(nam%flt_type)/='none') then
      deallocate(hdata%nn_nc2_index)
      deallocate(hdata%nn_nc2_dist)
   end if
end if
if (nam%local_diag.or.nam%displ_diag) then
   do il0=1,geom%nl0
      call linop_dealloc(hdata%h(il0))
      call linop_dealloc(hdata%s(il0))
   end do
   deallocate(hdata%h)
   deallocate(hdata%s)
end if
if (nam%local_diag.and.(nam%nldwv>0)) deallocate(hdata%nn_ldwv_index)

! End associate
end associate

end subroutine hdata_dealloc

!----------------------------------------------------------------------
! Subroutine: hdata_read
!> Purpose: read hdata object
!----------------------------------------------------------------------
integer function hdata_read(hdata)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata !< HDIAG data

! Local variables
integer :: il0,il0i,ic1,ic,ic2
integer :: nl0_1_test,nl0_2_test,nc_test,nc1_test,nc2_test
integer :: info,ncid,nl0_1_id,nl0_2_id,nc_id,nc1_id,nc2_id
integer :: ic1_to_ic0_id,ic1il0_log_id,ic1icil0_to_ic0_id,ic1icil0_log_id,swgt_id
integer :: ic2_to_ic1_id,ic2_to_ic0_id,local_mask_id,displ_mask_id,nn_nc2_index_id,nn_nc2_dist_id
integer :: ic1il0_logint(hdata%nam%nc1,hdata%geom%nl0),ic1icil0_logint(hdata%nam%nc1,hdata%nam%nc,hdata%geom%nl0)
integer,allocatable :: local_maskint(:,:),displ_maskint(:,:)
character(len=3) :: il0ichar
character(len=1024) :: subr = 'hdata_read'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Initialization
hdata_read = 0

! Open file
info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',nf90_nowrite,ncid)
if (info/=nf90_noerr) then
   call msgwarning('cannot find sampling data to read, recomputing sampling')
   nam%sam_write = .true.
   hdata_read = 1
   return
end if

! Check dimensions
call ncerr(subr,nf90_inq_dimid(ncid,'nl0_1',nl0_1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl0_1_id,len=nl0_1_test))
call ncerr(subr,nf90_inq_dimid(ncid,'nl0_2',nl0_2_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl0_2_id,len=nl0_2_test))
call ncerr(subr,nf90_inq_dimid(ncid,'nc',nc_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc_id,len=nc_test))
call ncerr(subr,nf90_inq_dimid(ncid,'nc1',nc1_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nc1_id,len=nc1_test))
if (nam%local_diag.or.nam%displ_diag) then
   info = nf90_inq_dimid(ncid,'nc2_1',nc2_id)
   if (info==nf90_noerr) then
      call ncerr(subr,nf90_inquire_dimension(ncid,nc2_id,len=nc2_test))
   else
      call msgwarning('cannot find nc2 when reading sampling, recomputing sampling')
      nam%sam_write = .true.
      hdata_read = 2
   end if
end if
if ((geom%nl0/=nl0_1_test).or.(geom%nl0/=nl0_2_test).or.(nam%nc/=nc_test).or.(nam%nc1/=nc1_test)) then
   call msgwarning('wrong dimension when reading sampling, recomputing sampling')
   nam%sam_write = .true.
   call ncerr(subr,nf90_close(ncid))
   hdata_read = 1
   return
end if
if (nam%local_diag.or.nam%displ_diag) then
   if (hdata%nc2/=nc2_test) then
      call msgwarning('wrong dimension when reading sampling, recomputing sampling')
      nam%sam_write = .true.
      hdata_read = 2
   end if
end if

write(mpl%unit,'(a7,a)') '','Read sampling'

! Get arrays ID
call ncerr(subr,nf90_inq_varid(ncid,'ic1_to_ic0',ic1_to_ic0_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic1il0_log',ic1il0_log_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic1icil0_to_ic0',ic1icil0_to_ic0_id))
call ncerr(subr,nf90_inq_varid(ncid,'ic1icil0_log',ic1icil0_log_id))
call ncerr(subr,nf90_inq_varid(ncid,'swgt',swgt_id))
if ((hdata_read==0).and.(nam%local_diag.or.nam%displ_diag)) then
   call ncerr(subr,nf90_inq_varid(ncid,'ic2_to_ic1',ic2_to_ic1_id))
   call ncerr(subr,nf90_inq_varid(ncid,'ic2_to_ic0',ic2_to_ic0_id))
end if

! Read arrays
call ncerr(subr,nf90_get_var(ncid,ic1_to_ic0_id,hdata%ic1_to_ic0))
call ncerr(subr,nf90_get_var(ncid,ic1il0_log_id,ic1il0_logint))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (ic1il0_logint(ic1,il0)==0) then
         hdata%ic1il0_log(ic1,il0) = .false.
      else if (ic1il0_logint(ic1,il0)==1) then
         hdata%ic1il0_log(ic1,il0) = .true.
      end if
   end do
end do
call ncerr(subr,nf90_get_var(ncid,ic1icil0_to_ic0_id,hdata%ic1icil0_to_ic0))
call ncerr(subr,nf90_get_var(ncid,ic1icil0_log_id,ic1icil0_logint))
do il0=1,geom%nl0
   do ic=1,nam%nc
      do ic1=1,nam%nc1
         if (ic1icil0_logint(ic1,ic,il0)==0) then
            hdata%ic1icil0_log(ic1,ic,il0) = .false.
         else if (ic1icil0_logint(ic1,ic,il0)==1) then
            hdata%ic1icil0_log(ic1,ic,il0) = .true.
         end if
      end do
   end do
end do
call ncerr(subr,nf90_get_var(ncid,swgt_id,hdata%swgt))
if ((hdata_read==0).and.(nam%local_diag.or.nam%displ_diag)) then
   call ncerr(subr,nf90_get_var(ncid,ic2_to_ic1_id,hdata%ic2_to_ic1))
   call ncerr(subr,nf90_get_var(ncid,ic2_to_ic0_id,hdata%ic2_to_ic0))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

! Read nearest neighbors and interpolation
if ((hdata_read==0).and.(nam%local_diag.or.nam%displ_diag)) then
   ! Allocation
   allocate(local_maskint(nam%nc1,hdata%nc2))
   allocate(displ_maskint(nam%nc1,hdata%nc2))

   do il0i=1,geom%nl0i
      write(il0ichar,'(i3.3)') il0i
      info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc',nf90_nowrite,ncid)
      if (info/=nf90_noerr) then
         call msgwarning('cannot find nearest neighbors and interpolation data to read, recomputing sampling')
         nam%sam_write = .true.
         hdata_read = 3
         return
      end if
      call ncerr(subr,nf90_inq_varid(ncid,'local_mask',local_mask_id))
      call ncerr(subr,nf90_inq_varid(ncid,'displ_mask',displ_mask_id))
      call ncerr(subr,nf90_get_var(ncid,local_mask_id,local_maskint))
      call ncerr(subr,nf90_get_var(ncid,displ_mask_id,displ_maskint))
      do ic2=1,hdata%nc2
         do ic1=1,nam%nc1
            if (local_maskint(ic1,ic2)==1) then
               hdata%local_mask(ic1,ic2,il0i) = .true.
            elseif (local_maskint(ic1,ic2)==0) then
               hdata%local_mask(ic1,ic2,il0i) = .false.
            else
               call msgerror('wrong local_mask')
            end if
            if (displ_maskint(ic1,ic2)==1) then
               hdata%displ_mask(ic1,ic2,il0i) = .true.
            elseif (displ_maskint(ic1,ic2)==0) then
               hdata%displ_mask(ic1,ic2,il0i) = .false.
            else
               call msgerror('wrong displ_mask')
            end if
         end do
      end do
      if (trim(nam%flt_type)/='none') then
         info = nf90_inq_varid(ncid,'nn_nc2_index',nn_nc2_index_id)
         if (info==nf90_noerr) then
            call ncerr(subr,nf90_inq_varid(ncid,'nn_nc2_dist',nn_nc2_dist_id))
            call ncerr(subr,nf90_get_var(ncid,nn_nc2_index_id,hdata%nn_nc2_index(:,:,il0i)))
            call ncerr(subr,nf90_get_var(ncid,nn_nc2_dist_id,hdata%nn_nc2_dist(:,:,il0i)))
         else
            call msgwarning('cannot find nc2 nearest neighbors data to read, recomputing sampling')
            nam%sam_write = .true.
            hdata_read = 4
         end if
      end if
      call linop_read(ncid,'h',hdata%h(il0i))
      call linop_read(ncid,'s',hdata%s(il0i))
      call ncerr(subr,nf90_close(ncid))
   end do

   ! Release memory
   deallocate(local_maskint)
   deallocate(displ_maskint)
end if

! End associate
end associate

end function hdata_read

!----------------------------------------------------------------------
! Subroutine: hdata_write
!> Purpose: write hdata object
!----------------------------------------------------------------------
subroutine hdata_write(hdata)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata !< HDIAG data

! Local variables
integer :: il0,il0i,ic1,ic,ic2
integer :: ncid,nl0_1_id,nl0_2_id,nc1_id,nc2_1_id,nc2_2_id,nc_id
integer :: lat_id,lon_id,smax_id,ic1_to_ic0_id,ic1il0_log_id,ic1icil0_to_ic0_id,ic1icil0_log_id,swgt_id
integer :: ic2_to_ic1_id,ic2_to_ic0_id,local_mask_id,displ_mask_id,nn_nc2_index_id,nn_nc2_dist_id
integer :: ic1il0_logint(hdata%nam%nc1,hdata%geom%nl0),ic1icil0_logint(hdata%nam%nc1,hdata%nam%nc,hdata%geom%nl0)
integer,allocatable :: local_maskint(:,:),displ_maskint(:,:)
real(kind_real) :: lon(hdata%nam%nc1,hdata%nam%nc,hdata%geom%nl0),lat(hdata%nam%nc1,hdata%nam%nc,hdata%geom%nl0)
character(len=3) :: il0ichar
character(len=1024) :: subr = 'hdata_write'

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Create file
write(mpl%unit,'(a7,a)') '','Write sampling'
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(subr,nf90_def_dim(ncid,'nl0_1',geom%nl0,nl0_1_id))
call ncerr(subr,nf90_def_dim(ncid,'nl0_2',geom%nl0,nl0_2_id))
call ncerr(subr,nf90_def_dim(ncid,'nc',nam%nc,nc_id))
call ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
if (nam%local_diag.or.nam%displ_diag) call ncerr(subr,nf90_def_dim(ncid,'nc2_1',hdata%nc2,nc2_1_id))

! Define arrays
call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc1_id,nc_id,nl0_1_id/),lat_id))
call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc1_id,nc_id,nl0_1_id/),lon_id))
call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'smax',ncfloat,(/nc_id,nl0_1_id/),smax_id))
call ncerr(subr,nf90_put_att(ncid,smax_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'ic1_to_ic0',nf90_int,(/nc1_id/),ic1_to_ic0_id))
call ncerr(subr,nf90_put_att(ncid,ic1_to_ic0_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'ic1il0_log',nf90_int,(/nc1_id,nl0_1_id/),ic1il0_log_id))
call ncerr(subr,nf90_put_att(ncid,ic1il0_log_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'ic1icil0_to_ic0',nf90_int,(/nc1_id,nc_id,nl0_1_id/),ic1icil0_to_ic0_id))
call ncerr(subr,nf90_put_att(ncid,ic1icil0_to_ic0_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'ic1icil0_log',nf90_int,(/nc1_id,nc_id,nl0_1_id/),ic1icil0_log_id))
call ncerr(subr,nf90_put_att(ncid,ic1icil0_log_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'swgt',ncfloat,(/nc1_id,nc_id,nl0_1_id,nl0_2_id/),swgt_id))
call ncerr(subr,nf90_put_att(ncid,swgt_id,'_FillValue',msvalr))
if (nam%local_diag.or.nam%displ_diag) then
   call ncerr(subr,nf90_def_var(ncid,'ic2_to_ic1',nf90_int,(/nc2_1_id/),ic2_to_ic1_id))
   call ncerr(subr,nf90_put_att(ncid,ic2_to_ic1_id,'_FillValue',msvali))
   call ncerr(subr,nf90_def_var(ncid,'ic2_to_ic0',nf90_int,(/nc2_1_id/),ic2_to_ic0_id))
   call ncerr(subr,nf90_put_att(ncid,ic2_to_ic0_id,'_FillValue',msvali))
end if

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write arrays
call msr(lon)
call msr(lat)
do il0=1,geom%nl0
   do ic=1,nam%nc
      do ic1=1,nam%nc1
         if (hdata%ic1icil0_log(ic1,ic,il0)) then
            lon(ic1,ic,il0) = hdata%geom%lon(hdata%ic1icil0_to_ic0(ic1,ic,il0))*rad2deg
            lat(ic1,ic,il0) = hdata%geom%lat(hdata%ic1icil0_to_ic0(ic1,ic,il0))*rad2deg
         end if
      end do
   end do
end do
call ncerr(subr,nf90_put_var(ncid,lon_id,lon))
call ncerr(subr,nf90_put_var(ncid,lat_id,lat))
call ncerr(subr,nf90_put_var(ncid,smax_id,float(count(hdata%ic1icil0_log,dim=1))))
call ncerr(subr,nf90_put_var(ncid,ic1_to_ic0_id,hdata%ic1_to_ic0))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (hdata%ic1il0_log(ic1,il0)) then
         ic1il0_logint(ic1,il0) = 1
      else
         ic1il0_logint(ic1,il0) = 0
      end if
   end do
end do
call ncerr(subr,nf90_put_var(ncid,ic1il0_log_id,ic1il0_logint))
call ncerr(subr,nf90_put_var(ncid,ic1icil0_to_ic0_id,hdata%ic1icil0_to_ic0))
do il0=1,geom%nl0
   do ic=1,nam%nc
      do ic1=1,nam%nc1
         if (hdata%ic1icil0_log(ic1,ic,il0)) then
            ic1icil0_logint(ic1,ic,il0) = 1
         else
            ic1icil0_logint(ic1,ic,il0) = 0
         end if
      end do
   end do
end do
call ncerr(subr,nf90_put_var(ncid,ic1icil0_log_id,ic1icil0_logint))
call ncerr(subr,nf90_put_var(ncid,swgt_id,hdata%swgt))
if (nam%local_diag.or.nam%displ_diag) then
   call ncerr(subr,nf90_put_var(ncid,ic2_to_ic1_id,hdata%ic2_to_ic1))
   call ncerr(subr,nf90_put_var(ncid,ic2_to_ic0_id,hdata%ic2_to_ic0))
end if

! Close file
call ncerr(subr,nf90_close(ncid))

! Write nearest neighbors and interpolation
if (nam%local_diag.or.nam%displ_diag) then
   ! Allocation
   allocate(local_maskint(nam%nc1,hdata%nc2))
   allocate(displ_maskint(nam%nc1,hdata%nc2))

   do il0i=1,geom%nl0i
      write(il0ichar,'(i3.3)') il0i
      call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_sampling_'//il0ichar//'.nc', &
    & or(nf90_clobber,nf90_64bit_offset),ncid))
      call ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
      call ncerr(subr,nf90_def_dim(ncid,'nc2_1',hdata%nc2,nc2_1_id))
      call ncerr(subr,nf90_def_var(ncid,'local_mask',nf90_int,(/nc1_id,nc2_1_id/),local_mask_id))
      call ncerr(subr,nf90_def_var(ncid,'displ_mask',nf90_int,(/nc1_id,nc2_1_id/),displ_mask_id))
      call ncerr(subr,nf90_put_att(ncid,local_mask_id,'_FillValue',msvali))
      call ncerr(subr,nf90_put_att(ncid,displ_mask_id,'_FillValue',msvali))
      if (trim(nam%flt_type)/='none') then
         call ncerr(subr,nf90_def_dim(ncid,'nc2_2',hdata%nc2,nc2_2_id))
         call ncerr(subr,nf90_def_var(ncid,'nn_nc2_index',nf90_int,(/nc2_1_id,nc2_2_id/),nn_nc2_index_id))
         call ncerr(subr,nf90_put_att(ncid,nn_nc2_index_id,'_FillValue',msvali))
         call ncerr(subr,nf90_def_var(ncid,'nn_nc2_dist',ncfloat,(/nc2_1_id,nc2_2_id/),nn_nc2_dist_id))
         call ncerr(subr,nf90_put_att(ncid,nn_nc2_dist_id,'_FillValue',msvalr))
      end if
      call ncerr(subr,nf90_enddef(ncid))
      do ic2=1,hdata%nc2
         do ic1=1,nam%nc1
            if (hdata%local_mask(ic1,ic2,il0i)) then
               local_maskint(ic1,ic2) = 1
            else
               local_maskint(ic1,ic2) = 0
            end if
            if (hdata%displ_mask(ic1,ic2,il0i)) then
               displ_maskint(ic1,ic2) = 1
            else
               displ_maskint(ic1,ic2) = 0
            end if
         end do
      end do
      call ncerr(subr,nf90_put_var(ncid,local_mask_id,local_maskint))
      call ncerr(subr,nf90_put_var(ncid,displ_mask_id,displ_maskint))
      if (trim(nam%flt_type)/='none') then
         call ncerr(subr,nf90_put_var(ncid,nn_nc2_index_id,hdata%nn_nc2_index(:,:,il0i)))
         call ncerr(subr,nf90_put_var(ncid,nn_nc2_dist_id,hdata%nn_nc2_dist(:,:,il0i)))
      end if
      call linop_write(ncid,hdata%h(il0i))
      call linop_write(ncid,hdata%s(il0i))
      call ncerr(subr,nf90_close(ncid))
   end do

   ! Release memory
   deallocate(local_maskint)
   deallocate(displ_maskint)
end if

! End associate
end associate

end subroutine hdata_write

end module type_hdata
