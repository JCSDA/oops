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
use type_com, only: comtype
use type_geom, only: geomtype
use type_linop, only: linoptype,linop_dealloc,linop_read,linop_write
use type_mpl, only: mpl
use type_nam, only: namtype

implicit none

! Sampling data derived type
type hdatatype
   ! Namelist
   type(namtype),pointer :: nam                     !< Namelist

   ! Geometry
   type(geomtype),pointer :: geom                   !< Geometry

   ! Block parameters
   type(bpartype),pointer :: bpar                   !< Block parameters

   ! Sampling
   integer,allocatable :: c1_to_c0(:)               !< First sampling index
   logical,allocatable :: c1l0_log(:,:)             !< Log for the first sampling index
   integer,allocatable :: c1c3_to_c0(:,:)           !< Second horizontal sampling index
   logical,allocatable :: c1c3l0_log(:,:,:)         !< Log for the second horizontal sampling index
   integer :: nc2                                   !< Subgrid size
   integer,allocatable :: c2_to_c1(:)               !< Subgrid to diagnostic points
   integer,allocatable :: c2_to_c0(:)               !< Subgrid to grid

   ! MPI splitting
   integer :: nc2a                                  !< Number of points in subset Sc2, halo A
   integer,allocatable :: c2a_to_c2(:)              !< Subset Sc2, halo A to global
   integer,allocatable :: c2_to_c2a(:)              !< Subset Sc2, global to halo A
   integer,allocatable :: c2_to_proc(:)             !< Subset Sc2, global to processor
   integer,allocatable :: proc_to_nc2a(:)           !< Number of points in subset Sc2, halo A, for each processor

   ! Cover tree and nearest neighbors
   logical,allocatable ::  local_mask(:,:,:)        !< Local mask
   logical,allocatable ::  displ_mask(:,:,:)        !< Displacement mask
   integer,allocatable :: nn_c2_index(:,:,:)        !< Nearest diagnostic neighbors from diagnostic points
   real(kind_real),allocatable :: nn_c2_dist(:,:,:) !< Nearest diagnostic neighbors distance from diagnostic points
   integer,allocatable :: nn_ldwv_index(:)          !< Nearest diagnostic neighbors for local diagnostics profiles

   ! Sampling mesh
   integer :: nt                                    !< Number of triangles
   integer,allocatable :: ltri(:,:)                 !< Triangles indices
   integer,allocatable :: larc(:,:)                 !< Arcs indices
   real(kind_real),allocatable :: bdist(:)          !< Distance to the closest boundary arc

   ! Interpolations
   type(linoptype),allocatable :: h(:)              !< Horizontal interpolation from Sc2 to Sc0
   type(linoptype),allocatable :: s(:)              !< Horizontal interpolation from Sc2 to Sc1
   type(linoptype),allocatable :: d(:,:)            !< Displacement interpolation

   ! MPI distribution
   integer :: nc0c                                  !< Number of points in subset Sc0, halo C
   integer :: nc0d                                  !< Number of points in subset Sc0, halo D
   integer :: nc1a                                  !< Number of points in subset Sc1, halo A
   logical,allocatable :: lcheck_c0a(:)             !< Detection of halo A on subset Sc0
   logical,allocatable :: lcheck_c0c(:)             !< Detection of halo C on subset Sc0
   logical,allocatable :: lcheck_c0d(:)             !< Detection of halo D on subset Sc0
   logical,allocatable :: lcheck_c1a(:)             !< Detection of halo A on subset Sc1
   logical,allocatable :: lcheck_d(:,:,:)           !< Detection of displacement interpolation coefficients
   integer,allocatable :: c0c_to_c0(:)              !< Subset Sc0, halo C to global
   integer,allocatable :: c0_to_c0c(:)              !< Subset Sc0, global to halo C
   integer,allocatable :: c0a_to_c0c(:)             !< Subset Sc0, halo A to halo C
   integer,allocatable :: c0d_to_c0(:)              !< Subset Sc0, halo D to global
   integer,allocatable :: c0_to_c0d(:)              !< Subset Sc0, global to halo D
   integer,allocatable :: c0a_to_c0d(:)             !< Subset Sc0, halo A to halo D
   integer,allocatable :: c1a_to_c1(:)              !< Subset Sc1, halo A to global
   integer,allocatable :: c1_to_c1a(:)              !< Subset Sc1, global to halo A
   type(comtype) :: AC                              !< Communication between halos A and C
   type(comtype) :: AD                              !< Communication between halos A and D
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
allocate(hdata%c1_to_c0(nam%nc1))
allocate(hdata%c1l0_log(nam%nc1,geom%nl0))
allocate(hdata%c1c3_to_c0(nam%nc1,nam%nc3))
allocate(hdata%c1c3l0_log(nam%nc1,nam%nc3,geom%nl0))
if (nam%local_diag.or.nam%displ_diag) then
   allocate(hdata%c2_to_c1(hdata%nc2))
   allocate(hdata%c2_to_c0(hdata%nc2))
   allocate(hdata%c2_to_c2a(hdata%nc2))
   allocate(hdata%c2_to_proc(hdata%nc2))
   allocate(hdata%proc_to_nc2a(mpl%nproc))
end if
if (nam%local_diag.or.nam%displ_diag) then
   allocate(hdata%local_mask(nam%nc1,hdata%nc2,geom%nl0i))
   allocate(hdata%displ_mask(nam%nc1,hdata%nc2,geom%nl0i))
   if (trim(nam%flt_type)/='none') then
      allocate(hdata%nn_c2_index(hdata%nc2,hdata%nc2,geom%nl0i))
      allocate(hdata%nn_c2_dist(hdata%nc2,hdata%nc2,geom%nl0i))
   end if
   allocate(hdata%h(geom%nl0i))
   allocate(hdata%s(geom%nl0i))
end if

! Initialization
call msi(hdata%c1_to_c0)
hdata%c1l0_log = .false.
call msi(hdata%c1c3_to_c0)
hdata%c1c3l0_log = .false.
if (nam%local_diag.or.nam%displ_diag) then
   call msi(hdata%c2_to_c1)
   call msi(hdata%c2_to_c0)
   call msi(hdata%c2_to_c2a)
   call msi(hdata%c2_to_proc)
   call msi(hdata%proc_to_nc2a)
end if
if (nam%local_diag.or.nam%displ_diag) then
   hdata%local_mask = .false.
   hdata%displ_mask = .false.
   if (trim(nam%flt_type)/='none') then
      call msi(hdata%nn_c2_index)
      call msr(hdata%nn_c2_dist)
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

! Release memory
if (allocated(hdata%c1_to_c0)) deallocate(hdata%c1_to_c0)
if (allocated(hdata%c1l0_log)) deallocate(hdata%c1l0_log)
if (allocated(hdata%c1c3_to_c0)) deallocate(hdata%c1c3_to_c0)
if (allocated(hdata%c1c3l0_log)) deallocate(hdata%c1c3l0_log)
if (allocated(hdata%c2_to_c1)) deallocate(hdata%c2_to_c1)
if (allocated(hdata%c2_to_c0)) deallocate(hdata%c2_to_c0)
if (allocated(hdata%c2a_to_c2)) deallocate(hdata%c2a_to_c2)
if (allocated(hdata%c2_to_c2a)) deallocate(hdata%c2_to_c2a)
if (allocated(hdata%local_mask)) deallocate(hdata%local_mask)
if (allocated(hdata%displ_mask)) deallocate(hdata%displ_mask)
if (allocated(hdata%nn_c2_index)) deallocate(hdata%nn_c2_index)
if (allocated(hdata%nn_c2_dist)) deallocate(hdata%nn_c2_dist)
if (allocated(hdata%h)) then
   do il0=1,hdata%geom%nl0
      call linop_dealloc(hdata%h(il0))
   end do
   deallocate(hdata%h)
end if
if (allocated(hdata%s)) then
   do il0=1,hdata%geom%nl0
      call linop_dealloc(hdata%s(il0))
   end do
   deallocate(hdata%s)
end if
if (allocated(hdata%nn_ldwv_index)) deallocate(hdata%nn_ldwv_index)

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
integer :: il0,il0i,ic1,jc3,ic2
integer :: nl0_test,nl0r_test,nc_test,nc1_test,nc2_test
integer :: info,ncid,nl0_id,nc_id,nc1_id,nc2_id
integer :: c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id,local_mask_id,displ_mask_id,nn_c2_index_id,nn_c2_dist_id
integer :: c1l0_logint(hdata%nam%nc1,hdata%geom%nl0),c1c3l0_logint(hdata%nam%nc1,hdata%nam%nc3,hdata%geom%nl0)
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
call ncerr(subr,nf90_inq_dimid(ncid,'nl0',nl0_id))
call ncerr(subr,nf90_inquire_dimension(ncid,nl0_id,len=nl0_test))
call ncerr(subr,nf90_get_att(ncid,nf90_global,'nl0r',nl0r_test))
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
if ((geom%nl0/=nl0_test).or.(nam%nl0r/=nl0r_test).or.(nam%nc3/=nc_test).or.(nam%nc1/=nc1_test)) then
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
call ncerr(subr,nf90_inq_varid(ncid,'c1_to_c0',c1_to_c0_id))
call ncerr(subr,nf90_inq_varid(ncid,'c1l0_log',c1l0_log_id))
call ncerr(subr,nf90_inq_varid(ncid,'c1c3_to_c0',c1c3_to_c0_id))
call ncerr(subr,nf90_inq_varid(ncid,'c1c3l0_log',c1c3l0_log_id))
if ((hdata_read==0).and.(nam%local_diag.or.nam%displ_diag)) then
   call ncerr(subr,nf90_inq_varid(ncid,'c2_to_c1',c2_to_c1_id))
   call ncerr(subr,nf90_inq_varid(ncid,'c2_to_c0',c2_to_c0_id))
end if

! Read arrays
call ncerr(subr,nf90_get_var(ncid,c1_to_c0_id,hdata%c1_to_c0))
call ncerr(subr,nf90_get_var(ncid,c1l0_log_id,c1l0_logint))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (c1l0_logint(ic1,il0)==0) then
         hdata%c1l0_log(ic1,il0) = .false.
      else if (c1l0_logint(ic1,il0)==1) then
         hdata%c1l0_log(ic1,il0) = .true.
      end if
   end do
end do
call ncerr(subr,nf90_get_var(ncid,c1c3_to_c0_id,hdata%c1c3_to_c0))
call ncerr(subr,nf90_get_var(ncid,c1c3l0_log_id,c1c3l0_logint))
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (c1c3l0_logint(ic1,jc3,il0)==0) then
            hdata%c1c3l0_log(ic1,jc3,il0) = .false.
         else if (c1c3l0_logint(ic1,jc3,il0)==1) then
            hdata%c1c3l0_log(ic1,jc3,il0) = .true.
         end if
      end do
   end do
end do
if ((hdata_read==0).and.(nam%local_diag.or.nam%displ_diag)) then
   call ncerr(subr,nf90_get_var(ncid,c2_to_c1_id,hdata%c2_to_c1))
   call ncerr(subr,nf90_get_var(ncid,c2_to_c0_id,hdata%c2_to_c0))
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
         info = nf90_inq_varid(ncid,'nn_c2_index',nn_c2_index_id)
         if (info==nf90_noerr) then
            call ncerr(subr,nf90_inq_varid(ncid,'nn_c2_dist',nn_c2_dist_id))
            call ncerr(subr,nf90_get_var(ncid,nn_c2_index_id,hdata%nn_c2_index(:,:,il0i)))
            call ncerr(subr,nf90_get_var(ncid,nn_c2_dist_id,hdata%nn_c2_dist(:,:,il0i)))
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
integer :: il0,il0i,ic1,jc3,ic2
integer :: ncid,nl0_id,nc1_id,nc2_1_id,nc2_2_id,nc_id
integer :: lat_id,lon_id,smax_id,c1_to_c0_id,c1l0_log_id,c1c3_to_c0_id,c1c3l0_log_id
integer :: c2_to_c1_id,c2_to_c0_id,local_mask_id,displ_mask_id,nn_c2_index_id,nn_c2_dist_id
integer :: c1l0_logint(hdata%nam%nc1,hdata%geom%nl0),c1c3l0_logint(hdata%nam%nc1,hdata%nam%nc3,hdata%geom%nl0)
integer,allocatable :: local_maskint(:,:),displ_maskint(:,:)
real(kind_real) :: lon(hdata%nam%nc1,hdata%nam%nc3,hdata%geom%nl0),lat(hdata%nam%nc1,hdata%nam%nc3,hdata%geom%nl0)
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
call ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
call ncerr(subr,nf90_put_att(ncid,nf90_global,'nl0r',nam%nl0r))
call ncerr(subr,nf90_def_dim(ncid,'nc',nam%nc3,nc_id))
call ncerr(subr,nf90_def_dim(ncid,'nc1',nam%nc1,nc1_id))
if (nam%local_diag.or.nam%displ_diag) call ncerr(subr,nf90_def_dim(ncid,'nc2_1',hdata%nc2,nc2_1_id))

! Define arrays
call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc1_id,nc_id,nl0_id/),lat_id))
call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc1_id,nc_id,nl0_id/),lon_id))
call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'smax',ncfloat,(/nc_id,nl0_id/),smax_id))
call ncerr(subr,nf90_put_att(ncid,smax_id,'_FillValue',msvalr))
call ncerr(subr,nf90_def_var(ncid,'c1_to_c0',nf90_int,(/nc1_id/),c1_to_c0_id))
call ncerr(subr,nf90_put_att(ncid,c1_to_c0_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'c1l0_log',nf90_int,(/nc1_id,nl0_id/),c1l0_log_id))
call ncerr(subr,nf90_put_att(ncid,c1l0_log_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'c1c3_to_c0',nf90_int,(/nc1_id,nc_id/),c1c3_to_c0_id))
call ncerr(subr,nf90_put_att(ncid,c1c3_to_c0_id,'_FillValue',msvali))
call ncerr(subr,nf90_def_var(ncid,'c1c3l0_log',nf90_int,(/nc1_id,nc_id,nl0_id/),c1c3l0_log_id))
call ncerr(subr,nf90_put_att(ncid,c1c3l0_log_id,'_FillValue',msvali))
if (nam%local_diag.or.nam%displ_diag) then
   call ncerr(subr,nf90_def_var(ncid,'c2_to_c1',nf90_int,(/nc2_1_id/),c2_to_c1_id))
   call ncerr(subr,nf90_put_att(ncid,c2_to_c1_id,'_FillValue',msvali))
   call ncerr(subr,nf90_def_var(ncid,'c2_to_c0',nf90_int,(/nc2_1_id/),c2_to_c0_id))
   call ncerr(subr,nf90_put_att(ncid,c2_to_c0_id,'_FillValue',msvali))
end if

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write arrays
call msr(lon)
call msr(lat)
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (hdata%c1c3l0_log(ic1,jc3,il0)) then
            lon(ic1,jc3,il0) = hdata%geom%lon(hdata%c1c3_to_c0(ic1,jc3))*rad2deg
            lat(ic1,jc3,il0) = hdata%geom%lat(hdata%c1c3_to_c0(ic1,jc3))*rad2deg
         end if
      end do
   end do
end do
call ncerr(subr,nf90_put_var(ncid,lon_id,lon))
call ncerr(subr,nf90_put_var(ncid,lat_id,lat))
call ncerr(subr,nf90_put_var(ncid,smax_id,float(count(hdata%c1c3l0_log,dim=1))))
call ncerr(subr,nf90_put_var(ncid,c1_to_c0_id,hdata%c1_to_c0))
do il0=1,geom%nl0
   do ic1=1,nam%nc1
      if (hdata%c1l0_log(ic1,il0)) then
         c1l0_logint(ic1,il0) = 1
      else
         c1l0_logint(ic1,il0) = 0
      end if
   end do
end do
call ncerr(subr,nf90_put_var(ncid,c1l0_log_id,c1l0_logint))
call ncerr(subr,nf90_put_var(ncid,c1c3_to_c0_id,hdata%c1c3_to_c0))
do il0=1,geom%nl0
   do jc3=1,nam%nc3
      do ic1=1,nam%nc1
         if (hdata%c1c3l0_log(ic1,jc3,il0)) then
            c1c3l0_logint(ic1,jc3,il0) = 1
         else
            c1c3l0_logint(ic1,jc3,il0) = 0
         end if
      end do
   end do
end do
call ncerr(subr,nf90_put_var(ncid,c1c3l0_log_id,c1c3l0_logint))
if (nam%local_diag.or.nam%displ_diag) then
   call ncerr(subr,nf90_put_var(ncid,c2_to_c1_id,hdata%c2_to_c1))
   call ncerr(subr,nf90_put_var(ncid,c2_to_c0_id,hdata%c2_to_c0))
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
         call ncerr(subr,nf90_def_var(ncid,'nn_c2_index',nf90_int,(/nc2_1_id,nc2_2_id/),nn_c2_index_id))
         call ncerr(subr,nf90_put_att(ncid,nn_c2_index_id,'_FillValue',msvali))
         call ncerr(subr,nf90_def_var(ncid,'nn_c2_dist',ncfloat,(/nc2_1_id,nc2_2_id/),nn_c2_dist_id))
         call ncerr(subr,nf90_put_att(ncid,nn_c2_dist_id,'_FillValue',msvalr))
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
         call ncerr(subr,nf90_put_var(ncid,nn_c2_index_id,hdata%nn_c2_index(:,:,il0i)))
         call ncerr(subr,nf90_put_var(ncid,nn_c2_dist_id,hdata%nn_c2_dist(:,:,il0i)))
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
