!----------------------------------------------------------------------
! Module: type_io
!> Purpose: I/O derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module type_io

use netcdf
use tools_const, only: pi,deg2rad,rad2deg,reqkm,msvalr
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isanymsi
use tools_nc, only: ncfloat
use tools_qsort, only: qsort
use type_com, only: com_type
use type_geom, only: geom_type
use type_linop, only: linop_type
use type_mpl, only: mpl_type
use type_nam, only: nam_type
use type_rng, only: rng_type

implicit none

! I/O derived type
type io_type
   integer :: nlon                        !< Number of output grid longitudes
   integer :: nlat                        !< Number of output grid latitudes
   real(kind_real),allocatable :: lon(:)  !< Output grid longitudes
   real(kind_real),allocatable :: lat(:)  !< Output grid latitudes
   integer,allocatable :: og_to_lon(:)    !< Output grid to longitude index
   integer,allocatable :: og_to_lat(:)    !< Output grid to latitude index
   integer :: nog                         !< Output grid size
   integer :: noga                        !< Output grid size, halo A
   integer :: nc0b                        !< Subset Sc0 size, halo B
   integer,allocatable :: og_to_proc(:)   !< Output grid to processor
   integer,allocatable :: proc_to_noga(:) !< Number of points in output grid, halo A, for each processor
   integer,allocatable :: oga_to_og(:)    !< Output grid, halo A to global
   integer,allocatable :: og_to_oga(:)    !< Output grid, global to halo A
   integer,allocatable :: c0b_to_c0(:)    !< Subset Sc0, halo B to global
   integer,allocatable :: c0_to_c0b(:)    !< Subset Sc0, global to halo B
   integer,allocatable :: c0a_to_c0b(:)   !< Subset Sc0, halo A to halo B
   type(linop_type) :: og                 !< Subset Sc0 to grid interpolation
   type(com_type) :: com_AB               !< Communication between halos A and B
   logical,allocatable :: mask(:,:)       !< Mask on output grid
contains
   procedure :: dealloc => io_dealloc
   procedure :: fld_read => io_fld_read
   procedure :: fld_write => io_fld_write
   procedure :: grid_init => io_grid_init
   procedure :: grid_write => io_grid_write
end type io_type

private
public :: io_type

contains

!----------------------------------------------------------------------
! Subroutine: io_dealloc
!> Purpose: deallocate I/O
!----------------------------------------------------------------------
subroutine io_dealloc(io)

implicit none

! Passed variables
class(io_type),intent(inout) :: io !< I/O

! Release memory
if (allocated(io%lon)) deallocate(io%lon)
if (allocated(io%lat)) deallocate(io%lat)
if (allocated(io%og_to_lon)) deallocate(io%og_to_lon)
if (allocated(io%og_to_lat)) deallocate(io%og_to_lat)
if (allocated(io%og_to_proc)) deallocate(io%og_to_proc)
if (allocated(io%proc_to_noga)) deallocate(io%proc_to_noga)
if (allocated(io%oga_to_og)) deallocate(io%oga_to_og)
if (allocated(io%og_to_oga)) deallocate(io%og_to_oga)
if (allocated(io%c0b_to_c0)) deallocate(io%c0b_to_c0)
if (allocated(io%c0_to_c0b)) deallocate(io%c0_to_c0b)
if (allocated(io%c0a_to_c0b)) deallocate(io%c0a_to_c0b)
call io%og%dealloc
call io%com_AB%dealloc
if (allocated(io%mask)) deallocate(io%mask)

end subroutine io_dealloc

!----------------------------------------------------------------------
! Subroutine: io_fld_read
!> Purpose: write field
!----------------------------------------------------------------------
subroutine io_fld_read(io,mpl,nam,geom,filename,varname,fld)

implicit none

! Passed variables
class(io_type),intent(in) :: io                        !< I/O
type(mpl_type),intent(in) :: mpl                       !< MPI data
type(nam_type),intent(in) :: nam                       !< Namelist
type(geom_type),intent(in) :: geom                     !< Geometry
character(len=*),intent(in) :: filename                !< File name
character(len=*),intent(in) :: varname                 !< Variable name
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
integer :: ncid,fld_id,il0,dum
real(kind_real) :: fld_c0(geom%nc0,geom%nl0)
character(len=1024) :: filename_proc
character(len=1024) :: subr = 'fld_read'

if (nam%field_io) then
   if (nam%split_io.and.(mpl%nproc>1)) then
      ! Open file
      write(filename_proc,'(a,a,a,a,i4.4,a,i4.4,a)') trim(nam%datadir),'/',trim(filename),'_',mpl%nproc,'-',mpl%myproc,'.nc'
      call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename),nf90_nowrite,ncid))

      ! Get variable id
      call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

      ! Get data
      call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld))

      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   else
      if (mpl%main) then

         ! Open file
         call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_nowrite,ncid))

         ! Get variable id
         call mpl%ncerr(subr,nf90_inq_varid(ncid,trim(varname),fld_id))

         ! Get data
         call mpl%ncerr(subr,nf90_get_var(ncid,fld_id,fld_c0))

         ! Close file
         call mpl%ncerr(subr,nf90_close(ncid))
      end if

      ! Global to local
      do il0=1,geom%nl0
         call mpl%scatterv(geom%proc_to_nc0a,geom%nc0,fld_c0(:,il0),geom%nc0a,fld(:,il0))
      end do
   end if
else
   ! No field I/O
   call mpl%abort('field/variable not read: '//trim(filename)//'/'//trim(varname))
end if

! Dummy call to avoid compilation warnings
dum = io%nog

end subroutine io_fld_read

!----------------------------------------------------------------------
! Subroutine: io_fld_write
!> Purpose: write field
!----------------------------------------------------------------------
subroutine io_fld_write(io,mpl,nam,geom,filename,varname,fld)

implicit none

! Passed variables
class(io_type),intent(in) :: io                       !< I/O
type(mpl_type),intent(inout) :: mpl                   !< MPI data
type(nam_type),intent(in) :: nam                      !< Namelist
type(geom_type),intent(in) :: geom                    !< Geometry
character(len=*),intent(in) :: filename               !< File name
character(len=*),intent(in) :: varname                !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
integer :: ic0a,ic0,il0,info
integer :: ncid,nc0a_id,nc0_id,nl0_id,fld_id,lon_id,lat_id
real(kind_real) :: fld_c0a(geom%nc0a,geom%nl0)
real(kind_real) :: fld_c0(geom%nc0,geom%nl0)
character(len=1024) :: filename_proc
character(len=1024) :: subr = 'io_fld_write'

if (nam%field_io) then
   ! Apply mask
   do il0=1,geom%nl0
      do ic0a=1,geom%nc0a
         ic0 = geom%c0a_to_c0(ic0a)
         if (geom%mask(ic0,il0)) then
            fld_c0a(ic0a,il0) = fld(ic0a,il0)
         else
            call msr(fld_c0a(ic0a,il0))
         end if
      end do
   end do

   if (nam%split_io.and.(mpl%nproc>1)) then
      ! Check if the file exists
      write(filename_proc,'(a,a,a,a,i4.4,a,i4.4,a)') trim(nam%datadir),'/',trim(filename),'_',mpl%nproc,'-',mpl%myproc,'.nc'
      info = nf90_create(trim(filename_proc),or(nf90_noclobber,nf90_64bit_offset),ncid)
      if (info==nf90_noerr) then
         ! Write namelist parameters
          call nam%ncwrite(mpl,ncid)

         ! Define attribute
         call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))

         ! End definition mode
         call mpl%ncerr(subr,nf90_enddef(ncid))
      else
         ! Open file
         call mpl%ncerr(subr,nf90_open(trim(filename_proc),nf90_write,ncid))
      end if

      ! Get variable id
      info = nf90_inq_varid(ncid,trim(varname),fld_id)

      ! Define dimensions and variable if necessary
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_redef(ncid))
         info = nf90_inq_dimid(ncid,'nc0a',nc0a_id)
         if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0a',geom%nc0a,nc0a_id))
         info = nf90_inq_dimid(ncid,'nl0',nl0_id)
         if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
         call mpl%ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nc0a_id,nl0_id/),fld_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
         call mpl%ncerr(subr,nf90_enddef(ncid))
      end if

      ! Write data
      call mpl%ncerr(subr,nf90_put_var(ncid,fld_id,fld_c0a))

      ! Write coordinates
      info = nf90_inq_varid(ncid,'lon',lon_id)
      if (info/=nf90_noerr) then
         call mpl%ncerr(subr,nf90_redef(ncid))
         call mpl%ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0a_id/),lon_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
         call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
         call mpl%ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0a_id/),lat_id))
         call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
         call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
         call mpl%ncerr(subr,nf90_enddef(ncid))
         call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon(geom%c0a_to_c0)*rad2deg))
         call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat(geom%c0a_to_c0)*rad2deg))
      end if

      ! Close file
      call mpl%ncerr(subr,nf90_close(ncid))
   else
      ! Local to global
      do il0=1,geom%nl0
         call mpl%gatherv(geom%nc0a,fld_c0a(:,il0),geom%proc_to_nc0a,geom%nc0,fld_c0(:,il0))
      end do

      if (mpl%main) then
         ! Check if the file exists
         info = nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_noclobber,nf90_64bit_offset),ncid)
         if (info==nf90_noerr) then
            ! Write namelist parameters
            call nam%ncwrite(mpl,ncid)

            ! Define attribute
            call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))

            ! End definition mode
            call mpl%ncerr(subr,nf90_enddef(ncid))
         else
            ! Open file
            call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_write,ncid))
         end if

         ! Get variable id
         info = nf90_inq_varid(ncid,trim(varname),fld_id)

         ! Define dimensions and variable if necessary
         if (info/=nf90_noerr) then
            call mpl%ncerr(subr,nf90_redef(ncid))
            info = nf90_inq_dimid(ncid,'nc0',nc0_id)
            if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
            info = nf90_inq_dimid(ncid,'nl0',nl0_id)
            if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'nl0',geom%nl0,nl0_id))
            call mpl%ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nc0_id,nl0_id/),fld_id))
            call mpl%ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
            call mpl%ncerr(subr,nf90_enddef(ncid))
         end if

         ! Write data
          call mpl%ncerr(subr,nf90_put_var(ncid,fld_id,fld_c0))

         ! Write coordinates
         info = nf90_inq_varid(ncid,'lon',lon_id)
         if (info/=nf90_noerr) then
            call mpl%ncerr(subr,nf90_redef(ncid))
            call mpl%ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
            call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
            call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
            call mpl%ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
            call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
            call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
            call mpl%ncerr(subr,nf90_enddef(ncid))
            call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
            call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
         end if

         ! Close file
         call mpl%ncerr(subr,nf90_close(ncid))
      end if
   end if

   ! Regridded field output
   if (nam%grid_output) call io%grid_write(mpl,nam,geom,filename,trim(varname)//'_regridded',fld)
else
   ! No field I/O
   call mpl%warning('field/variable not written: '//trim(filename)//'/'//trim(varname))
end if

end subroutine io_fld_write

!----------------------------------------------------------------------
! Subroutine: io_grid_init
!> Purpose: initialize fields regridding
!----------------------------------------------------------------------
subroutine io_grid_init(io,mpl,rng,nam,geom)

implicit none

! Passed variables
class(io_type),intent(inout) :: io  !< I/O
type(mpl_type),intent(inout) :: mpl !< MPI data
type(rng_type),intent(inout) :: rng !< Random number generator
type(nam_type),intent(in) :: nam    !< Namelist
type(geom_type),intent(in) :: geom  !< Geometry

! Local variables
integer ::ilon,ilat,i_s,i_s_loc,iog,iproc,ic0,ic0a,ic0b,ioga,il0
integer,allocatable :: order(:),order_inv(:),interpg_lg(:)
real(kind_real) :: dlon,dlat
real(kind_real),allocatable :: lon_og(:),lat_og(:)
logical :: mask_c0(geom%nc0)
logical,allocatable :: mask_lonlat(:,:),mask_og(:),lcheck_og(:),lcheck_c0b(:)
type(linop_type) :: ogfull

! Grid size
io%nlat = nint(pi/nam%grid_resol)
io%nlon = 2*io%nlat
dlon = 2.0*pi/real(io%nlon,kind_real)
dlat = pi/real(io%nlat,kind_real)

! Allocation
allocate(io%lon(io%nlon))
allocate(io%lat(io%nlat))
allocate(mask_lonlat(io%nlon,io%nlat))

! Setup output grid
do ilat=1,io%nlat
   do ilon=1,io%nlon
      io%lon(ilon) = (-pi+dlon/2)+real(ilon-1,kind_real)*dlon
      io%lat(ilat) = (-pi/2+dlat/2)+real(ilat-1,kind_real)*dlat
      call geom%mesh%inside(io%lon(ilon),io%lat(ilat),mask_lonlat(ilon,ilat))
   end do
end do
io%nog = count(mask_lonlat)

! Print results
write(mpl%unit,'(a7,a)') '','Output grid:'
write(mpl%unit,'(a10,a,f7.2,a,f5.2,a)') '','Effective resolution: ',0.5*(dlon+dlat)*reqkm,' km (', &
 & 0.5*(dlon+dlat)*rad2deg,' deg.)'
write(mpl%unit,'(a10,a,i4,a,i4)') '',      'Size (nlon x nlat):   ',io%nlon,' x ',io%nlat

! Allocation
allocate(io%og_to_lon(io%nog))
allocate(io%og_to_lat(io%nog))
allocate(lon_og(io%nog))
allocate(lat_og(io%nog))
allocate(mask_og(io%nog))
allocate(io%og_to_proc(io%nog))
allocate(order(io%nog))
allocate(order_inv(io%nog))
allocate(io%proc_to_noga(mpl%nproc))
allocate(io%og_to_oga(io%nog))

! Setup packed output grid
iog = 0
do ilon=1,io%nlon
   do ilat=1,io%nlat
      if (mask_lonlat(ilon,ilat)) then
         iog = iog+1
         lon_og(iog) = io%lon(ilon)
         lat_og(iog) = io%lat(ilat)
         io%og_to_lon(iog) = ilon
         io%og_to_lat(iog) = ilat
      end if
   end do
end do
mask_og = .true.

! Interpolation setup
mask_c0 = .true.
call ogfull%interp(mpl,rng,geom%nc0,geom%lon,geom%lat,mask_c0,io%nog,lon_og,lat_og,mask_og,nam%grid_interp)

! Define a processor for each output grid point
call msi(io%og_to_proc)
do i_s=1,ogfull%n_s
   ic0 = ogfull%col(i_s)
   iproc = geom%c0_to_proc(ic0)
   iog = ogfull%row(i_s)
   io%og_to_proc(iog) = iproc
end do
if (isanymsi(io%og_to_proc)) call mpl%abort('some output grid points are not interpolated')
do iproc=1,mpl%nproc
   io%proc_to_noga(iproc) = count(io%og_to_proc==iproc)
end do
io%noga = io%proc_to_noga(mpl%myproc)

! Reorder points
call qsort(io%nog,io%og_to_proc,order)
do iog=1,io%nog
   order_inv(order(iog)) = iog
end do
io%og_to_lon = io%og_to_lon(order)
io%og_to_lat = io%og_to_lat(order)
do i_s=1,ogfull%n_s
   ogfull%row(i_s) = order_inv(ogfull%row(i_s))
end do

! Conversions
allocate(io%oga_to_og(io%noga))
iog=0
do iproc=1,mpl%nproc
   do ioga=1,io%proc_to_noga(iproc)
      iog = iog+1
      io%og_to_oga(iog) = ioga
      if (iproc==mpl%myproc) io%oga_to_og(ioga) = iog
   end do
end do

! Allocation
allocate(lcheck_og(ogfull%n_s))
allocate(lcheck_c0b(geom%nc0))

! Halo definitions

! Halo B
lcheck_og = .false.
lcheck_c0b = .false.
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   lcheck_c0b(ic0) = .true.
end do
do i_s=1,ogfull%n_s
   iog = ogfull%row(i_s)
   iproc = io%og_to_proc(iog)
   if (iproc==mpl%myproc) then
      ic0 = ogfull%col(i_s)
      lcheck_og(i_s) = .true.
      lcheck_c0b(ic0) = .true.
   end if
end do

! Halo sizes
io%nc0b = count(lcheck_c0b)
io%og%n_s = count(lcheck_og)

! Halo B
allocate(io%c0b_to_c0(io%nc0b))
allocate(io%c0_to_c0b(geom%nc0))
ic0b = 0
do ic0=1,geom%nc0
   if (lcheck_c0b(ic0)) then
      ic0b = ic0b+1
      io%c0b_to_c0(ic0b) = ic0
      io%c0_to_c0b(ic0) = ic0b
   end if
end do

! Inter-halo conversions
allocate(io%c0a_to_c0b(geom%nc0a))
do ic0a=1,geom%nc0a
   ic0 = geom%c0a_to_c0(ic0a)
   ic0b = io%c0_to_c0b(ic0)
   io%c0a_to_c0b(ic0a) = ic0b
end do

! Global <-> local conversions for data
allocate(interpg_lg(io%og%n_s))
i_s_loc = 0
do i_s=1,ogfull%n_s
   if (lcheck_og(i_s)) then
      i_s_loc = i_s_loc+1
      interpg_lg(i_s_loc) = i_s
   end if
end do

! Local data

! Horizontal interpolation
write(io%og%prefix,'(a,i3.3)') 'og'
io%og%n_src = io%nc0b
io%og%n_dst = io%noga
call io%og%alloc
do i_s_loc=1,io%og%n_s
   i_s = interpg_lg(i_s_loc)
   io%og%row(i_s_loc) = io%og_to_oga(ogfull%row(i_s))
   io%og%col(i_s_loc) = io%c0_to_c0b(ogfull%col(i_s))
   io%og%S(i_s_loc) = ogfull%S(i_s)
end do
call io%og%reorder(mpl)

! Setup communications
call io%com_AB%setup(mpl,'com_AB',geom%nc0,geom%nc0a,io%nc0b,io%c0b_to_c0,io%c0a_to_c0b,geom%c0_to_proc,geom%c0_to_c0a)

! Compute output grid mask
allocate(io%mask(io%noga,geom%nl0))
io%mask = .true.
do i_s=1,io%og%n_s
   ic0b = io%og%col(i_s)
   ioga = io%og%row(i_s)
   ic0 = io%c0b_to_c0(ic0b)
   do il0=1,geom%nl0
      if (.not.geom%mask(ic0,il0)) io%mask(ioga,il0) = .false.
   end do
end do

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0 =    ',geom%nc0
write(mpl%unit,'(a10,a,i8)') '','nc0a =   ',geom%nc0a
write(mpl%unit,'(a10,a,i8)') '','nc0b =   ',io%nc0b
write(mpl%unit,'(a10,a,i8)') '','nog =    ',io%nog
write(mpl%unit,'(a10,a,i8)') '','noga =   ',io%noga
write(mpl%unit,'(a10,a,i8)') '','og%n_s = ',io%og%n_s
call flush(mpl%unit)

end subroutine io_grid_init

!----------------------------------------------------------------------
! Subroutine: io_grid_write
!> Purpose: interpolate and write field
!----------------------------------------------------------------------
subroutine io_grid_write(io,mpl,nam,geom,filename,varname,fld)

implicit none

! Passed variables
class(io_type),intent(in) :: io                       !< I/O
type(mpl_type),intent(inout) :: mpl                   !< MPI data
type(nam_type),intent(in) :: nam                      !< Namelist
type(geom_type),intent(in) :: geom                    !< Geometry
character(len=*),intent(in) :: filename               !< File name
character(len=*),intent(in) :: varname                !< Variable name
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
integer :: il0,info,ilon,ilat,iog,ioga,i,iproc
integer :: ncid,nlon_id,nlat_id,nlev_id,fld_id,lon_id,lat_id,lev_id
integer,allocatable :: oga_to_og(:)
real(kind_real) :: fld_c0b(io%nc0b,geom%nl0)
real(kind_real) :: fld_oga(io%noga,geom%nl0)
real(kind_real),allocatable :: sbuf(:),rbuf(:),fld_grid(:,:,:)
character(len=1024) :: subr = 'io_grid_write'

! Halo extension and interpolation
call io%com_AB%ext(mpl,geom%nl0,fld,fld_c0b)
do il0=1,geom%nl0
   call io%og%apply(mpl,fld_c0b(:,il0),fld_oga(:,il0),mssrc=.true.)
end do

! Apply mask
do il0=1,geom%nl0
   do ioga=1,io%noga
      if (.not.io%mask(ioga,il0)) call msr(fld_oga(ioga,il0))
   end do
end do

! Allocation
allocate(sbuf(io%noga*geom%nl0))

! Prepare buffer
i = 1
do il0=1,geom%nl0
   do ioga=1,io%noga
      sbuf(i) = fld_oga(ioga,il0)
      i = i+1
   end do
end do

if (mpl%main) then
   ! Allocation
   allocate(fld_grid(io%nlon,io%nlat,geom%nl0))

   ! Initialization
   call msr(fld_grid)

   do iproc=1,mpl%nproc
      ! Allocation
      allocate(oga_to_og(io%proc_to_noga(iproc)))
      allocate(rbuf(io%proc_to_noga(iproc)*geom%nl0))

      if (iproc==mpl%ioproc) then
         ! Copy buffer
         oga_to_og = io%oga_to_og
         rbuf = sbuf
      else
         ! Receive data on ioproc
         call mpl%recv(io%proc_to_noga(iproc),oga_to_og,iproc,mpl%tag)
         call mpl%recv(io%proc_to_noga(iproc)*geom%nl0,rbuf,iproc,mpl%tag+1)
      end if

      ! Write data
      i = 1
      do il0=1,geom%nl0
         do ioga=1,io%proc_to_noga(iproc)
            iog = oga_to_og(ioga)
            ilon = io%og_to_lon(iog)
            ilat = io%og_to_lat(iog)
            fld_grid(ilon,ilat,il0) = rbuf(i)
            i = i+1
         end do
      end do

      ! Release memory
      deallocate(oga_to_og)
      deallocate(rbuf)
   end do
else
   ! Send data to ioproc
   call mpl%send(io%noga,io%oga_to_og,mpl%ioproc,mpl%tag)
   call mpl%send(io%noga*geom%nl0,sbuf,mpl%ioproc,mpl%tag+1)
end if
mpl%tag = mpl%tag+2

! Release memory
deallocate(sbuf)

if (mpl%main) then
   ! Check if the file exists
   info = nf90_create(trim(nam%datadir)//'/'//trim(filename)//'.nc',or(nf90_noclobber,nf90_64bit_offset),ncid)
   if (info==nf90_noerr) then
      ! Write namelist parameters
      call nam%ncwrite(mpl,ncid)

      ! Define attribute
      call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,'_FillValue',msvalr))

      ! End definition mode
      call mpl%ncerr(subr,nf90_enddef(ncid))
   else
      ! Open file
      call mpl%ncerr(subr,nf90_open(trim(nam%datadir)//'/'//trim(filename)//'.nc',nf90_write,ncid))
   end if

   ! Get variable id
   info = nf90_inq_varid(ncid,trim(varname),fld_id)

   ! Define dimensions and variable if necessary
   if (info/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_redef(ncid))
      info = nf90_inq_dimid(ncid,'lon_gridded',nlon_id)
      if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'lon_gridded',io%nlon,nlon_id))
      info = nf90_inq_dimid(ncid,'lat_gridded',nlat_id)
      if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'lat_gridded',io%nlat,nlat_id))
      info = nf90_inq_dimid(ncid,'lev',nlev_id)
      if (info/=nf90_noerr) call mpl%ncerr(subr,nf90_def_dim(ncid,'lev',geom%nl0,nlev_id))
      call mpl%ncerr(subr,nf90_def_var(ncid,trim(varname),ncfloat,(/nlon_id,nlat_id,nlev_id/),fld_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,fld_id,'_FillValue',msvalr))
      call mpl%ncerr(subr,nf90_enddef(ncid))
   end if

   ! Write data
   call mpl%ncerr(subr,nf90_put_var(ncid,fld_id,fld_grid))

   ! Write coordinates
   info = nf90_inq_varid(ncid,'lon_gridded',lon_id)
   if (info/=nf90_noerr) then
      call mpl%ncerr(subr,nf90_redef(ncid))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lon_gridded',ncfloat,(/nlon_id/),lon_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lon_id,'unit','degrees_north'))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lat_gridded',ncfloat,(/nlat_id/),lat_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lat_id,'unit','degrees_east'))
      call mpl%ncerr(subr,nf90_def_var(ncid,'lev',ncfloat,(/nlev_id/),lev_id))
      call mpl%ncerr(subr,nf90_put_att(ncid,lev_id,'_FillValue',msvalr))
      call mpl%ncerr(subr,nf90_put_att(ncid,lev_id,'unit','layer'))
      call mpl%ncerr(subr,nf90_enddef(ncid))
      call mpl%ncerr(subr,nf90_put_var(ncid,lon_id,io%lon*rad2deg))
      call mpl%ncerr(subr,nf90_put_var(ncid,lat_id,io%lat*rad2deg))
      do il0=1,geom%nl0
         call mpl%ncerr(subr,nf90_put_var(ncid,lev_id,real(il0,kind_real),(/il0/)))
      end do
   end if

   ! Close file
   call mpl%ncerr(subr,nf90_close(ncid))

   ! Release memory
   deallocate(fld_grid)
end if

end subroutine io_grid_write

end module type_io
