! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2017-2019 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_obsdb_mod

use atlas_module
use datetime_mod
use duration_mod
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use netcdf
use qg_constants_mod
use qg_locs_mod
use qg_obsvec_mod
use qg_projection_mod
use qg_tools_mod
use random_mod
use string_f_c_mod

implicit none

private
public :: qg_obsdb
public :: qg_obsdb_registry
public :: qg_obsdb_setup,qg_obsdb_delete,qg_obsdb_save,qg_obsdb_get,qg_obsdb_put,qg_obsdb_locations,qg_obsdb_generate,qg_obsdb_nobs
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 1 !< Random seed (for reproducibility)

type column_data
  character(len=50) :: colname                !< Column name
  type(column_data),pointer :: next => null() !< Next column
  integer :: nlev                             !< Number of levels
  real(kind_real),allocatable :: values(:,:)  !< Values
end type column_data

type group_data
  character(len=50) :: grpname                   !< Group name
  type(group_data),pointer :: next => null()     !< Next group
  integer :: nobs                                !< Number of observations
  type(datetime),allocatable :: times(:)         !< Time-slots
  type(column_data),pointer :: colhead => null() !< Head column
end type group_data

type qg_obsdb
  integer :: ngrp                               !< Number of groups
  character(len=1024) :: filein                 !< Input filename
  character(len=1024) :: fileout                !< Output filename
  type(group_data),pointer :: grphead => null() !< Head group
end type qg_obsdb

#define LISTED_TYPE qg_obsdb

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_obsdb_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup observation data
subroutine qg_obsdb_setup(self,f_conf,winbgn,winend)
use string_utils

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self           !< Observation data
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(datetime),intent(in) :: winbgn            !< Start of window
type(datetime),intent(in) :: winend            !< End of window

! Local variables
character(len=1024) :: fin,fout
character(len=:),allocatable :: str

! Input file
if (f_conf%has("obsdatain")) then
  call f_conf%get_or_die("obsdatain.obsfile",str)
  fin = str
else
  fin = ''
endif
call fckit_log%info('qg_obsdb_setup: file in = '//trim(fin))

! Output file
if (f_conf%has("obsdataout")) then
  call f_conf%get_or_die("obsdataout.obsfile",str)
  call swap_name_member(f_conf, str, 6)

  fout = str
  call fckit_log%info('qg_obsdb_setup: file out = '//trim(fout))
else
  fout = ''
endif

! Set attributes
self%ngrp = 0
self%filein = fin
self%fileout = fout

! Read observation data
if (self%filein/='') call qg_obsdb_read(self,winbgn,winend)

end subroutine qg_obsdb_setup
! ------------------------------------------------------------------------------
!> Delete observation data
subroutine qg_obsdb_delete(self)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs

! Release memory
do while (associated(self%grphead))
  jgrp => self%grphead
  self%grphead => jgrp%next
  do jobs=1,jgrp%nobs
    call datetime_delete(jgrp%times(jobs))
  enddo
  deallocate(jgrp%times)
  do while (associated(jgrp%colhead))
    jcol => jgrp%colhead
    jgrp%colhead => jcol%next
    deallocate(jcol%values)
    deallocate(jcol)
  enddo
  deallocate(jgrp)
enddo

end subroutine qg_obsdb_delete
! ------------------------------------------------------------------------------
!> Delete observation data
subroutine qg_obsdb_save(self)
implicit none
type(qg_obsdb),intent(in) :: self !< Observation data

if (self%fileout/='') call qg_obsdb_write(self)

end subroutine qg_obsdb_save
! ------------------------------------------------------------------------------
!> Get observation data
subroutine qg_obsdb_get(self,grp,col,ovec)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self     !< Observation data
character(len=*),intent(in) :: grp    !< Group
character(len=*),intent(in) :: col    !< Column
type(qg_obsvec),intent(inout) :: ovec !< Observation vector

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs,jlev

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) then
  jgrp => self%grphead
  do while (associated(jgrp))
    call fckit_log%info('qg_obsdb_get: group '//trim(jgrp%grpname)//' exists')
    jgrp => jgrp%next
  enddo
  call fckit_log%error('qg_obsdb_get: cannot find '//trim(grp))
  call abor1_ftn('qg_obsdb_get: obs group not found')
endif

! Find observation column
call qg_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) then
  call fckit_log%error('qg_obsdb_get: cannot find '//trim(col))
  call abor1_ftn('qg_obsdb_get: obs column not found')
endif

! Get observation data
if (allocated(ovec%values)) deallocate(ovec%values)
ovec%nlev = jcol%nlev
! get all the obs
ovec%nobs = jgrp%nobs
allocate(ovec%values(ovec%nlev,ovec%nobs))
do jobs=1,jgrp%nobs
  do jlev=1,jcol%nlev
    ovec%values(jlev,jobs) = jcol%values(jlev,jobs)
  enddo
enddo

end subroutine qg_obsdb_get
! ------------------------------------------------------------------------------
!> Put observations data
subroutine qg_obsdb_put(self,grp,col,ovec)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data
character(len=*),intent(in) :: grp   !< Group
character(len=*),intent(in) :: col   !< Column
type(qg_obsvec),intent(in) :: ovec   !< Observation vector

! Local variables
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
integer :: jobs,jlev

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) then
  jgrp => self%grphead
  do while (associated(jgrp))
    call fckit_log%info('qg_obsdb_put: group '//trim(jgrp%grpname)//' exists')
    jgrp => jgrp%next
  enddo
  call fckit_log%error('qg_obsdb_put: cannot find '//trim(grp))
  call abor1_ftn('qg_obsdb_put: obs group not found')
endif

! Find observation column (and add it if not there)
call qg_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) then
  if (.not.associated(jgrp%colhead)) call abor1_ftn('qg_obsdb_put: no locations')
  jcol => jgrp%colhead
  do while (associated(jcol%next))
    jcol => jcol%next
  enddo
  allocate(jcol%next)
  jcol => jcol%next
  jcol%colname = col
  jcol%nlev = ovec%nlev
  allocate(jcol%values(jcol%nlev,jgrp%nobs))
endif

! Put observation data
if (ovec%nobs/=jgrp%nobs) call abor1_ftn('qg_obsdb_put: error obs number')
if (ovec%nlev/=jcol%nlev) call abor1_ftn('qg_obsdb_put: error col number')
do jobs=1,jgrp%nobs
  do jlev=1,jcol%nlev
    jcol%values(jlev,jobs) = ovec%values(jlev,jobs)
  enddo
enddo

end subroutine qg_obsdb_put
! ------------------------------------------------------------------------------
!> Get locations from observation data
subroutine qg_obsdb_locations(self,grp,fields,c_times)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self   !< Observation data
character(len=*),intent(in) :: grp  !< Group
type(atlas_fieldset), intent(inout) :: fields !< Locations FieldSet
type(c_ptr), intent(in), value :: c_times !< pointer to times array in C++

! Local variables
integer :: nlocs, jo
character(len=8),parameter :: col = 'Location'
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
type(atlas_field) :: field_z, field_lonlat
real(kind_real), pointer :: z(:), lonlat(:,:)

! Find observation group
call qg_obsdb_find_group(self,grp,jgrp)
if (.not.associated(jgrp)) call abor1_ftn('qg_obsdb_locations: obs group not found')
nlocs = jgrp%nobs

! Find observation column
call qg_obsdb_find_column(jgrp,col,jcol)
if (.not.associated(jcol)) call abor1_ftn('qg_obsdb_locations: obs column not found')

! Set number of observations

field_lonlat = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=[2,nlocs])
field_z = atlas_field(name="altitude", kind=atlas_real(kind_real), shape=[nlocs])

call field_lonlat%data(lonlat)
call field_z%data(z)

! Copy coordinates
do jo = 1, nlocs
  lonlat(1,jo) = jcol%values(1,jo)
  lonlat(2,jo) = jcol%values(2,jo)
  z(jo) = jcol%values(3,jo)
  call f_c_push_to_datetime_vector(c_times, jgrp%times(jo))
enddo

call fields%add(field_lonlat)
call fields%add(field_z)

! release pointers
call field_lonlat%final()
call field_z%final()

end subroutine qg_obsdb_locations
! ------------------------------------------------------------------------------
!> Generate observation data
subroutine qg_obsdb_generate(self,grp,f_conf,bgn,step,ktimes,kobs)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self           !< Observation data
character(len=*),intent(in) :: grp             !< Group
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(datetime),intent(in) :: bgn               !< Start time
type(duration),intent(in) :: step              !< Time-step
integer,intent(in) :: ktimes                   !< Number of time-slots
integer,intent(inout) :: kobs                  !< Number of observations

! Local variables
integer :: nlev,nlocs, jobs
real(kind_real),allocatable :: x(:), y(:), lon(:), lat(:), z(:)
real(kind_real) :: err
type(datetime),allocatable :: times(:)
type(qg_obsvec) :: obsloc,obserr

! Get number of observations
if (f_conf%has("obs density")) then
  call f_conf%get_or_die("obs density",nlocs)
  allocate(x(nlocs), y(nlocs), z(nlocs), lon(nlocs), lat(nlocs))
  ! Generate random locations
  call uniform_distribution(x,0.0_kind_real,domain_zonal,rseed)
  call uniform_distribution(y,0.0_kind_real,domain_meridional,rseed)
  call uniform_distribution(z,0.0_kind_real,domain_depth,rseed)
  ! Convert to lon/lat
  do jobs=1,nlocs
    call xy_to_lonlat(x(jobs),y(jobs),lon(jobs),lat(jobs))
  enddo
  deallocate(x, y)
else
  nlocs = f_conf%get_size("obs locations.lon")
  allocate(lon(nlocs), lat(nlocs), z(nlocs))
  call f_conf%get_or_die("obs locations.lon", lon)
  call f_conf%get_or_die("obs locations.lat", lat)
  call f_conf%get_or_die("obs locations.z", z)
endif

! Allocation
kobs = nlocs*ktimes
allocate(times(kobs))
! Generate locations
call qg_obsdb_generate_locations(nlocs,lon,lat,z,ktimes,bgn,step,times,obsloc)

! Create observations data
call qg_obsdb_create(self,trim(grp),times,obsloc)

! Create observation error
call f_conf%get_or_die("obs error",err)
call f_conf%get_or_die("nval",nlev)
call qg_obsvec_setup(obserr,nlev,kobs)
obserr%values(:,:) = err
call qg_obsdb_put(self,trim(grp),'ObsError',obserr)

! Release memory
deallocate(times)
deallocate(obsloc%values)
deallocate(obserr%values)
deallocate(lon,lat,z)

end subroutine qg_obsdb_generate
! ------------------------------------------------------------------------------
!> Get observation data size
subroutine qg_obsdb_nobs(self,grp,kobs)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self  !< Observation data
character(len=*),intent(in) :: grp !< Group
integer,intent(inout) :: kobs      !< Number of observations

! Local variables
type(group_data),pointer :: jgrp

! Find group
call qg_obsdb_find_group(self,grp,jgrp)

! Get observation data size
if (associated(jgrp)) then
  kobs = jgrp%nobs
else
  kobs = 0
endif

end subroutine qg_obsdb_nobs
! ------------------------------------------------------------------------------
!  Private
! ------------------------------------------------------------------------------
!> Read observation data
subroutine qg_obsdb_read(self,winbgn,winend)

implicit none

! Passed variables
type(qg_obsdb),intent(inout) :: self !< Observation data
type(datetime),intent(in) :: winbgn  !< Start of window
type(datetime),intent(in) :: winend  !< End of window

! Local variables
integer,parameter :: ngrpmax = 3
integer,parameter :: ncolmax = 7
integer :: grp_ids(ngrpmax),igrp,nobs_in_grp,iobs,jobs,ncol,col_ids(ncolmax),icol,nlev_id,values_id
integer :: ncid,nobs_id,times_id
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
character(len=50) :: stime
logical, allocatable :: inwindow(:)
type(datetime) :: tobs
type(datetime), allocatable :: alltimes(:)
real(kind_real),allocatable :: readbuf(:,:)

! Open NetCDF file
call ncerr(nf90_open(trim(self%filein),nf90_nowrite,ncid))

! Get groups ids
call ncerr(nf90_inq_grps(ncid,self%ngrp,grp_ids))

do igrp=1,self%ngrp
  ! Allocation
  if (igrp==1) then
    allocate(self%grphead)
    jgrp => self%grphead
  else
    allocate(jgrp%next)
    jgrp => jgrp%next
  endif

  ! Get group name
  call ncerr(nf90_inq_grpname(grp_ids(igrp),jgrp%grpname))

  ! Get dimension id
  call ncerr(nf90_inq_dimid(grp_ids(igrp),'nobs',nobs_id))

  ! Get dimension
  call ncerr(nf90_inquire_dimension(grp_ids(igrp),nobs_id,len=nobs_in_grp))

  ! Get variable id
  call ncerr(nf90_inq_varid(grp_ids(igrp),'times',times_id))

  ! Allocation
  allocate(inwindow(nobs_in_grp))
  allocate(alltimes(nobs_in_grp))

  ! Read in times
  jgrp%nobs = 0
  do iobs=1,nobs_in_grp
    call ncerr(nf90_get_var(grp_ids(igrp),times_id,stime,(/1,iobs/),(/50,1/)))
    call datetime_create(stime,tobs)
    if ((tobs > winbgn).and.(tobs <= winend)) then
      inwindow(iobs) = .true.
      alltimes(iobs) = tobs
      jgrp%nobs = jgrp%nobs+1
    else
      inwindow(iobs) = .false.
    endif
    call datetime_delete(tobs)
  enddo

  ! Allocation
  allocate(jgrp%times(jgrp%nobs))

  ! Copy times
  jobs=0
  do iobs=1,nobs_in_grp
    if (inwindow(iobs)) then
      jobs = jobs+1
      jgrp%times(jobs) = alltimes(iobs)
    endif
  end do

  ! Count columns
  call ncerr(nf90_inq_grps(grp_ids(igrp),ncol,col_ids))

  ! Loop over columns
  do icol=1,ncol
    ! Allocation
    if (icol==1) then
      allocate(jgrp%colhead)
      jcol => jgrp%colhead
    else
      allocate(jcol%next)
      jcol => jcol%next
    endif

    ! Get column name
    call ncerr(nf90_inq_grpname(col_ids(icol),jcol%colname))

    ! Get dimension id
    call ncerr(nf90_inq_dimid(col_ids(icol),'nlev',nlev_id))

    ! Get dimension
    call ncerr(nf90_inquire_dimension(col_ids(icol),nlev_id,len=jcol%nlev))

    ! Get variable id
    call ncerr(nf90_inq_varid(col_ids(icol),'values',values_id))

    ! Allocation
    allocate(readbuf(jcol%nlev,nobs_in_grp))
    allocate(jcol%values(jcol%nlev,jgrp%nobs))

    ! Get values
    call ncerr(nf90_get_var(col_ids(icol),values_id,readbuf(1:jcol%nlev,:),(/1,1/),(/jcol%nlev,nobs_in_grp/)))

    ! Copy values
    jobs = 0
    do iobs=1,nobs_in_grp
      if (inwindow(iobs)) then
        jobs = jobs+1
        jcol%values(:,jobs) = readbuf(:,iobs)
      endif
    enddo

    ! Release memory
    deallocate(readbuf)
  enddo

  ! Release memory
  do iobs=1,nobs_in_grp
    if (inwindow(iobs)) call datetime_delete(alltimes(iobs))
  enddo
  deallocate(alltimes)
  deallocate(inwindow)
enddo

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_obsdb_read
! ------------------------------------------------------------------------------
!> Write observation data
subroutine qg_obsdb_write(self)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self !< Observation data

! Local variables
integer :: iobs
integer :: ncid,nstrmax_id,grp_id,nobs_id,times_id,col_id,nlev_id,values_id
type(group_data),pointer :: jgrp
type(column_data),pointer :: jcol
character(len=50) :: stime

! Create NetCDF file
call ncerr(nf90_create(trim(self%fileout),ior(nf90_clobber,nf90_netcdf4),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nstrmax',50,nstrmax_id))

! Loop over groups
jgrp => self%grphead
do while (associated(jgrp))
  if (jgrp%nobs > 0) then
    ! Create group
    call ncerr(nf90_def_grp(ncid,jgrp%grpname,grp_id))

    ! Define dimension
    call ncerr(nf90_def_dim(grp_id,'nobs',jgrp%nobs,nobs_id))

    ! Define variable
    call ncerr(nf90_def_var(grp_id,'times',nf90_char,(/nstrmax_id,nobs_id/),times_id))

    ! Put variable
    do iobs=1,jgrp%nobs
      call datetime_to_string(jgrp%times(iobs),stime)
      call ncerr(nf90_put_var(grp_id,times_id,stime,(/1,iobs/),(/50,1/)))
    end do

    ! Loop over columns
    jcol => jgrp%colhead
    do while (associated(jcol))
      ! Create subgroup
      call ncerr(nf90_def_grp(grp_id,jcol%colname,col_id))

      ! Define dimension
      call ncerr(nf90_def_dim(col_id,'nlev',jcol%nlev,nlev_id))

      ! Define variable
      call ncerr(nf90_def_var(col_id,'values',nf90_double,(/nlev_id,nobs_id/),values_id))

      ! Put variable
      call ncerr(nf90_put_var(col_id,values_id,jcol%values(1:jcol%nlev,:),(/1,1/),(/jcol%nlev,jgrp%nobs/)))

      ! Update
      jcol => jcol%next
    enddo
  endif

  ! Update
  jgrp=>jgrp%next
end do

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_obsdb_write
! ------------------------------------------------------------------------------
!> Find observation data group
subroutine qg_obsdb_find_group(self,grp,find)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: self              !< Observation data
character(len=*),intent(in) :: grp             !< Group
type(group_data),pointer,intent(inout) :: find !< Result

! Initialization
find => self%grphead

! Loop
do while (associated(find))
  if (find%grpname==grp) exit
  find => find%next
enddo

end subroutine qg_obsdb_find_group
! ------------------------------------------------------------------------------
!> Find observation data column
subroutine qg_obsdb_find_column(grp,col,find)

implicit none

! Passed variables
type(group_data),intent(in) :: grp              !< Observation data
character(len=*),intent(in) :: col              !< Column
type(column_data),pointer,intent(inout) :: find !< Result

! Initialization
find=>grp%colhead

! Loop
do while (associated(find))
  if (find%colname==col) exit
  find => find%next
enddo

end subroutine qg_obsdb_find_column
! ------------------------------------------------------------------------------
!> Generate random locations
subroutine qg_obsdb_generate_locations(nlocs,lon,lat,z,ntimes,bgn,step,times,obsloc)

implicit none

! Passed variables
integer,intent(in) :: nlocs                         !< Number of locations
real(kind_real),intent(in) :: lon(nlocs),lat(nlocs),z(nlocs) !< Locations
integer,intent(in) :: ntimes                        !< Number of time-slots
type(datetime),intent(in) :: bgn                    !< Start time
type(duration),intent(in) :: step                   !< Time-step
type(datetime),intent(inout) :: times(nlocs*ntimes) !< Time-slots
type(qg_obsvec),intent(inout) :: obsloc             !< Observation locations

! Local variables
integer :: jobs,iobs,jstep
type(datetime) :: now

! Setup observation vector
call qg_obsvec_setup(obsloc,3,nlocs*ntimes)

! Set observation locations
now = bgn
iobs=0
do jstep=1,ntimes
  do jobs=1,nlocs
    iobs = iobs+1
    times(iobs) = now
    obsloc%values(:,iobs) = (/lon(jobs),lat(jobs),z(jobs)/)
  enddo
  call datetime_update(now,step)
enddo

! Release memory
call datetime_delete(now)

end subroutine qg_obsdb_generate_locations
! ------------------------------------------------------------------------------
!> Create observation data
subroutine qg_obsdb_create(self,grp,times,locs)

implicit none

! Passed varaibles
type(qg_obsdb),intent(inout) :: self  !< Observation data
character(len=*),intent(in) :: grp    !< Group
type(datetime),intent(in) :: times(:) !< Time-slots
type(qg_obsvec),intent(in) :: locs    !< Locations

! Local variables
type(group_data),pointer :: igrp
integer :: jobs,jlev

! Find observation group
call qg_obsdb_find_group(self,grp,igrp)
if (associated(igrp)) call abor1_ftn('qg_obsdb_create: obs group already exists')
if (associated(self%grphead)) then
  igrp => self%grphead
  do while (associated(igrp%next))
    igrp => igrp%next
  enddo
  allocate(igrp%next)
  igrp => igrp%next
else
  allocate(self%grphead)
  igrp => self%grphead
endif

! Create observation data
igrp%grpname = grp
igrp%nobs = size(times)
allocate(igrp%times(igrp%nobs))
igrp%times(:) = times(:)
allocate(igrp%colhead)
igrp%colhead%colname = 'Location'
igrp%colhead%nlev = 3
allocate(igrp%colhead%values(3,igrp%nobs))
if (locs%nlev/=3) call abor1_ftn('qg_obsdb_create: error locations not 3D')
if (locs%nobs/=igrp%nobs) call abor1_ftn('qg_obsdb_create: error locations number')
do jobs=1,igrp%nobs
  do jlev=1,3
    igrp%colhead%values(jlev,jobs) = locs%values(jlev,jobs)
  enddo
enddo
self%ngrp = self%ngrp+1

end subroutine qg_obsdb_create
! ------------------------------------------------------------------------------
end module qg_obsdb_mod
