! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module unstructured_interpolation_mod

use atlas_module, only: atlas_geometry, atlas_indexkdtree

use netcdf

use kinds
use string_utils, only: replace_string
use netcdf_utils_mod, only: nccheck

use fckit_mpi_module,      only: fckit_mpi_comm, fckit_mpi_sum

implicit none
private
public unstrc_interp
public unstrc_interp_registry

!---------------------------------------------------------------------------------------------------

type unstrc_interp
  integer :: nn                                       ! Number of neighbours
  integer :: ngrid_in                                 ! Number of input grid points
  integer :: ngrid_out                                ! Number of output grid points
  integer :: ngrid_in_glo                             ! Number of global input grid points
  type(fckit_mpi_comm) :: comm                        ! Communicator
  integer, allocatable :: displs(:)                   ! Displacement for global gather
  integer, allocatable :: rcvcnt(:)                   ! Receive count for global gather
  real(kind=kind_real), allocatable :: interp_w(:,:)  ! Interpolation weights
  integer             , allocatable :: interp_i(:,:)  ! Interpolation index
  contains
    generic, public :: create =>  create_new, create_read
    procedure, private :: create_new, create_read
    procedure, public :: delete
    procedure, public :: apply
    procedure, public :: apply_ad
    procedure, public :: write
endtype unstrc_interp

! ------------------------------------------------------------------------------
!> Registry for unstrc_interpo objects

#define LISTED_TYPE unstrc_interp

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: unstrc_interp_registry

!---------------------------------------------------------------------------------------------------

contains
!-------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

!---------------------------------------------------------------------------------------------------

subroutine create_new( self, comm, nn, wtype, ngrid_in,  lats_in,  lons_in, &
                                              ngrid_out, lats_out, lons_out )
class(unstrc_interp), intent(inout) :: self                ! Myself
type(fckit_mpi_comm), intent(in)    :: comm                ! Communicator
integer,              intent(in)    :: nn                  ! Number of neighbours
character(len=*),     intent(in)    :: wtype               ! Weight types to return
integer,              intent(in)    :: ngrid_in            ! Number of grid points on input grid
real(kind=kind_real), intent(in)    :: lats_in(ngrid_in)   ! Input grid latitude in degrees
real(kind=kind_real), intent(in)    :: lons_in(ngrid_in)   ! Input grid longitude in degrees
integer,              intent(in)    :: ngrid_out           ! Number of grid points on output grid
real(kind=kind_real), intent(in)    :: lats_out(ngrid_out) ! Output grid latitide in degrees
real(kind=kind_real), intent(in)    :: lons_out(ngrid_out) ! Output grid longitude in degrees

!Locals
type(atlas_indexkdtree) :: kd
type(atlas_geometry) :: ageometry
integer :: j, n, jj, kk, ngrid_in_glo
integer :: index1, index2, nindex
real(kind=kind_real) :: dist, wprod, bsw
real(kind=kind_real) :: lons_in_loc(ngrid_in), lons_out_loc(ngrid_out)
real(kind=kind_real), allocatable :: lats_in_glo(:), lons_in_glo(:)
real(kind=kind_real), allocatable :: nn_dist(:,:)
real(kind=kind_real), allocatable :: bw(:)

! wtype options
! -------------
! 1. (none)     No weights needed
! 2. (distance) Inverse of distance to nearest neighbours
! 3. (barycent) Barycentric weights

! Allocate self type
self%nn = nn
self%ngrid_in = ngrid_in
self%ngrid_out = ngrid_out
self%comm = comm
allocate(self%interp_w(self%nn,self%ngrid_out))
allocate(self%interp_i(self%nn,self%ngrid_out))

! Get input grid proc counts and displacement
call input_grid_share(self%comm, self%ngrid_in, self%rcvcnt, self%displs)

! Global number of grid points

! ----------------------------
ngrid_in_glo = sum(self%rcvcnt)
self%ngrid_in_glo = ngrid_in_glo

! Diplacement for each processor
! ------------------------------
self%displs(1) = 0
do j = 2,comm%size()
   self%displs(j) = self%displs(j-1) + self%rcvcnt(j-1)
enddo

! Rotate longitudes less than 0
lons_in_loc  = lons_in
lons_out_loc = lons_out
where (lons_in_loc  < 0.0) lons_in_loc  = lons_in_loc  + 360.0_kind_real
where (lons_out_loc < 0.0) lons_out_loc = lons_out_loc + 360.0_kind_real

! Gather the lat/lons to all
allocate(lats_in_glo(ngrid_in_glo))
allocate(lons_in_glo(ngrid_in_glo))
call comm%allgather(lats_in    ,lats_in_glo,ngrid_in,self%rcvcnt,self%displs)
call comm%allgather(lons_in_loc,lons_in_glo,ngrid_in,self%rcvcnt,self%displs)

! Create UnitSphere geometry
ageometry = atlas_geometry("UnitSphere")

! Create KDTree
kd = atlas_indexkdtree(ageometry)
call kd%reserve(ngrid_in_glo)
call kd%build(ngrid_in_glo,lons_in_glo,lats_in_glo)

! Loop over observations
! ----------------------
allocate(nn_dist(nn,ngrid_out))

do n = 1,ngrid_out

  ! Get nearest neighbours
  call kd%closestPoints(lons_out_loc(n), lats_out(n), self%nn, self%interp_i(:,n))

  ! Compute distances
  do kk = 1, nn
    nindex = self%interp_i(kk,n)
    nn_dist(kk,n) = ageometry%distance(lons_out_loc(n), lats_out(n), lons_in_glo(nindex), lats_in_glo(nindex))
  enddo

enddo


! Set weights based on user choice
! --------------------------------
select case (wtype)
  case ('none')
    self%interp_w = 0.0_kind_real

  case ('distance')
    do n = 1,ngrid_out
      bsw = 0.0_kind_real
      do jj = 1,nn
        bsw = bsw + nn_dist(jj,n)
      enddo
      self%interp_w(:,n) = nn_dist(:,n) / bsw
    enddo

  case ('barycent')
    allocate(bw(self%nn))
    do n = 1,ngrid_out
      !Barycentric weights formula
      bw(:) = 0.0_kind_real
      do jj = 1, nn
        wprod = 1.0_kind_real
        index1 = self%interp_i(jj,n)
        do kk = 1, nn
          if (jj.ne.kk) then
            index2 = self%interp_i(kk,n)
            dist = ageometry%distance(lons_in_glo(index1),lats_in_glo(index1),&
                                      lons_in_glo(index2),lats_in_glo(index2))
            wprod = wprod * max(dist, 1e-10)
          endif
        enddo
        bw(jj) = 1.0_kind_real / wprod
      enddo

      !Barycentric weights
      self%interp_w(:,n) = 0.0_kind_real
      if (minval(nn_dist(:,n)) < 1e-10) then

        ! special case if very close to one grid point
        jj = minloc(nn_dist(:,n),dim=1)
        self%interp_w(jj,n) = 1.0_kind_real

      else

        !otherwise continue with the normal algorithm
        bsw = 0.0_kind_real
        do jj = 1,nn
          bsw = bsw + (bw(jj) / nn_dist(jj,n))
        enddo

        do jj = 1,nn
          self%interp_w(jj,n) = ( bw(jj) / nn_dist(jj,n) ) / bsw
        enddo
      end if

    enddo

    deallocate(bw)

  case default
    call abor1_ftn("unstrc_interpolation.create: wrong choice of weight output")

end select

!Deallocate
call kd%final()
deallocate(nn_dist)
deallocate(lats_in_glo, lons_in_glo)

end subroutine create_new

!---------------------------------------------------------------------------------------------------

subroutine create_read(self, comm, filename_in)
class(unstrc_interp), intent(inout) :: self
type(fckit_mpi_comm), intent(in)    :: comm
character(len=*),     intent(in)    :: filename_in

character(len=2048) :: filename
integer :: ncid, varid

self%comm = comm

! One file per processor
filename = filename_in
call getfilename(self%comm%rank(),filename)

! Open file
call nccheck ( nf90_open( trim(filename), NF90_NOWRITE, ncid ), "nf90_open"//trim(filename) )

! Get the dimensions
! ------------------
call nccheck ( nf90_inq_dimid(ncid, "ngrid_in", varid), "nf90_inq_dimid ngrid_in" )
call nccheck ( nf90_inquire_dimension(ncid, varid, len = self%ngrid_in), "nf90_inquire_dimension ngrid_in" )

call nccheck ( nf90_inq_dimid(ncid, "ngrid_out", varid), "nf90_inq_dimid ngrid_out" )
call nccheck ( nf90_inquire_dimension(ncid, varid, len = self%ngrid_out), "nf90_inquire_dimension ngrid_out" )

call nccheck ( nf90_inq_dimid(ncid, "nn", varid), "nf90_inq_dimid nn" )
call nccheck ( nf90_inquire_dimension(ncid, varid, len = self%nn), "nf90_inquire_dimension nn" )

! Allocate arrays
allocate(self%interp_w(self%nn,self%ngrid_out))
allocate(self%interp_i(self%nn,self%ngrid_out))

! Get the interpolation weights and indices
! ----------------------------------------
call nccheck ( nf90_inq_varid (ncid, "interp_weights", varid), "nf90_inq_varid interp_weights" )
call nccheck ( nf90_get_var   (ncid, varid, self%interp_w), "nf90_get_var interp_weights" )

call nccheck ( nf90_inq_varid (ncid, "interp_indices", varid), "nf90_inq_varid interp_indices" )
call nccheck ( nf90_get_var   (ncid, varid, self%interp_i), "nf90_get_var interp_indices" )

! Close the file
call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine create_read

!---------------------------------------------------------------------------------------------------

subroutine delete( self )
class(unstrc_interp), intent(inout) :: self

if (allocated(self%interp_w)) deallocate(self%interp_w)
if (allocated(self%interp_i)) deallocate(self%interp_i)
if (allocated(self%rcvcnt)) deallocate(self%rcvcnt)
if (allocated(self%displs)) deallocate(self%displs)

end subroutine delete

!---------------------------------------------------------------------------------------------------

subroutine apply(self, field_in, field_out, field_nn_out )
class(unstrc_interp),           intent(in)  :: self                                 ! Myself
real(kind=kind_real),           intent(in)  :: field_in(self%ngrid_in)              ! Input field
real(kind=kind_real),           intent(out) :: field_out(self%ngrid_out)            ! Result of interpolation
real(kind=kind_real), optional, intent(out) :: field_nn_out(self%nn,self%ngrid_out) ! Neighbours

!Locals
integer :: n, kk
real(kind=kind_real), allocatable :: field_nn(:,:), field_in_glo(:)


! Gather the field
allocate(field_in_glo(sum(self%rcvcnt)))
call self%comm%allgather(field_in,field_in_glo,self%ngrid_in,self%rcvcnt,self%displs)

! Get output neighbours
allocate(field_nn(self%nn,self%ngrid_out))
do n = 1, self%ngrid_out
  do kk = 1, self%nn
    field_nn(kk,n) = field_in_glo(self%interp_i(kk,n))
  enddo
enddo

! Deallocate global field
deallocate(field_in_glo)

! Apply weights
! -------------
do n = 1, self%ngrid_out
  field_out(n) = 0.0_kind_real
  do kk = 1, self%nn
    field_out(n) = field_out(n) + self%interp_w(kk,n) * field_nn(kk,n)
  enddo
enddo

! Return neighbours if requested
if (present(field_nn_out)) field_nn_out = field_nn

deallocate(field_nn)

end subroutine apply

!---------------------------------------------------------------------------------------------------

subroutine apply_ad(self, field_in, field_out )
class(unstrc_interp),           intent(in) :: self
real(kind=kind_real),          intent(out) :: field_in(self%ngrid_in)
real(kind=kind_real), optional, intent(in) :: field_out(self%ngrid_out)

integer :: n, kk, nstart, nfinish
real(kind=kind_real), allocatable :: field_in_glo(:), field_in_glo_all(:)


! Apply backward interpolation
! ----------------------------
allocate(field_in_glo(self%ngrid_in_glo))
field_in_glo = 0.0
do n = 1, self%ngrid_out
  do kk = 1, self%nn
     field_in_glo(self%interp_i(kk,n)) = field_in_glo(self%interp_i(kk,n)) +  self%interp_w(kk,n)*field_out(n)
  enddo
enddo

! Global sum
! ----------
allocate(field_in_glo_all(self%ngrid_in_glo))
call self%comm%allreduce(field_in_glo, field_in_glo_all, fckit_mpi_sum())

! Return local field
! ------------------
nstart = self%displs(self%comm%rank()+1) + 1
nfinish = self%displs(self%comm%rank()+1) + self%rcvcnt(self%comm%rank()+1)
field_in(1:self%ngrid_in) = field_in_glo_all(nstart:nfinish)
deallocate(field_in_glo)
deallocate(field_in_glo_all)

end subroutine apply_ad

!---------------------------------------------------------------------------------------------------

subroutine write(self,filename_in)
class(unstrc_interp), intent(in) :: self
character(len=*),     intent(in) :: filename_in

integer :: ncid, ni_dimid, no_dimid, nn_dimid, vc(2)
character(len=2048) :: filename

! One file per processor
filename = filename_in
call getfilename(self%comm%rank(),filename)

! Create file
call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_CLOBBER), ncid), "nf90_create" )

! Define dimensions
call nccheck( nf90_def_dim(ncid, "ngrid_in",  self%ngrid_in,  ni_dimid), "nf90_def_dim ng"  )
call nccheck( nf90_def_dim(ncid, "ngrid_out", self%ngrid_out, no_dimid), "nf90_def_dim ng"  )
call nccheck( nf90_def_dim(ncid, "nn",        self%nn,        nn_dimid), "nf90_def_dim nn"  )

! Define variables
call nccheck( nf90_def_var(ncid, "interp_weights", NF90_DOUBLE, (/ nn_dimid, no_dimid /), vc(1)), "nf90_def_var interp_weights" );
call nccheck( nf90_def_var(ncid, "interp_indices", NF90_INT,    (/ nn_dimid, no_dimid /), vc(2)), "nf90_def_var interp_indices" );

! End define
call nccheck( nf90_enddef(ncid), "nf90_enddef" )

! Write variables
call nccheck( nf90_put_var( ncid, vc(1), self%interp_w ), "nf90_put_var interp_weights" );
call nccheck( nf90_put_var( ncid, vc(2), self%interp_I ), "nf90_put_var interp_indices" );

! Close file
call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine write

!---------------------------------------------------------------------------------------------------

subroutine input_grid_share(comm, ngrid_in, rcvcnt, displs)
type(fckit_mpi_comm), intent(in)    :: comm
integer,              intent(in)    :: ngrid_in
integer, allocatable, intent(inout) :: rcvcnt(:)
integer, allocatable, intent(inout) :: displs(:)

integer :: j

! Gather receive count from each processor
allocate(rcvcnt(comm%size()))
call comm%allgather(ngrid_in, rcvcnt)

! Displacement for each processor
allocate(displs(comm%size()))
displs(1) = 0
do j = 2,comm%size()
   displs(j) = displs(j-1) + rcvcnt(j-1)
enddo

end subroutine input_grid_share

!---------------------------------------------------------------------------------------------------

subroutine getfilename(procrank, filename)
integer,          intent(in)    :: procrank
character(len=*), intent(inout) :: filename

character(len=6) :: procrankstr

! Error check
if (procrank>999999) call abor1_ftn("unstructured_interpolation_mod not ready for more than 999999 processors")

! String processor number
write(procrankstr,'(I0.6)') procrank

! Get new filename
filename = replace_string(filename,'%PROC',procrankstr)

end subroutine getfilename

! --------------------------------------------------------------------------------------------------

end module unstructured_interpolation_mod
