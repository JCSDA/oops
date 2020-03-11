! (C) Copyright 2019-2020 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module bump_interpolation_mod

  use kinds
  use type_bump, only : bump_type

  implicit none
  private
  
  integer, parameter :: max_string = 400

  ! ------------------------------------------------------------------------------
  !> Unstructured latlon grid
  !! This defines the grid as a 1D vector of latitudes and a 1D vector of longitudes

  type, public :: bump_grid
     private
     integer :: ngrid
     real(kind_real), allocatable :: lat(:)
     real(kind_real), allocatable :: lon(:)
   contains
     private
     procedure :: latlon_2d_to_1d_
     generic, public :: create => latlon_2d_to_1d_
     procedure, public :: delete => delete_bump_grid_
  end type

  ! ------------------------------------------------------------------------------
  !> one or more fields on an unstructured grid 

  type, public :: bump_field
     private
     integer :: ngrid
     real(kind_real), allocatable :: data(:)
   contains
     private
     !< make these generic as other interfaces are added
     procedure, public :: create => create_fields_2d_to_1d_ 
     procedure, public :: result => result_fields_1d_to_2d_
     procedure, public :: delete => delete_bump_field_
  end type

  ! ------------------------------------------------------------------------------
  !!
  !! bint = shorthand for bump interpolator
  !!
  type, public :: bump_interpolator  
  private
    type(bump_type) :: bump
  contains
    private

    ! after we implement bump, we can add another init and apply interface
    ! that takes 1D unstructured grids as input and then passes them directly
    ! to bump, thus avoiding a copy
    procedure :: bint_grid_init_
    procedure :: bint_latlon_init_
    generic, public :: init => bint_grid_init_, bint_latlon_init_
    procedure :: bint_grid_apply_
    procedure :: bint_latlon_apply_
    generic, public :: apply => bint_grid_apply_, bint_latlon_apply_
    procedure, public :: delete => bint_delete_
    
 end type bump_interpolator

contains

! -----------------------------------------------------------------------------
subroutine latlon_2d_to_1d_(self, lat, lon)
  use kinds
  class(bump_grid), intent(inout) :: self
  real(kind_real), intent(in) :: lat(:,:), lon(:,:)

  ! for now do a copy.  But, if efficiency is a concern we can avoid a copy in
  ! the future by declaring self%lat and self%lon as pointers and setting them
  ! as follows:
  ! self%lat => lat(1,1)
  ! self%lon => lon(1,1)
 
  
  if (allocated(self%lat)) deallocate(self%lat)
  if (allocated(self%lon)) deallocate(self%lon)

  self%ngrid = size(lat)
  allocate(self%lat(self%ngrid),self%lon(self%ngrid))
  self%lat = reshape(lat, (/ self%ngrid /))
  self%lon = reshape(lon, (/ self%ngrid /))

end subroutine latlon_2d_to_1d_

! -----------------------------------------------------------------------------
!> this is an interface for a single field
!! we can add an interface for multiple fields if needed
subroutine create_fields_2d_to_1d_(self, infield, allocate_only)
  use kinds
  class(bump_field), intent(inout) :: self
  real(kind_real), intent(in) :: infield(:,:)
  logical, optional, intent(in) :: allocate_only
    
  ! for now do a copy.  But, if efficiency is a concern we can avoid a copy in
  ! the future by declaring self%data as a pointers and setting it as follows:
  ! self%data => infield(1,1)

  if (allocated(self%data)) deallocate(self%data)
    
  self%ngrid = size(infield)
  allocate(self%data(self%ngrid))

  if (present(allocate_only)) then
     if (allocate_only) return
  endif
    
  self%data = reshape(infield, (/ self%ngrid /))

end subroutine create_fields_2d_to_1d_

! ------------------------------------------------------------------------------
subroutine result_fields_1d_to_2d_(self, outfield)
  use kinds
  class(bump_field), intent(inout) :: self
  real(kind_real), intent(inout) :: outfield(:,:)
  character(len=max_string) :: err_msg
  integer nx, ny

  if (size(outfield) .lt. self%ngrid) then
     write(err_msg,'(a)') &
          "oops::bump_interpolation_mod::result_fields_1d_to_2d: "// &
          "Output field is not properly allocated"
     call abor1_ftn(trim(err_msg))
  endif

  nx = size(outfield,1)
  ny = size(outfield,2)
  outfield = reshape(self%data, (/ nx, ny /), order = (/ 1, 2 /) )
    
end subroutine result_fields_1d_to_2d_

! -----------------------------------------------------------------------------
subroutine delete_bump_grid_(self)
  class(bump_grid), intent(inout) :: self

  self%ngrid=0
  if (allocated(self%lat)) deallocate(self%lat)
  if (allocated(self%lon)) deallocate(self%lon)
    
end subroutine delete_bump_grid_

! -----------------------------------------------------------------------------
subroutine delete_bump_field_(self)
  class(bump_field), intent(inout) :: self

  self%ngrid=0
  if (allocated(self%data)) deallocate(self%data)
    
end subroutine delete_bump_field_
    
! ------------------------------------------------------------------------------
!> init that takes bump_grid objects as input
!!
  
subroutine bint_grid_init_(self, comm, config, in_grid, out_grid)
  use fckit_configuration_module, only : fckit_configuration
  use fckit_mpi_module, only: fckit_mpi_comm
  use iso_c_binding, only : c_char
  class(bump_interpolator), intent(inout) :: self
  class(bump_grid), intent(in)  :: in_grid, out_grid
  class(fckit_configuration), intent(in) :: config
  class(fckit_mpi_comm), intent(in) :: comm
  character(kind=c_char,len=:), allocatable :: string_buffer

  call self%bump%nam%init(comm%size())

  if (config%has("prefix")) then
     call config%get_or_die("prefix",string_buffer)
     self%bump%nam%prefix = string_buffer
     if (allocated(string_buffer)) deallocate(string_buffer)
  else
     self%bump%nam%prefix = 'bump_data'
  endif

  ! include other namelist settings next - the above is an example
  
  ! here is where the mean of the bump init() method goes
  write(*,*) "Main init routine"  

end subroutine bint_grid_init_

! -----------------------------------------------------------------------------
!> latlon init
!!
!! init option that takes input and output grids as 2D latitude and
!! longitude arrays
subroutine bint_latlon_init_(self, comm, config, lat_in, lon_in, lat_out, lon_out)
  use fckit_configuration_module, only : fckit_configuration
  use fckit_mpi_module, only : fckit_mpi_comm
  class(bump_interpolator), intent(inout) :: self
  class(fckit_configuration), intent(in) :: config
  class(fckit_mpi_comm), intent(in) :: comm
  real(kind_real), intent(in) :: lat_in(:,:), lon_in(:,:)
  real(kind_real), intent(in) :: lat_out(:,:), lon_out(:,:)
  type(bump_grid) :: in_grid, out_grid

  call in_grid%create(lat_in, lon_in)
  call out_grid%create(lat_out, lon_out)

  call self%init(comm, config, in_grid, out_grid)

  call in_grid%delete()
  call out_grid%delete()

end subroutine bint_latlon_init_
  
  ! ------------------------------------------------------------------------------
  !> apply to fields on a 2D latlon grid
  !!
  !! This assumes the input and output fields are defined on 2D lat, lon
  !! grids as used in bint_latlon_init_.
  !!
  !! It is assumedd that out_fields is allocated by the client before this
  !! routine is called.
  !!

  subroutine bint_latlon_apply_(self, in_fields, out_fields)
    class(bump_interpolator), intent(inout) :: self
    real(kind_real), intent(in)    :: in_fields(:,:)
    real(kind_real), intent(inout) :: out_fields(:,:)
    type(bump_field) :: fields_in, fields_out

     call fields_in%create(in_fields)
     call fields_out%create(out_fields, allocate_only=.true.)
 
     call self%apply(fields_in, fields_out)

     call fields_out%result(out_fields)

     call fields_in%delete()
     call fields_out%delete()

  end subroutine bint_latlon_apply_

! ------------------------------------------------------------------------------
subroutine bint_grid_apply_(self, in_fields, out_fields)
  class(bump_interpolator), intent(inout) :: self
  class(bump_field), intent(in)     :: in_fields
  class(bump_field), intent(inout)  :: out_fields

  ! here is where the meat of the bump apply() method goes
  write(*,*) "Main apply routine"

end subroutine bint_grid_apply_

! ------------------------------------------------------------------------------

subroutine bint_delete_(self)
  class(bump_interpolator), intent(inout) :: self

  call self%bump%dealloc()
  
end subroutine bint_delete_
   
end module bump_interpolation_mod
