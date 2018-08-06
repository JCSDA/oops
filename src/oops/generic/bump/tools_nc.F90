!----------------------------------------------------------------------
! Module: tools_nc
!> Purpose: NetCDF routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module tools_nc

use netcdf
use tools_kinds, only: kind_real
use type_mpl, only: mpl_type

implicit none

interface put_att
  module procedure put_att_integer
  module procedure put_att_integer_array
  module procedure put_att_real
  module procedure put_att_real_array
  module procedure put_att_logical
  module procedure put_att_logical_array
  module procedure put_att_string
  module procedure put_att_string_array
end interface

integer,parameter :: ncfloat = nf90_double !< NetCDF type for real

private
public :: ncfloat
public :: put_att

contains

!----------------------------------------------------------------------
! Subroutine: put_att_integer
!> Purpose: write namelist integer as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_integer(mpl,ncid,varname,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: var              !< Integer

! Local variables
character(len=1024) :: subr='put_att_integer'

! Write integer
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))

end subroutine put_att_integer

!----------------------------------------------------------------------
! Subroutine: put_att_integer_array
!> Purpose: write namelist integer array as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_integer_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< Integer array size
integer,intent(in) :: var(n)           !< Integer array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024) :: subr='put_att_integer_array'

! Write integer array as a string
if (n>0) then
   write(fullstr,'(i3.3)') var(1)
   do i=2,n
      write(str,'(i3.3)') var(i)
      fullstr = trim(fullstr)//':'//trim(str)
   end do
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine put_att_integer_array

!----------------------------------------------------------------------
! Subroutine: put_att_real
!> Purpose: write namelist real as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_real(mpl,ncid,varname,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
real(kind_real),intent(in) :: var      !< Real

! Local variables
character(len=1024) :: subr='put_att_real'

! Write real
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),var))

end subroutine put_att_real

!----------------------------------------------------------------------
! Subroutine: put_att_real_array
!> Purpose: write namelist real array as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_real_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< Real array size
real(kind_real),intent(in) :: var(n)   !< Real array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024) :: subr='put_att_real_array'

! Write real array as a string
if (n>0) then
   write(fullstr,'(e10.3)') var(1)
   do i=2,n
      write(str,'(e10.3)') var(i)
      fullstr = trim(fullstr)//':'//trim(str)
   end do
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine put_att_real_array

!----------------------------------------------------------------------
! Subroutine: put_att_logical
!> Purpose: write namelist logical as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_logical(mpl,ncid,varname,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
logical,intent(in) :: var              !< Logical

! Local variables
character(len=1024) :: subr='put_att_logical'

! Write logical as a string
if (var) then
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.true.'))
else
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),'.false.'))
end if

end subroutine put_att_logical

!----------------------------------------------------------------------
! Subroutine: put_att_logical_array
!> Purpose: write namelist logical array as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_logical_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< Real array size
logical,intent(in) :: var(n)           !< Logical array

! Local variables
integer :: i
character(len=1024) :: str,fullstr
character(len=1024) :: subr='put_att_logical_array'

! Write real array as a string
if (n>0) then
   if (var(1)) then
      write(fullstr,'(a6)') '.true.'
   else
      write(fullstr,'(a7)') '.false.'
   end if
   do i=2,n
      if (var(i)) then
         write(str,'(a6)') '.true.'
      else
         write(str,'(a7)') '.false.'
      end if
      fullstr = trim(fullstr)//':'//trim(str)
   end do
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine put_att_logical_array

!----------------------------------------------------------------------
! Subroutine: put_att_string
!> Purpose: write namelist string as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_string(mpl,ncid,varname,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
character(len=*),intent(in) :: var     !< String

! Local variables
character(len=1024) :: subr='put_att_string'

! Write string
call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(var)))

end subroutine put_att_string

!----------------------------------------------------------------------
! Subroutine: put_att_string_array
!> Purpose: write namelist string array as NetCDF attribute
!----------------------------------------------------------------------
subroutine put_att_string_array(mpl,ncid,varname,n,var)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
integer,intent(in) :: ncid             !< NetCDF file id
character(len=*),intent(in) :: varname !< Variable name
integer,intent(in) :: n                !< String array size
character(len=*),intent(in) :: var(n)  !< String array

! Local variables
integer :: i
character(len=1024) :: fullstr
character(len=1024) :: subr='put_att_string_array'

! Write string array
if (n>0) then
   fullstr = trim(var(1))
   do i=2,n
      fullstr = trim(fullstr)//':'//trim(var(i))
   end do
   call mpl%ncerr(subr,nf90_put_att(ncid,nf90_global,trim(varname),trim(fullstr)))
end if

end subroutine put_att_string_array

end module tools_nc
