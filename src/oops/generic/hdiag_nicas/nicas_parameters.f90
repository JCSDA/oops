!----------------------------------------------------------------------
! Module: nicas_parameters.f90
!> Purpose: compute NICAS parameters
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_parameters

use nicas_parameters_convol, only: compute_convol_network,compute_convol_distance
use nicas_parameters_interp, only: compute_interp_h,compute_interp_v,compute_interp_s
use nicas_parameters_mpi, only: compute_mpi_ab,compute_mpi_c
use nicas_parameters_normalization, only: compute_normalization
use nicas_parameters_sampling, only: compute_sampling
use tools_display, only: msgerror
use tools_kinds,only: kind_real
use type_bdata, only: bdatatype
use type_mpl, only: mpl
use type_ndata, only: ndatatype

implicit none

private
public :: compute_parameters

contains

!----------------------------------------------------------------------
! Subroutine: compute_parameters
!> Purpose: compute NICAS parameters
!----------------------------------------------------------------------
subroutine compute_parameters(bdata,ndata)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata    !< B data
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: il0i,il1

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Compute adaptive sampling
write(mpl%unit,'(a7,a)') '','Compute adaptive sampling'
call compute_sampling(bdata,ndata)

! Compute horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute horizontal interpolation data'
call compute_interp_h(ndata)

! Compute vertical interpolation data
write(mpl%unit,'(a7,a)') '','Compute vertical interpolation data'
call compute_interp_v(ndata)

! Compute subsampling horizontal interpolation data
write(mpl%unit,'(a7,a)') '','Compute subsampling horizontal interpolation data'
call compute_interp_s(ndata)

! Compute MPI distribution, halos A-B
write(mpl%unit,'(a7,a)') '','Compute MPI distribution, halos A-B'
call compute_mpi_ab(ndata)

! Compute convolution data
write(mpl%unit,'(a7,a)') '','Compute convolution data'
if (nam%network) then
   call compute_convol_network(bdata,ndata)
else
   call compute_convol_distance(bdata,ndata)
end if

! Compute MPI distribution, halo C
write(mpl%unit,'(a7,a)') '','Compute MPI distribution, halo C'
call compute_mpi_c(ndata)

! Compute normalization
write(mpl%unit,'(a7,a)') '','Compute normalization'
call compute_normalization(ndata)

! Print results
write(mpl%unit,'(a7,a,i4)') '','Parameters for processor #',mpl%myproc
write(mpl%unit,'(a10,a,i8)') '','nc0 =        ',geom%nc0
write(mpl%unit,'(a10,a,i8)') '','nc0a =       ',geom%nc0a
write(mpl%unit,'(a10,a,i8)') '','nl0 =        ',geom%nl0
write(mpl%unit,'(a10,a,i8)') '','nc1 =        ',ndata%nc1
write(mpl%unit,'(a10,a,i8)') '','nc1a =       ',ndata%nc1a
write(mpl%unit,'(a10,a,i8)') '','nc1b =       ',ndata%nc1b
write(mpl%unit,'(a10,a,i8)') '','nl1 =        ',ndata%nl1
do il1=1,ndata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','nc2(',il1,') =   ',ndata%nc2(il1)
end do
write(mpl%unit,'(a10,a,i8)') '','ns =         ',ndata%ns
write(mpl%unit,'(a10,a,i8)') '','nsa =        ',ndata%nsa
write(mpl%unit,'(a10,a,i8)') '','nsb =        ',ndata%nsb
write(mpl%unit,'(a10,a,i8)') '','nsc =        ',ndata%nsc
do il0i=1,geom%nl0i
   write(mpl%unit,'(a10,a,i3,a,i8)') '','h(',il0i,')%n_s = ',ndata%h(il0i)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','v%n_s =      ',ndata%v%n_s
do il1=1,ndata%nl1
   write(mpl%unit,'(a10,a,i3,a,i8)') '','s(',il1,')%n_s = ',ndata%s(il1)%n_s
end do
write(mpl%unit,'(a10,a,i8)') '','c%n_s =      ',ndata%c%n_s

! End associate
end associate

end subroutine compute_parameters

end module nicas_parameters
