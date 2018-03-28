!----------------------------------------------------------------------
! Module: hdiag_transform.f90
!> Purpose: transform routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module hdiag_transform

use tools_display, only: msgerror
use tools_jacobi_eigenvalue, only: jacobi_eigenvalue
use tools_kinds, only: kind_real
use type_avg, only: avgtype
use type_mpl, only: mpl
use type_hdata, only: hdatatype
implicit none

real(kind_real),parameter :: rvth = 0.5       !< Theoretical vertical support radius
real(kind_real),parameter :: egvmin = 1.0e-12 !< Minimum eigenvalue

private
public :: compute_transform

contains

!----------------------------------------------------------------------
! Subroutine: compute_transform
!> Purpose: compute transform
!----------------------------------------------------------------------
subroutine compute_transform(hdata,ib,avg,trans,transinv)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                                    !< HDIAG data
integer,intent(in) :: ib                                               !< Block index
type(avgtype),intent(in) :: avg                                        !< Averaged statistics
real(kind_real),intent(out) :: trans(hdata%geom%nl0,hdata%geom%nl0)    !< Direct transform
real(kind_real),intent(out) :: transinv(hdata%geom%nl0,hdata%geom%nl0) !< Inverse transform

! Local variables
integer :: jl0r,jl0,il0,it_num,rot_num
real(kind_real) :: cor(hdata%geom%nl0,hdata%geom%nl0),corth(hdata%geom%nl0,hdata%geom%nl0)
real(kind_real) :: corsqrt(hdata%geom%nl0,hdata%geom%nl0),corsqrtinv(hdata%geom%nl0,hdata%geom%nl0)
real(kind_real) :: corthsqrt(hdata%geom%nl0,hdata%geom%nl0),corthsqrtinv(hdata%geom%nl0,hdata%geom%nl0)
real(kind_real) :: v(hdata%geom%nl0,hdata%geom%nl0),d(hdata%geom%nl0)
real(kind_real) :: dd(hdata%geom%nl0,hdata%geom%nl0),ddinv(hdata%geom%nl0,hdata%geom%nl0)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Copy correlation
cor = 0.0
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
      cor(jl0,il0) = avg%cor(1,jl0r,il0)
   end do
end do

! Compute eigenvalues
call jacobi_eigenvalue(geom%nl0,cor,500,v,d,it_num,rot_num)

! Eigenvalues thresholding
d = max(d,egvmin)

! Inverse correlation square-root
dd = 0.0
ddinv = 0.0
do il0=1,geom%nl0
   dd(il0,il0) = sqrt(d(il0))
   ddinv(il0,il0) = 1.0/sqrt(d(il0))
end do
corsqrt = matmul(v,matmul(dd,transpose(v)))
corsqrtinv = matmul(v,matmul(ddinv,transpose(v)))

! Theoretical correlation
do il0=1,geom%nl0
   do jl0r=1,bpar%nl0r(ib)
      jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
      corth(jl0,il0) = exp(-0.5*(geom%distv(jl0,il0)/rvth)**2)
   end do
end do

! Compute eigenvalues
call jacobi_eigenvalue(geom%nl0,corth,500,v,d,it_num,rot_num)

! Eigenvalues thresholding
d = max(d,egvmin)

! Inverse correlation square-root
dd = 0.0
ddinv = 0.0
do il0=1,geom%nl0
   dd(il0,il0) = sqrt(d(il0))
   ddinv(il0,il0) = 1.0/sqrt(d(il0))
end do
corthsqrt = matmul(v,matmul(dd,transpose(v)))
corthsqrtinv = matmul(v,matmul(ddinv,transpose(v)))

! Compute direct and inverse transform
trans = matmul(corsqrt,corthsqrtinv)
transinv = matmul(corthsqrt,corsqrtinv)

! Print results
write(mpl%unit,'(a7,a,e15.8,a,e15.8)') '','Direct transform (min/max):  ',minval(trans),' / ',maxval(trans)
write(mpl%unit,'(a7,a,e15.8,a,e15.8)') '','Inverse transform (min/max): ',minval(transinv),' / ',maxval(transinv)

! End associate
end associate

end subroutine compute_transform

end module hdiag_transform
