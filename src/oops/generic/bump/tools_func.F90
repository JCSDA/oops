!----------------------------------------------------------------------
! Module: tools_func
!> Purpose: usual functions
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module tools_func

use tools_const, only: pi
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr

implicit none

real(kind_real),parameter :: Hmin = 1.0e-12_kind_real !< Minimum tensor diagonal value
integer,parameter :: M = 0                            !< Number of implicit itteration for the Matern function (Gaussian function if M = 0)
real(kind_real),parameter :: eta = 1.0e-9_kind_real   !< Small parameter for the Cholesky decomposition

private
public :: lonlatmod,sphere_dist,reduce_arc,vector_product,vector_triple_product,add,divide,fac,fit_diag,gc99,fit_lct,cholesky

contains

!----------------------------------------------------------------------
! Subroutine: lonlatmod
!> Purpose: set latitude between -pi/2 and pi/2 and longitude between -pi and pi
!----------------------------------------------------------------------
subroutine lonlatmod(lon,lat)

implicit none

! Passed variables
real(kind_real),intent(inout) :: lon !< Longitude
real(kind_real),intent(inout) :: lat !< Latitude

! Check latitude bounds
if (lat>0.5*pi) then
   lat = pi-lat
   lon = lon+pi
elseif (lat<-0.5*pi) then
   lat = -pi-lat
   lon = lon+pi
end if

! Check longitude bounds
if (lon>pi) then
   lon = lon-2.0*pi
elseif (lon<-pi) then
   lon = lon+2.0*pi
end if

end subroutine lonlatmod

!----------------------------------------------------------------------
! Function: sphere_dist
!> Purpose: compute the great-circle distance between two points
!----------------------------------------------------------------------
subroutine sphere_dist(lon_i,lat_i,lon_f,lat_f,dist)

implicit none

! Passed variable
real(kind_real),intent(in) :: lon_i !< Initial point longitude (radian)
real(kind_real),intent(in) :: lat_i !< Initial point latitude (radian)
real(kind_real),intent(in) :: lon_f !< Final point longitude (radian)
real(kind_real),intent(in) :: lat_f !< Final point longilatitudetude (radian)
real(kind_real),intent(out) :: dist !< Great-circle distance

! Check that there is no missing value
if (isnotmsr(lon_i).and.isnotmsr(lat_i).and.isnotmsr(lon_f).and.isnotmsr(lat_f)) then
   ! Great-circle distance using Vincenty formula on the unit sphere
   dist = atan2(sqrt((cos(lat_f)*sin(lon_f-lon_i))**2 &
        & +(cos(lat_i)*sin(lat_f)-sin(lat_i)*cos(lat_f)*cos(lon_f-lon_i))**2), &
        & sin(lat_i)*sin(lat_f)+cos(lat_i)*cos(lat_f)*cos(lon_f-lon_i))
else
   call msr(dist)
end if

end subroutine sphere_dist

!----------------------------------------------------------------------
! Subroutine: reduce_arc
!> Purpose: reduce arc to a given distance
!----------------------------------------------------------------------
subroutine reduce_arc(lon_i,lat_i,lon_f,lat_f,maxdist,dist)

implicit none

! Passed variables
real(kind_real),intent(in) :: lon_i    !< Initial point longitude
real(kind_real),intent(in) :: lat_i    !< Initial point latitude
real(kind_real),intent(inout) :: lon_f !< Final point longitude
real(kind_real),intent(inout) :: lat_f !< Final point latitude
real(kind_real),intent(in) :: maxdist  !< Maximum distance
real(kind_real),intent(out) :: dist    !< Effective distance

! Local variable
real(kind_real) :: theta

! Compute distance
call sphere_dist(lon_i,lat_i,lon_f,lat_f,dist)

! Check with the maximum distance
if (dist>maxdist) then
   ! Compute bearing
   theta = atan2(sin(lon_f-lon_i)*cos(lat_f),cos(lat_i)*sin(lat_f)-sin(lat_i)*cos(lat_f)*cos(lon_f-lon_i))

   ! Reduce distance
   dist = maxdist

   ! Compute new point
   lat_f = asin(sin(lat_i)*cos(dist)+cos(lat_i)*sin(dist)*cos(theta))
   lon_f = lon_i+atan2(sin(theta)*sin(dist)*cos(lat_i),cos(dist)-sin(lat_i)*sin(lat_f))
end if

end subroutine reduce_arc

!----------------------------------------------------------------------
! Subroutine: vector_product
!> Purpose: compute normalized vector product
!----------------------------------------------------------------------
subroutine vector_product(v1,v2,vp)

implicit none

! Passed variables
real(kind_real),intent(in) :: v1(3)  !< First vector
real(kind_real),intent(in) :: v2(3)  !< Second vector
real(kind_real),intent(out) :: vp(3) !< Vector product

! Local variable
real(kind_real) :: r

! Vector product
vp(1) = v1(2)*v2(3)-v1(3)*v2(2)
vp(2) = v1(3)*v2(1)-v1(1)*v2(3)
vp(3) = v1(1)*v2(2)-v1(2)*v2(1)

! Normalization
r = sqrt(sum(vp**2))
if (r>0.0) vp = vp/r

end subroutine vector_product

!----------------------------------------------------------------------
! Subroutine: vector_triple_product
!> Purpose: compute vector triple product
!----------------------------------------------------------------------
subroutine vector_triple_product(v1,v2,v3,p)

implicit none

! Passed variables
real(kind_real),intent(in) :: v1(3) !< First vector
real(kind_real),intent(in) :: v2(3) !< Second vector
real(kind_real),intent(in) :: v3(3) !< Third vector
real(kind_real),intent(out) :: p    !< Triple product

! Local variable
real(kind_real) :: vp(3)

! Vector product
vp(1) = v1(2)*v2(3)-v1(3)*v2(2)
vp(2) = v1(3)*v2(1)-v1(1)*v2(3)
vp(3) = v1(1)*v2(2)-v1(2)*v2(1)

! Scalar product
p = sum(vp*v3)

end subroutine vector_triple_product

!----------------------------------------------------------------------
! Subroutine: add
!> Purpose: check if missing and add
!----------------------------------------------------------------------
subroutine add(value,cumul,num,wgt)

implicit none

! Passed variables
real(kind_real),intent(in) :: value        !< Value to add
real(kind_real),intent(inout) :: cumul     !< Cumul
real(kind_real),intent(inout) :: num       !< Number of values
real(kind_real),intent(in),optional :: wgt !< Weight

! Local variables
real(kind_real) :: lwgt

! Initialize weight
lwgt = 1.0
if (present(wgt)) lwgt = wgt

! Add value to cumul
if (isnotmsr(value)) then
   cumul = cumul+lwgt*value
   num = num+1.0
end if

end subroutine add

!----------------------------------------------------------------------
! Subroutine: divide
!> Purpose: check if missing and divide
!----------------------------------------------------------------------
subroutine divide(value,num)

implicit none

! Passed variables
real(kind_real),intent(inout) :: value !< Value to divide
real(kind_real),intent(in) :: num      !< Divider

! Divide cumul by num
if (abs(num)>0.0) then
   value = value/num
else
   call msr(value)
end if

end subroutine divide

!----------------------------------------------------------------------
! Function: fac
!> Purpose: factorial
!----------------------------------------------------------------------
integer function fac(n)

implicit none

! Passed variables
integer,intent(in) :: n !< Argument

! Local variables
integer :: j

if (n<0) call msgerror('factorial requires a non-negative argument')

fac = 1
do j=2,n
   fac = fac*j
end do

end function fac

!----------------------------------------------------------------------
! Subroutine: fit_diag
!> Purpose: diagnostic fit
!----------------------------------------------------------------------
subroutine fit_diag(nc3,nl0r,nl0,l0rl0_to_l0,disth,distvr,rh,rv,fit)

implicit none

! Passed variables
integer,intent(in) :: nc3                        !< Number of classes
integer,intent(in) :: nl0r                       !< Reduced number of levels
integer,intent(in) :: nl0                        !< Number of levels
integer,intent(in) :: l0rl0_to_l0(nl0r,nl0)      !< Reduced level to level
real(kind_real),intent(in) :: disth(nc3)         !< Horizontal distance
real(kind_real),intent(in) :: distvr(nl0r,nl0)   !< Vertical distance
real(kind_real),intent(in) :: rh(nl0)            !< Horizontal support radius
real(kind_real),intent(in) :: rv(nl0)            !< Vertical support radius
real(kind_real),intent(out) :: fit(nc3,nl0r,nl0) !< Fit

! Local variables
integer :: jl0r,jl0,il0,kl0r,kl0,jc3,kc3,ip,jp,np,np_new
integer,allocatable :: plist(:,:),plist_new(:,:)
real(kind_real) :: rhsq,rvsq,distnorm,disttest
real(kind_real),allocatable :: dist(:,:)
logical :: add_to_front

! Initialization
fit = 0.0

!$omp parallel do schedule(static) private(il0,np,jl0r,np_new,ip,jc3,jl0,kc3,kl0r,kl0,rhsq,rvsq,distnorm,disttest,add_to_front), &
!$omp&                             firstprivate(plist,plist_new,dist)
do il0=1,nl0
   ! Allocation
   allocate(plist(nc3*nl0r,2))
   allocate(plist_new(nc3*nl0r,2))
   allocate(dist(nc3,nl0r))

   ! Initialize the front
   np = 1
   call msi(plist)
   plist(1,1) = 1
   do jl0r=1,nl0r
      if (l0rl0_to_l0(jl0r,il0)==il0) plist(1,2) = jl0r
   end do
   dist = 1.0
   dist(plist(1,1),plist(1,2)) = 0.0

   do while (np>0)
      ! Propagate the front
      np_new = 0

      do ip=1,np
         ! Indices of the central point
         jc3 = plist(ip,1)
         jl0r = plist(ip,2)
         jl0 = l0rl0_to_l0(jl0r,il0)

         ! Loop over neighbors
         do kc3=max(jc3-1,1),min(jc3+1,nc3)
            do kl0r=max(jl0r-1,1),min(jl0r+1,nl0r)
               kl0 = l0rl0_to_l0(kl0r,il0)
               if (isnotmsr(rh(jl0)).and.isnotmsr(rh(kl0))) then
                  rhsq = 0.5*(rh(jl0)**2+rh(kl0)**2)
               else
                  rhsq = 0.0
               end if
               if (isnotmsr(rv(jl0)).and.isnotmsr(rv(kl0))) then
                  rvsq = 0.5*(rv(jl0)**2+rv(kl0)**2)
               else
                  rvsq = 0.0
               end if
               distnorm = 0.0
               if (rhsq>0.0) then
                  distnorm = distnorm+(disth(kc3)-disth(jc3))**2/rhsq
               elseif (kc3/=jc3) then
                  distnorm = distnorm+0.5*huge(1.0)
               end if
               if (rvsq>0.0) then
                  distnorm = distnorm+distvr(kl0r,jl0)**2/rvsq
               elseif (kl0r/=jl0r) then
                  distnorm = distnorm+0.5*huge(1.0)
               end if
               disttest = dist(jc3,jl0r)+sqrt(distnorm)
               if (disttest<1.0) then
                  ! Point is inside the support
                  if (disttest<dist(kc3,kl0r)) then
                     ! Update distance
                     dist(kc3,kl0r) = disttest

                     ! Check if the point should be added to the front (avoid duplicates)
                     add_to_front = .true.
                     do jp=1,np_new
                        if ((plist_new(jp,1)==kc3).and.(plist_new(jp,2)==kl0r)) then
                           add_to_front = .false.
                           exit
                        end if
                     end do

                     if (add_to_front) then
                        ! Add point to the front
                        np_new = np_new+1
                        plist_new(np_new,1) = kc3
                        plist_new(np_new,2) = kl0r
                     end if
                  end if
               end if
            end do
         end do
      end do

      ! Copy new front
      np = np_new
      plist(1:np,:) = plist_new(1:np,:)
   end do

   do jl0r=1,nl0r
      do jc3=1,nc3
         ! Gaspari-Cohn (1999) function
         distnorm = dist(jc3,jl0r)
         if (distnorm<1.0) fit(jc3,jl0r,il0) = gc99(distnorm)
      end do
   end do

   ! Release memory
   deallocate(plist)
   deallocate(plist_new)
   deallocate(dist)
end do
!$omp end parallel do

end subroutine fit_diag

!----------------------------------------------------------------------
! Function: gc99
!> Purpose: Gaspari and Cohn (1999) function, with the support radius as a parameter
!----------------------------------------------------------------------
function gc99(distnorm)

! Passed variables
real(kind_real),intent(in) :: distnorm !< Normalized distance

! Returned variable
real(kind_real) :: gc99

! Distance check bound
if (distnorm<0.0) call msgerror('negative normalized distance')

if (.true.) then
   ! Gaspari and Cohn (1999) function
   if (distnorm<0.5) then
      gc99 = 1.0-8.0*distnorm**5+8.0*distnorm**4+5.0*distnorm**3-20.0/3.0*distnorm**2
   else if (distnorm<1.0) then
      gc99 = 8.0/3.0*distnorm**5-8.0*distnorm**4+5.0*distnorm**3+20.0/3.0*distnorm**2-10.0*distnorm+4.0-1.0/(3.0*distnorm)
   else
      gc99 = 0.0
   end if
else
   ! Gaussian equivalent
   if (distnorm<1.0) then
      gc99 = exp(-0.5*(3.53*distnorm)**2)
   else
      gc99 = 0.0
   end if
end if

return

end function gc99

!----------------------------------------------------------------------
! Subroutine: fit_lct
!> Purpose: LCT fit
!----------------------------------------------------------------------
subroutine fit_lct(nc,nl0,dx,dy,dz,dmask,nscales,ncomp,H,coef,fit)

implicit none

! Passed variables
integer,intent(in) :: nc                    !< Number of classes
integer,intent(in) :: nl0                   !< Number of levels
real(kind_real),intent(in) :: dx(nc,nl0)    !< Zonal separation
real(kind_real),intent(in) :: dy(nc,nl0)    !< Meridian separation
real(kind_real),intent(in) :: dz(nl0)       !< Vertical separation
logical,intent(in) :: dmask(nc,nl0)         !< Mask
integer,intent(in) :: nscales               !< Number of LCT scales
integer,intent(in) :: ncomp(nscales)        !< Number of LCT components
real(kind_real),intent(in) :: H(sum(ncomp)) !< LCT components
real(kind_real),intent(in) :: coef(nscales) !< LCT coefficients
real(kind_real),intent(out) :: fit(nc,nl0)  !< Fit

! Local variables
integer :: jl0,jc3,iscales,offset
real(kind_real) :: H11,H22,H33,Hc12,rsq

! Initialization
offset = 0
call msr(fit)

do iscales=1,nscales
   ! Ensure positive-definiteness
   H11 = max(Hmin,H(offset+1))
   H22 = max(Hmin,H(offset+2))
   H33 = max(Hmin,H(offset+3))
   call msr(Hc12)
   if (ncomp(iscales)==4) Hc12 = max(-1.0_kind_real+Hmin,min(H(offset+4),1.0_kind_real-Hmin))

   ! Homogeneous anisotropic approximation
   !$omp parallel do schedule(static) private(jl0,jc3,rsq)
   do jl0=1,nl0
      do jc3=1,nc
         if (dmask(jc3,jl0)) then
            ! Initialization
            if (iscales==1) fit(jc3,jl0) = 0.0

            ! Squared distance
            rsq = H11*dx(jc3,jl0)**2+H22*dy(jc3,jl0)**2+H33*dz(jl0)**2
            if (ncomp(iscales)==4) rsq = rsq+2.0*sqrt(H11*H22)*Hc12*dx(jc3,jl0)*dy(jc3,jl0)

            if (M==0) then
               ! Gaussian function
               fit(jc3,jl0) = fit(jc3,jl0)+coef(iscales)*exp(-0.5*rsq)
            else
               ! Matern function
               fit(jc3,jl0) = fit(jc3,jl0)+coef(iscales)*matern(M,sqrt(rsq))
            end if
         end if
      end do
   end do
   !$omp end parallel do

   ! Update offset
   offset = offset+ncomp(iscales)
end do

end subroutine fit_lct

!----------------------------------------------------------------------
! Function: matern
!> Purpose: compute the normalized diffusion function from eq. (55) of Mirouze and Weaver (2013), for the 3d case (d = 3)
!----------------------------------------------------------------------
real(kind_real) function matern(M,x)

implicit none

! Passed variables
integer,intent(in) :: M         !< Matern function order
real(kind_real),intent(in) :: x !< Argument

! Local variables
integer :: j
real(kind_real) :: xtmp,beta

! Check
if (M<2) call msgerror('M should be larger than 2')
if (mod(M,2)>0) call msgerror('M should be even')

! Initialization
matern = 0.0
beta = 1.0
xtmp = x*sqrt(real(2*M-5,kind_real))

do j=0,M-3
   ! Update sum
   matern = matern+beta*(xtmp)**(M-2-j)

   ! Update beta
   beta = beta*real((j+1+M-2)*(-j+M-2),kind_real)/real(2*(j+1),kind_real)
end do

! Last term and normalization
matern = matern/beta+1.0

! Exponential factor
matern = matern*exp(-xtmp)

end function matern

!----------------------------------------------------------------------
! Function: cholesky
!> Purpose: compute cholesky decomposition
!> Author: Original FORTRAN77 version by Michael Healy, modifications by AJ Miller, FORTRAN90 version by John Burkardt.
!----------------------------------------------------------------------
subroutine cholesky(n,nn,a,u)

implicit none

! Passed variables
integer,intent(in) :: n              !< Matrix rank
integer,intent(in) :: nn             !< Half-matrix size (n*(n-1)/2)
real(kind_real),intent(in) :: a(nn)  !< Matrix
real(kind_real),intent(out) :: u(nn) !< Matrix square-root

! Local variables
integer :: i,icol,ii,irow,j,k,kk,l,m
real(kind_real) :: w,x

! Initialization
ii = 0
j = 1
k = 0
if (nn/=(n*(n+1))/2) then
   call msgerror('wrong size in Cholesky decomposition')
end if

! Factorize column by column, ICOL = column number
do icol=1,n
   ii = ii+icol
   x = eta**2*a(ii)
   l = 0
   kk = 0

   ! IROW = row number within column ICOL
   do irow=1,icol
      kk = kk+irow
      k = k+1
      w = a(k)
      m = j
      do i=1,irow-1
        l = l+1
        w = w-u(l)*u(m)
        m = m+1
      end do
      l = l+1
      if (irow==icol) exit
      if (abs(u(l))>0.0) then
         u(k) = w/u(l)
      else
         if (abs(x*a(k))<w**2) call msgerror('A is not positive semi-definite')
      end if
   end do

   ! End of row, estimate relative accuracy of diagonal element
   if (abs(w)<=abs(eta*a(k))) then
      u(k) = 0.0
   else
      if (w<0.0) call msgerror('A is not positive semi-definite')
   end if
   j = j+icol
end do

end subroutine cholesky

end module tools_func
