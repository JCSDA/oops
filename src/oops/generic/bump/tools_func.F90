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

use tools_asa007, only: asa007_cholesky,asa007_syminv
use tools_const, only: pi
use tools_kinds, only: kind_real
use tools_missing, only: msi,msr,isnotmsr
use type_mpl, only: mpl_type

implicit none

real(kind_real),parameter :: gc_gau = 0.28            !< Gaussian to GC99 factor (empirical)
real(kind_real),parameter :: Dmin = 1.0e-12_kind_real !< Minimum tensor diagonal value
integer,parameter :: M = 0                            !< Number of implicit itteration for the Matern function (GC 99 function if M = -1 and Gaussian function if M = 0)

private
public :: lonlatmod,sphere_dist,reduce_arc,vector_product,vector_triple_product,add,divide, &
        & fit_diag,fit_diag_dble,gc99,fit_lct,cholesky,syminv

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
! Subroutine: sphere_dist
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
! Subroutine: fit_diag
!> Purpose: diagnostic fit
!----------------------------------------------------------------------
subroutine fit_diag(mpl,nc3,nl0r,nl0,l0rl0_to_l0,disth,distv,rh,rv,fit)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl                 !< MPI data
integer,intent(in) :: nc3                        !< Number of classes
integer,intent(in) :: nl0r                       !< Reduced number of levels
integer,intent(in) :: nl0                        !< Number of levels
integer,intent(in) :: l0rl0_to_l0(nl0r,nl0)      !< Reduced level to level
real(kind_real),intent(in) :: disth(nc3)         !< Horizontal distance
real(kind_real),intent(in) :: distv(nl0,nl0)     !< Vertical distance
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
                  distnorm = distnorm+distv(kl0,jl0)**2/rvsq
               elseif (kl0/=jl0) then
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
         fit(jc3,jl0r,il0) = gc99(mpl,distnorm)
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
! Subroutine: fit_diag_dble
!> Purpose: diagnostic fit, double-fit
!----------------------------------------------------------------------
subroutine fit_diag_dble(mpl,nc3,nl0r,nl0,l0rl0_to_l0,disth,distv,rh,rv,rv_rfac,rv_coef,fit)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl                 !< MPI data
integer,intent(in) :: nc3                        !< Number of classes
integer,intent(in) :: nl0r                       !< Reduced number of levels
integer,intent(in) :: nl0                        !< Number of levels
integer,intent(in) :: l0rl0_to_l0(nl0r,nl0)      !< Reduced level to level
real(kind_real),intent(in) :: disth(nc3)         !< Horizontal distance
real(kind_real),intent(in) :: distv(nl0,nl0)     !< Vertical distance
real(kind_real),intent(in) :: rh(nl0)            !< Horizontal support radius
real(kind_real),intent(in) :: rv(nl0)            !< Vertical support radius
real(kind_real),intent(in) :: rv_rfac(nl0)       !< Vertical fit support radius ratio for the positive component
real(kind_real),intent(in) :: rv_coef(nl0)       !< Vertical fit coefficient
real(kind_real),intent(out) :: fit(nc3,nl0r,nl0) !< Fit

! Local variables
integer :: jl0r,jl0,il0,kl0r,kl0,jc3,kc3,ip,jp,np,np_new
integer,allocatable :: plist(:,:),plist_new(:,:)
real(kind_real) :: rhsq,rvsq,distnorm,disttest,rfac,coef,distnormv,distnormh
real(kind_real),allocatable :: dist(:,:)
logical :: add_to_front

! Initialization
fit = 0.0

!$omp parallel do schedule(static) private(il0,np,jl0r,np_new,ip,jc3,jl0,kc3,kl0r,kl0,rhsq,rvsq,distnorm,disttest,add_to_front), &
!$omp&                             private(coef,distnormv,distnormh) firstprivate(plist,plist_new,dist)
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
                  distnorm = distnorm+distv(kl0,jl0)**2/rvsq
               elseif (kl0/=jl0) then
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
      jl0 = l0rl0_to_l0(jl0r,il0)
      distnormv = dist(1,jl0r)
      if ((abs(rv_rfac(il0))<0.0).or.(abs(rv_rfac(jl0))<0.0)) then
          rfac = 0.0
      else
          rfac = sqrt(rv_rfac(il0)*rv_rfac(jl0))
      end if
      if ((abs(rv_coef(il0))<0.0).or.(abs(rv_coef(jl0))<0.0)) then
          coef = 0.0
      else
          coef = sqrt(rv_coef(il0)*rv_coef(jl0))
      end if
      do jc3=1,nc3
         ! Double Gaspari-Cohn (1999) function
         distnorm = dist(jc3,jl0r)
         distnormh = sqrt(distnorm**2-distnormv**2)
         fit(jc3,jl0r,il0) = gc99(mpl,distnormh)*((1.0+coef)*gc99(mpl,distnormv/rfac)-coef*gc99(mpl,distnormv))
      end do
   end do

   ! Release memory
   deallocate(plist)
   deallocate(plist_new)
   deallocate(dist)
end do
!$omp end parallel do

end subroutine fit_diag_dble

!----------------------------------------------------------------------
! Function: gc99
!> Purpose: Gaspari and Cohn (1999) function, with the support radius as a parameter
!----------------------------------------------------------------------
function gc99(mpl,distnorm)

! Passed variables
type(mpl_type),intent(in) :: mpl       !< MPI data
real(kind_real),intent(in) :: distnorm !< Normalized distance

! Returned variable
real(kind_real) :: gc99

! Distance check bound
if (distnorm<0.0) call mpl%abort('negative normalized distance')

! Gaspari and Cohn (1999) function
if (distnorm<0.5) then
   gc99 = -8.0*distnorm**5+8.0*distnorm**4+5.0*distnorm**3-20.0/3.0*distnorm**2+1.0
else if (distnorm<1.0) then
   gc99 = 8.0/3.0*distnorm**5-8.0*distnorm**4+5.0*distnorm**3+20.0/3.0*distnorm**2-10.0*distnorm+4.0-1.0/(3.0*distnorm)
else
   gc99 = 0.0
end if

! Enforce positivity
gc99 = max(gc99,0.0)

end function gc99

!----------------------------------------------------------------------
! Subroutine: fit_lct
!> Purpose: LCT fit
!----------------------------------------------------------------------
subroutine fit_lct(mpl,nc,nl0,dx,dy,dz,dmask,nscales,ncomp,D,coef,fit)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl            !< MPI data
integer,intent(in) :: nc                    !< Number of classes
integer,intent(in) :: nl0                   !< Number of levels
real(kind_real),intent(in) :: dx(nc,nl0)    !< Zonal separation
real(kind_real),intent(in) :: dy(nc,nl0)    !< Meridian separation
real(kind_real),intent(in) :: dz(nl0)       !< Vertical separation
logical,intent(in) :: dmask(nc,nl0)         !< Mask
integer,intent(in) :: nscales               !< Number of LCT scales
integer,intent(in) :: ncomp(nscales)        !< Number of LCT components
real(kind_real),intent(in) :: D(sum(ncomp)) !< LCT components
real(kind_real),intent(in) :: coef(nscales) !< LCT coefficients
real(kind_real),intent(out) :: fit(nc,nl0)  !< Fit

! Local variables
integer :: jl0,jc3,iscales,offset
real(kind_real) :: Hcoef(nscales),D11,D22,D33,D12,H11,H22,H33,H12,rsq,det,distnorm

! Initialization
offset = 0
call msr(fit)

! Coefficients
Hcoef = max(Dmin,min(coef,1.0_kind_real))
Hcoef = Hcoef/sum(Hcoef)

do iscales=1,nscales
   ! Ensure positive-definiteness of D
   D11 = max(Dmin,D(offset+1))
   D22 = max(Dmin,D(offset+2))
   call msr(D33)
   if (nl0>1) D33 = max(Dmin,D(offset+3))
   call msr(D12)
   if (ncomp(iscales)==4) D12 = sqrt(D11*D22)*max(-1.0_kind_real+Dmin,min(D(offset+4),1.0_kind_real-Dmin))

   ! Inverse D to get H
   if (ncomp(iscales)==3) then
      det = D11*D22
   else
      det = D11*D22-D12**2
   end if
   H11 = D22/det
   H22 = D11/det
   call msr(H33)
   if (nl0>1) H33 = 1.0/D33
   call msr(H12)
   if (ncomp(iscales)==4) H12 = -D12/det

   ! Homogeneous anisotropic approximation
   !$omp parallel do schedule(static) private(jl0,jc3,rsq)
   do jl0=1,nl0
      do jc3=1,nc
         if (dmask(jc3,jl0)) then
            ! Initialization
            if (iscales==1) fit(jc3,jl0) = 0.0

            ! Squared distance
            rsq = H11*dx(jc3,jl0)**2+H22*dy(jc3,jl0)**2
            if (nl0>1) rsq = rsq+H33*dz(jl0)**2
            if (ncomp(iscales)==4) rsq = rsq+2.0*H12*dx(jc3,jl0)*dy(jc3,jl0)

            if (M==-1) then
               ! Gaspari-Cohn 1999 function
               distnorm = sqrt(rsq)*gc_gau
               fit(jc3,jl0) = fit(jc3,jl0)+Hcoef(iscales)*gc99(mpl,distnorm)
            elseif (M==0) then
               ! Gaussian function
               if (rsq<40.0) fit(jc3,jl0) = fit(jc3,jl0)+Hcoef(iscales)*exp(-0.5*rsq)
            else
               ! Matern function
               fit(jc3,jl0) = fit(jc3,jl0)+Hcoef(iscales)*matern(mpl,M,sqrt(rsq))
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
real(kind_real) function matern(mpl,M,x)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl !< MPI data
integer,intent(in) :: M          !< Matern function order
real(kind_real),intent(in) :: x  !< Argument

! Local variables
integer :: j
real(kind_real) :: xtmp,beta

! Check
if (M<2) call mpl%abort('M should be larger than 2')
if (mod(M,2)>0) call mpl%abort('M should be even')

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
! Subroutine: cholesky
!> Purpose: compute cholesky decomposition
!> Author: Original FORTRAN77 version by Michael Healy, modifications by AJ Miller, FORTRAN90 version by John Burkardt.
!----------------------------------------------------------------------
subroutine cholesky(mpl,n,a,u)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl      !< MPI data
integer,intent(in) :: n               !< Matrix rank
real(kind_real),intent(in) :: a(n,n)  !< Matrix
real(kind_real),intent(out) :: u(n,n) !< Matrix square-root

! Local variables
integer :: nn,i,j,ij
real(kind_real),allocatable :: apack(:),upack(:)

! Allocation
nn = (n*(n+1))/2
allocate(apack(nn))
allocate(upack(nn))

! Pack matrix
ij = 0
do i=1,n
   do j=1,i
      ij = ij+1
      apack(ij) = a(i,j)
   end do
end do

! Cholesky decomposition
call asa007_cholesky(mpl,n,nn,apack,upack)

! Unpack matrix
ij = 0
u = 0.0
do i=1,n
   do j=1,i
      ij = ij+1
      u(i,j) = upack(ij)
   end do
end do

end subroutine cholesky

!----------------------------------------------------------------------
! Subroutine: syminv
!> Purpose: compute inverse of a symmetric matrix
!> Author: Original FORTRAN77 version by Michael Healy, modifications by AJ Miller, FORTRAN90 version by John Burkardt.
!----------------------------------------------------------------------
subroutine syminv(mpl,n,a,c)

implicit none

! Passed variables
type(mpl_type),intent(in) :: mpl      !< MPI data
integer,intent(in) :: n               !< Matrix rank
real(kind_real),intent(in) :: a(n,n)  !< Matrix
real(kind_real),intent(out) :: c(n,n) !< Matrix inverse

! Local variables
integer :: nn,i,j,ij
real(kind_real),allocatable :: apack(:),cpack(:)

! Allocation
nn = (n*(n+1))/2
allocate(apack(nn))
allocate(cpack(nn))

! Pack matrix
ij = 0
do i=1,n
   do j=1,i
      ij = ij+1
      apack(ij) = a(i,j)
   end do
end do

! Matrix inversion
call asa007_syminv(mpl,n,nn,apack,cpack)

! Unpack matrix
ij = 0
do i=1,n
   do j=1,i
      ij = ij+1
      c(i,j) = cpack(ij)
      c(j,i) = c(i,j)
   end do
end do

end subroutine syminv

end module tools_func
