!----------------------------------------------------------------------
! Module: obsop_test.f90
!> Purpose: test observation operator parameters
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module obsop_test

use obsop_apply, only: apply_obsop,apply_obsop_ad
use tools_const, only: pi,reqkm,rad2deg,sphere_dist
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr
use type_geom, only: geomtype,fld_com_gl,fld_com_lg
use type_mpl, only: mpl,mpl_dot_prod
use type_odata, only: odatatype,yobs_com_gl,yobs_com_lg
use type_randgen, only: rand_real

implicit none

private
public :: test_adjoint,test_accuracy

contains

!----------------------------------------------------------------------
! Subroutine: test_adjoint
!> Purpose: test observation operator adjoints accuracy
!----------------------------------------------------------------------
subroutine test_adjoint(odata)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata !< Observation operator data

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real) :: fld(odata%geom%nc0a,odata%geom%nl0),fld_save(odata%geom%nc0a,odata%geom%nl0)
real(kind_real) :: yobs(odata%nobsa,odata%geom%nl0),yobs_save(odata%nobsa,odata%geom%nl0)

! Associate
associate(geom=>odata%geom)

! Generate random fields
call rand_real(0.0_kind_real,1.0_kind_real,fld_save)
call rand_real(0.0_kind_real,1.0_kind_real,yobs_save)

! Apply direct and adjoint obsservation operators
call apply_obsop(geom,odata,fld_save,yobs)
call apply_obsop_ad(geom,odata,yobs_save,fld)

! Compute adjoint test
call mpl_dot_prod(fld,fld_save,sum1)
call mpl_dot_prod(yobs,yobs_save,sum2)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Observation operator adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

! End associate
end associate

end subroutine test_adjoint

!----------------------------------------------------------------------
! Subroutine: test_accuracy
!> Purpose: test observation operator accuracy
!----------------------------------------------------------------------
subroutine test_accuracy(odata)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata !< Observation operator data

! Local variables
integer :: il0,iobs,ind(1)
real(kind_real) :: dist(odata%nobs)
real(kind_real),allocatable :: lon(:,:),lat(:,:)
real(kind_real),allocatable :: ylon(:,:),ylat(:,:)

! Associate
associate(geom=>odata%geom)

! Allocation
allocate(ylon(odata%nobs,geom%nl0))
allocate(ylat(odata%nobs,geom%nl0))
allocate(lon(geom%nc0,geom%nl0))
allocate(lat(geom%nc0,geom%nl0))

if (mpl%main) then
   ! Initialization
   do il0=1,geom%nl0
      lon(:,il0) = geom%lon
      lat(:,il0) = geom%lat
   end do
end if

! Global to local
call fld_com_gl(geom,lon)
call fld_com_gl(geom,lat)

! Apply obsop
call apply_obsop(geom,odata,lon,ylon)
call apply_obsop(geom,odata,lat,ylat)

! Local to global
call yobs_com_lg(odata,ylon)
call yobs_com_lg(odata,ylat)

! Print difference
if (mpl%main) then
   ! Initialization
   call msr(dist)

   do iobs=1,odata%nobs
      ! Remove points close to the longitude discontinuity and to the poles
      if ((abs(odata%lonobs(iobs))<0.8*pi).and.(abs(odata%latobs(iobs))<0.4*pi)) then
         call sphere_dist(ylon(iobs,1),ylat(iobs,1),odata%lonobs(iobs),odata%latobs(iobs),dist(iobs))
         dist(iobs) = dist(iobs)*reqkm
      end if
   end do
   if (count(isnotmsr(dist))>0) then
      write(mpl%unit,'(a7,a,f8.2,a,f8.2,a,f8.2,a)') '','Interpolation error (min/mean/max): ',minval(dist,mask=isnotmsr(dist)), &
    & ' km / ',sum(dist,mask=isnotmsr(dist))/float(count(isnotmsr(dist))),' km / ',maxval(dist,mask=isnotmsr(dist)),' km'
      ind = maxloc(dist,mask=isnotmsr(dist))
      write(mpl%unit,'(a7,a)') '','Max. interpolation error location (lon/lat): '
      write(mpl%unit,'(a10,a14,f8.2,a,f8.2,a)') '','Observation:  ',odata%lonobs(ind(1))*rad2deg, &
    & ' deg. / ' ,odata%latobs(ind(1))*rad2deg,' deg.'
      write(mpl%unit,'(a10,a14,f8.2,a,f8.2,a)') '','Interpolation:',ylon(ind(1),1)*rad2deg, &
    & ' deg. / ' ,ylat(ind(1),1)*rad2deg,' deg.'
   else
      call msgerror('observations are out of the test windows')
   end if
end if

! End associate
end associate

end subroutine test_accuracy

end module obsop_test
