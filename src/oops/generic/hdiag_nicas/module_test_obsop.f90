!----------------------------------------------------------------------
! Module: module_test_obsop.f90
!> Purpose: test observation operator parameters
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_test_obsop

use module_apply_obsop, only: apply_obsop,apply_obsop_ad
use tools_const, only: pi,reqkm,rad2deg,sphere_dist
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr
use type_geom, only: geomtype,fld_com_gl,fld_com_lg
use type_mpl, only: mpl,mpl_allreduce_sum
use type_odata, only: odatatype,odataloctype,yobs_com_gl,yobs_com_lg
use type_randgen, only: rand_real

implicit none

private
public :: test_adjoint_obsop,test_mpi_obsop,test_mpi_obsop_ad,test_obsop

contains

!----------------------------------------------------------------------
! Subroutine: test_adjoint_obsop
!> Purpose: test observation operator adjoints accuracy
!----------------------------------------------------------------------
subroutine test_adjoint_obsop(odata)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata !< Observation operator data

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real) :: fld(odata%geom%nc0,odata%geom%nl0),fld_save(odata%geom%nc0,odata%geom%nl0)
real(kind_real) :: yobs(odata%nobs,odata%geom%nl0),yobs_save(odata%nobs,odata%geom%nl0)

! Generate random fields
call rand_real(0.0_kind_real,1.0_kind_real,fld_save)
call rand_real(0.0_kind_real,1.0_kind_real,yobs_save)

! Apply direct and adjoint obsservation operators
call apply_obsop(odata,fld_save,yobs)
call apply_obsop_ad(odata,yobs_save,fld)

! Compute adjoint test
sum1 = sum(fld*fld_save)
sum2 = sum(yobs*yobs_save)
write(mpl%unit,'(a7,a,e15.8,a,e15.8,a,e15.8)') '','Observation operator adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

end subroutine test_adjoint_obsop

!----------------------------------------------------------------------
! Subroutine: test_mpi_obsop
!> Purpose: test observation operator global/local equivalence
!----------------------------------------------------------------------
subroutine test_mpi_obsop(odata,odataloc)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata       !< Observation operator data
type(odataloctype),intent(inout) :: odataloc !< Observation operator data, local

! Local variables
real(kind_real),allocatable :: fld(:,:),fldloc(:,:)
real(kind_real),allocatable :: yobs(:,:),yobsloc(:,:)

! Associate
associate(geom=>odata%geom)

! Allocation
if (mpl%main) then
   ! Allocation
   allocate(fld(geom%nc0,geom%nl0))
   allocate(fldloc(geom%nc0,geom%nl0))
   allocate(yobs(odata%nobs,geom%nl0))

   ! Initialization
   call rand_real(0.0_kind_real,1.0_kind_real,fld)
   fldloc = fld
end if
allocate(yobsloc(odataloc%nobsa,geom%nl0))

! Global to local
call fld_com_gl(geom,fldloc)

! Global
if (mpl%main) call apply_obsop(odata,fld,yobs)

! Local
call apply_obsop(odataloc,fldloc,yobsloc)

! Local to global
call yobs_com_lg(odata,yobsloc)

! Print difference
if (mpl%main) write(mpl%unit,'(a7,a,e15.8)') '','RMSE between single-proc and multi-procs executions, direct:  ', &
 & sqrt(sum((yobs-yobsloc)**2)/float(odata%nobs*geom%nl0))

! End associate
end associate

end subroutine test_mpi_obsop

!----------------------------------------------------------------------
! Subroutine: test_mpi_obsop_ad
!> Purpose: test adjoint observation operator global/local equivalence
!----------------------------------------------------------------------
subroutine test_mpi_obsop_ad(odata,odataloc)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata       !< Observation operator data
type(odataloctype),intent(inout) :: odataloc !< Observation operator data, local

! Local variables
real(kind_real),allocatable :: fld(:,:),fldloc(:,:)
real(kind_real),allocatable :: yobs(:,:),yobsloc(:,:)

! Associate
associate(geom=>odata%geom)

! Allocation
if (mpl%main) then
   ! Allocation
   allocate(yobs(odata%nobs,geom%nl0))
   allocate(yobsloc(odata%nobs,geom%nl0))
   allocate(fld(geom%nc0,geom%nl0))

   ! Initialization
   call rand_real(0.0_kind_real,1.0_kind_real,yobs)
   yobsloc = yobs
end if
allocate(fldloc(odataloc%nc0a,geom%nl0))

! Global to local
call yobs_com_gl(odata,yobsloc)

! Global
if (mpl%main) call apply_obsop_ad(odata,yobs,fld)

! Local
call apply_obsop_ad(odataloc,yobsloc,fldloc)

! Local to global
call fld_com_lg(geom,fldloc)

! Print difference
if (mpl%main) write(mpl%unit,'(a7,a,e15.8)') '','RMSE between single-proc and multi-procs executions, adjoint: ', &
 & sqrt(sum((fld-fldloc)**2)/float(geom%nc0*geom%nl0))

! End associate
end associate

end subroutine test_mpi_obsop_ad

!----------------------------------------------------------------------
! Subroutine: test_obsop
!> Purpose: test observation operator precision
!----------------------------------------------------------------------
subroutine test_obsop(odata,odataloc)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata       !< Observation operator data
type(odataloctype),intent(inout) :: odataloc !< Observation operator data, local

! Local variables
integer :: il0,iobs,ind(1)
real(kind_real) :: dist(odata%nobs)
real(kind_real),allocatable :: lon(:,:),lat(:,:)
real(kind_real),allocatable :: ylon(:,:),ylat(:,:)

! Associate
associate(geom=>odata%geom)

! Allocation
if (mpl%main) then
   ! Allocation
   allocate(ylon(odata%nobs,geom%nl0))
   allocate(ylat(odata%nobs,geom%nl0))
   allocate(lon(geom%nc0,geom%nl0))
   allocate(lat(geom%nc0,geom%nl0))

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
call apply_obsop(odataloc,lon,ylon)
call apply_obsop(odataloc,lat,ylat)

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
   write(mpl%unit,'(a7,a,f8.2,a,f8.2,a,f8.2,a)') '','Interpolation error (min/mean/max): ',minval(dist,mask=isnotmsr(dist)), &
 & ' km / ',sum(dist,mask=isnotmsr(dist))/float(count(isnotmsr(dist))),' km / ',maxval(dist,mask=isnotmsr(dist)),' km'
   ind = maxloc(dist,mask=isnotmsr(dist))
   write(mpl%unit,'(a7,a)') '','Max. interpolation error location: '
   write(mpl%unit,'(a10,a32,f8.2,a,f8.2,a)') '','observation location (lon/lat):',odata%lonobs(ind(1))*rad2deg, &
 & ' deg. / ' ,odata%latobs(ind(1))*rad2deg,' deg.'
   write(mpl%unit,'(a10,a32,f8.2,a,f8.2,a)') '','interpolation location (lon/lat):',ylon(ind(1),1)*rad2deg, &
 & ' deg. / ' ,ylat(ind(1),1)*rad2deg,' deg.'
end if

! End associate
end associate

end subroutine test_obsop

end module module_test_obsop
