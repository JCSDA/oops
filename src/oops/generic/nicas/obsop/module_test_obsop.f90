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
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_geom, only: geomtype,fld_com_gl,fld_com_lg
use type_mpl, only: mpl,mpl_allreduce_sum
use type_odata, only: odatatype,odataloctype,yobs_com_gl,yobs_com_lg
use type_randgen, only: rng,rand_real

implicit none

private
public :: test_adjoint_obsop,test_mpi_obsop,test_mpi_obsop_ad

contains

!----------------------------------------------------------------------
! Subroutine: test_adjoint_obsop
!> Purpose: test observation operator adjoints accuracy
!----------------------------------------------------------------------
subroutine test_adjoint_obsop(odata)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata

! Local variables
real(kind_real) :: sum1,sum2
real(kind_real) :: fld(odata%nc0,odata%nl0),fld_save(odata%nc0,odata%nl0)
real(kind_real) :: yobs(odata%nobs,odata%nl0),yobs_save(odata%nobs,odata%nl0)

! Generate random fields
call rand_real(rng,0.0_kind_real,1.0_kind_real,.true.,fld_save)
call rand_real(rng,0.0_kind_real,1.0_kind_real,.true.,yobs_save)

! Apply direct and adjoint obsservation operators
call apply_obsop(odata,fld_save,yobs)
call apply_obsop_ad(odata,yobs_save,fld)

! Compute adjoint test
sum1 = sum(fld*fld_save)
sum2 = sum(yobs*yobs_save)
write(mpl%unit,'(a7,a,e14.8,a,e14.8,a,e14.8)') '','Observation operator adjoint test: ', &
 & sum1,' / ',sum2,' / ',2.0*abs(sum1-sum2)/abs(sum1+sum2)

end subroutine test_adjoint_obsop

!----------------------------------------------------------------------
! Subroutine: test_mpi_obsop
!> Purpose: test observation operator global/local equivalence
!----------------------------------------------------------------------
subroutine test_mpi_obsop(odata,odataloc)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata
type(odataloctype),intent(inout) :: odataloc

! Local variables
real(kind_real),allocatable :: fld(:,:),fldloc(:,:)
real(kind_real),allocatable :: yobs(:,:),yobsloc(:,:)

! Allocation
if (mpl%main) then
   ! Allocation
   allocate(fld(odata%nc0,odata%nl0))
   allocate(fldloc(odata%nc0,odata%nl0))
   allocate(yobs(odata%nobs,odata%nl0))

   ! Initialization
   call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,fld)
   fldloc = fld
end if
allocate(yobsloc(odataloc%nobsa,odataloc%nl0))

! Global to local
call fld_com_gl(odata%geom,fldloc)

! Global
if (mpl%main) call apply_obsop(odata,fld,yobs)

! Local
call apply_obsop(odataloc,fldloc,yobsloc)

! Local to global
call yobs_com_lg(odata,yobsloc)

! Print difference
if (mpl%main) write(mpl%unit,'(a7,a,e14.8)') '','RMSE between single-proc and multi-procs executions, direct:  ', &
 & sqrt(sum((yobs-yobsloc)**2)/float(odata%nobs*odata%nl0))

end subroutine test_mpi_obsop

!----------------------------------------------------------------------
! Subroutine: test_mpi_obsop_ad
!> Purpose: test adjoint observation operator global/local equivalence
!----------------------------------------------------------------------
subroutine test_mpi_obsop_ad(odata,odataloc)

implicit none

! Passed variables
type(odatatype),intent(inout) :: odata
type(odataloctype),intent(inout) :: odataloc

! Local variables
real(kind_real),allocatable :: fld(:,:),fldloc(:,:)
real(kind_real),allocatable :: yobs(:,:),yobsloc(:,:)

! Allocation
if (mpl%main) then
   ! Allocation
   allocate(yobs(odata%nobs,odata%nl0))
   allocate(yobsloc(odata%nobs,odata%nl0))
   allocate(fld(odata%nc0,odata%nl0))

   ! Initialization
   call rand_real(rng,0.0_kind_real,1.0_kind_real,.false.,yobs)
   yobsloc = yobs
end if
allocate(fldloc(odataloc%nc0a,odataloc%nl0))

! Global to local
call yobs_com_gl(odata,yobsloc)

! Global
if (mpl%main) call apply_obsop_ad(odata,yobs,fld)

! Local
call apply_obsop_ad(odataloc,yobsloc,fldloc)

! Local to global
call fld_com_lg(odata%geom,fldloc)

! Print difference
if (mpl%main) write(mpl%unit,'(a7,a,e14.8)') '','RMSE between single-proc and multi-procs executions, adjoint: ', &
 & sqrt(sum((fld-fldloc)**2)/float(odata%nc0*odata%nl0))

end subroutine test_mpi_obsop_ad

end module module_test_obsop
