!----------------------------------------------------------------------
! Module: module_apply_com.f90
!> Purpose: communication routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_com

use tools_display, only: msgerror
use tools_kinds, only: kind_real
use type_fields, only: fldtype,alphatype
use type_mpl, only: mpl,mpl_send,mpl_recv,mpl_alltoallv
use type_ndata, only: ndatatype,ndataloctype

implicit none

private
public :: fld_com_gl,fld_com_lg
public :: alpha_com_AB,alpha_com_BA,alpha_com_CA,alpha_copy_AB,alpha_copy_AC,alpha_copy_BC

contains

!----------------------------------------------------------------------
! Subroutine: fld_com_gl
!> Purpose: communicate full field from global to local distribution
!----------------------------------------------------------------------
subroutine fld_com_gl(ndata,ndataloc,fld)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata       !< Sampling data
type(ndataloctype),intent(in) :: ndataloc !< Sampling data, local
type(fldtype),intent(inout) :: fld        !< Field

! Local variables
integer :: ic0,ic0a,iproc
real(kind_real),allocatable :: sbuf(:,:),rbuf(:)

! Allocation
if (mpl%main) allocate(sbuf(ndata%nc0amax*ndata%nl0,mpl%nproc))
allocate(rbuf(ndataloc%nc0a*ndata%nl0))
allocate(fld%vala(ndataloc%nc0a,ndataloc%nl0))

if (mpl%main) then
   ! Sending buffer
   do ic0=1,ndata%nc0
      iproc = ndata%ic0_to_iproc(ic0)
      ic0a = ndata%ic0_to_ic0a(ic0)
      sbuf((ic0a-1)*ndata%nl0+1:ic0a*ndata%nl0,iproc) = fld%val(ic0,:)
   end do
end if

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf = sbuf(1:ndataloc%nc0a*ndata%nl0,iproc)
      else
         ! Send data to iproc
         call mpl_send(ndataloc%nc0a*ndata%nl0,sbuf(1:ndataloc%nc0a*ndata%nl0,iproc),iproc,mpl%tag)
      end if
   end do
else
   ! Receive data from ioproc
   call mpl_recv(ndataloc%nc0a*ndata%nl0,rbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

! Received buffer
do ic0a=1,ndataloc%nc0a
   fld%vala(ic0a,:) = rbuf((ic0a-1)*ndata%nl0+1:ic0a*ndata%nl0)
end do

! Release memory
if (mpl%main) deallocate(sbuf)
deallocate(rbuf)
if (mpl%main) deallocate(fld%val)

end subroutine fld_com_gl

!----------------------------------------------------------------------
! Subroutine: fld_com_lg
!> Purpose: communicate full field from local to global distribution
!----------------------------------------------------------------------
subroutine fld_com_lg(ndata,ndataloc,fld)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata       !< Sampling data
type(ndataloctype),intent(in) :: ndataloc !< Sampling data, local
type(fldtype),intent(inout) :: fld        !< Field

! Local variables
integer :: ic0,ic0a,iproc
real(kind_real),allocatable :: sbuf(:),rbuf(:,:)

! Allocation
allocate(sbuf(ndataloc%nc0a*ndata%nl0))
if (mpl%main) allocate(rbuf(ndata%nc0amax*ndata%nl0,mpl%nproc))
if (mpl%main) allocate(fld%val(ndata%nc0,ndata%nl0))

! Sending buffer
do ic0a=1,ndataloc%nc0a
   sbuf((ic0a-1)*ndata%nl0+1:ic0a*ndata%nl0) = fld%vala(ic0a,:)
end do

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc==mpl%ioproc) then
         ! Copy data
         rbuf(1:ndataloc%nc0a*ndata%nl0,iproc) = sbuf
      else
         ! Receive data from iproc
         call mpl_recv(ndataloc%nc0a*ndata%nl0,rbuf(1:ndataloc%nc0a*ndata%nl0,iproc),iproc,mpl%tag)
      end if
   end do
else
   ! Sending data to iproc
   call mpl_send(ndataloc%nc0a*ndata%nl0,sbuf,mpl%ioproc,mpl%tag)
end if
mpl%tag = mpl%tag+1

if (mpl%main) then
   ! Received buffer
   do ic0=1,ndata%nc0
      iproc = ndata%ic0_to_iproc(ic0)
      ic0a = ndata%ic0_to_ic0a(ic0)
      fld%val(ic0,:) = rbuf((ic0a-1)*ndata%nl0+1:ic0a*ndata%nl0,iproc)
   end do
end if

! Release memory
deallocate(sbuf)
if (mpl%main) deallocate(rbuf)
deallocate(fld%vala)

end subroutine fld_com_lg

!----------------------------------------------------------------------
! Subroutine: alpha_com_AB
!> Purpose: communicate reduced field from zone A to zone B
!----------------------------------------------------------------------
subroutine alpha_com_AB(ndataloc,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Local variables
real(kind_real) :: sbuf(ndataloc%AB%nexcl),rbuf(ndataloc%AB%nhalo)

! Allocation
allocate(alpha%valb(ndataloc%nsb))

! Prepare buffers to send
sbuf = alpha%vala(ndataloc%AB%excl)

! Communication
call mpl_alltoallv(ndataloc%AB%nexcl,sbuf,ndataloc%AB%jexclcounts,ndataloc%AB%jexcldispl, &
 & ndataloc%AB%nhalo,rbuf,ndataloc%AB%jhalocounts,ndataloc%AB%jhalodispl)

! Copy zone A into zone B
alpha%valb(ndataloc%isa_to_isb) = alpha%vala

! Copy halo into zone B
alpha%valb(ndataloc%AB%halo) = rbuf

! Release memory
deallocate(alpha%vala)

end subroutine alpha_com_AB

!----------------------------------------------------------------------
! Subroutine: alpha_com_BA
!> Purpose: communicate reduced field from zone B to zone A
!----------------------------------------------------------------------
subroutine alpha_com_BA(ndataloc,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Local variables
real(kind_real) :: sbuf(ndataloc%AB%nhalo),rbuf(ndataloc%AB%nexcl)

! Allocation
allocate(alpha%vala(ndataloc%nsa))

! Prepare buffers to send
sbuf = alpha%valb(ndataloc%AB%halo)

! Communication
call mpl_alltoallv(ndataloc%AB%nhalo,sbuf,ndataloc%AB%jhalocounts,ndataloc%AB%jhalodispl, &
 & ndataloc%AB%nexcl,rbuf,ndataloc%AB%jexclcounts,ndataloc%AB%jexcldispl)

! Copy zone B into zone A
alpha%vala = alpha%valb(ndataloc%isa_to_isb)

! Copy halo into zone B
alpha%vala(ndataloc%AB%excl) = alpha%vala(ndataloc%AB%excl)+rbuf

! Release memory
deallocate(alpha%valb)

end subroutine alpha_com_BA

!----------------------------------------------------------------------
! Subroutine: alpha_com_CA
!> Purpose: communicate reduced field from zone C to zone A
!----------------------------------------------------------------------
subroutine alpha_com_CA(ndataloc,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Local variables
real(kind_real) :: sbuf(ndataloc%AC%nhalo),rbuf(ndataloc%AC%nexcl)

! Allocation
allocate(alpha%vala(ndataloc%nsa))

! Prepare buffers to send
sbuf = alpha%valc(ndataloc%AC%halo)

! Communication
call mpl_alltoallv(ndataloc%AC%nhalo,sbuf,ndataloc%AC%jhalocounts,ndataloc%AC%jhalodispl, &
 & ndataloc%AC%nexcl,rbuf,ndataloc%AC%jexclcounts,ndataloc%AC%jexcldispl)

! Copy zone C into zone A
alpha%vala = alpha%valc(ndataloc%isa_to_isc)

! Copy halo into zone A
alpha%vala(ndataloc%AC%excl) = alpha%vala(ndataloc%AC%excl)+rbuf

! Release memory
deallocate(alpha%valc)

end subroutine alpha_com_CA

!----------------------------------------------------------------------
! Subroutine: alpha_copy_AB
!> Purpose: copy reduced field from zone A to zone B
!----------------------------------------------------------------------
subroutine alpha_copy_AB(ndataloc,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Allocation
allocate(alpha%valb(ndataloc%nsb))

! Initialize
alpha%valb = 0.0

! Copy zone A into zone B
alpha%valb(ndataloc%isa_to_isb) = alpha%vala

! Release memory
deallocate(alpha%vala)

end subroutine alpha_copy_AB

!----------------------------------------------------------------------
! Subroutine: alpha_copy_AC
!> Purpose: copy reduced field from zone A to zone C
!----------------------------------------------------------------------
subroutine alpha_copy_AC(ndataloc,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Allocation
allocate(alpha%valc(ndataloc%nsc))

! Initialize
alpha%valc = 0.0

! Copy zone A into zone C
alpha%valc(ndataloc%isa_to_isc) = alpha%vala

! Release memory
deallocate(alpha%vala)

end subroutine alpha_copy_AC

!----------------------------------------------------------------------
! Subroutine: alpha_copy_BC
!> Purpose: copy reduced field from zone B to zone C
!----------------------------------------------------------------------
subroutine alpha_copy_BC(ndataloc,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Allocation
allocate(alpha%valc(ndataloc%nsc))

! Initialization
alpha%valc = 0.0

! Copy zone B into zone C
alpha%valc(ndataloc%isb_to_isc) = alpha%valb

! Release memory
deallocate(alpha%valb)

end subroutine alpha_copy_BC

end module module_apply_com
