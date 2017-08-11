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
use type_fields, only: fldtype,alphatype,buftype
use type_sdata, only: sdatatype,sdatampitype
implicit none

private
public :: fld_com_gl,fld_com_lg
public :: alpha_com_AB,alpha_com_BA,alpha_com_CA,alpha_copy_AB,alpha_copy_AC,alpha_copy_BC

contains

!----------------------------------------------------------------------
! Subroutine: fld_com_gl
!> Purpose: communicate full field from global to local distribution
!----------------------------------------------------------------------
subroutine fld_com_gl(sdata,fld_global,fld_local)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                   !< Sampling data
type(fldtype),intent(inout) :: fld_global             !< Global field
type(fldtype),intent(inout) :: fld_local(sdata%nproc) !< Local field

! Local variables
integer :: ic0,iproc,ic0a

! Allocation
do iproc=1,sdata%nproc
   if (.not.allocated(fld_local(iproc)%vala)) allocate(fld_local(iproc)%vala(sdata%mpi(iproc)%nc0a,sdata%mpi(iproc)%nl0))
end do

! Loop over cells
do ic0=1,sdata%nc0
   iproc = sdata%ic0_to_iproc(ic0)
   ic0a = sdata%ic0_to_ic0a(ic0)
   fld_local(iproc)%vala(ic0a,:) = fld_global%val(ic0,:)
end do

end subroutine fld_com_gl

!----------------------------------------------------------------------
! Subroutine: fld_com_lg
!> Purpose: communicate full field from local to global distribution
!----------------------------------------------------------------------
subroutine fld_com_lg(sdata,fld_local,fld_global)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                   !< Sampling data
type(fldtype),intent(inout) :: fld_local(sdata%nproc) !< Local field
type(fldtype),intent(inout) :: fld_global             !< Global field

! Local variables
integer :: ic0,iproc,ic0a

! Allocation
if (.not.allocated(fld_global%val)) allocate(fld_global%val(sdata%nc0,sdata%nl0))

! Loop over cells
do ic0=1,sdata%nc0
   iproc = sdata%ic0_to_iproc(ic0)
   ic0a = sdata%ic0_to_ic0a(ic0)
   fld_global%val(ic0,:) = fld_local(iproc)%vala(ic0a,:)
end do

end subroutine fld_com_lg

!----------------------------------------------------------------------
! Subroutine: alpha_com_AB
!> Purpose: communicate reduced field from zone A to zone B
!----------------------------------------------------------------------
subroutine alpha_com_AB(sdata,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                 !< Sampling data
type(alphatype),intent(inout) :: alpha(sdata%nproc) !< Subgrid variable

! Local variables
integer :: iproc,jproc,ir_s,ir_e,is_s,is_e
type(buftype) :: bufs(sdata%nproc),bufr(sdata%nproc)

! Allocation
do iproc=1,sdata%nproc
   allocate(alpha(iproc)%valb(sdata%mpi(iproc)%nsb))
end do

if (sdata%nproc==1) then
   ! Copy data
   alpha(1)%valb = alpha(1)%vala
elseif (sdata%nproc>1) then
   ! Allocation
   do iproc=1,sdata%nproc
      allocate(bufs(iproc)%val(sdata%mpi(iproc)%AB%nexcl))
      allocate(bufr(iproc)%val(sdata%mpi(iproc)%AB%nhalo))
   end do

   ! Prepare buffers to send
   do iproc=1,sdata%nproc
      bufs(iproc)%val = alpha(iproc)%vala(sdata%mpi(iproc)%AB%excl)
   end do

   ! Emulate mpi_alltoallv
   do iproc=1,sdata%nproc
      do jproc=1,sdata%nproc
         ! Communication iproc->jproc
         is_s = sdata%mpi(iproc)%AB%jexcldispl(jproc)+1
         is_e = sdata%mpi(iproc)%AB%jexcldispl(jproc)+sdata%mpi(iproc)%AB%jexclcounts(jproc)
         ir_s = sdata%mpi(jproc)%AB%jhalodispl(iproc)+1
         ir_e = sdata%mpi(jproc)%AB%jhalodispl(iproc)+sdata%mpi(jproc)%AB%jhalocounts(iproc)
         if (is_e-is_s/=ir_e-ir_s) call msgerror('wrong mpi_alltoallv indices')
         if (ir_e>=ir_s) bufr(jproc)%val(ir_s:ir_e) = bufs(iproc)%val(is_s:is_e)
      end do
   end do

   do iproc=1,sdata%nproc
      ! Copy zone A into zone B
      alpha(iproc)%valb(sdata%mpi(iproc)%isa_to_isb) = alpha(iproc)%vala

      ! Copy halo into zone B
      alpha(iproc)%valb(sdata%mpi(iproc)%AB%halo) = bufr(iproc)%val
   end do
end if

! Release memory
do iproc=1,sdata%nproc
   deallocate(alpha(iproc)%vala)
end do

end subroutine alpha_com_AB

!----------------------------------------------------------------------
! Subroutine: alpha_com_BA
!> Purpose: communicate reduced field from zone B to zone A
!----------------------------------------------------------------------
subroutine alpha_com_BA(sdata,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                 !< Sampling data
type(alphatype),intent(inout) :: alpha(sdata%nproc) !< Subgrid variable

! Local variables
integer :: iproc,jproc,ir_s,ir_e,is_s,is_e,i
type(buftype) :: bufs(sdata%nproc),bufr(sdata%nproc)

! Allocation
do iproc=1,sdata%nproc
   allocate(alpha(iproc)%vala(sdata%mpi(iproc)%nsa))
end do

if (sdata%nproc==1) then
   ! Copy data
   alpha(1)%vala = alpha(1)%valb
elseif (sdata%nproc>1) then
   ! Allocation
   do iproc=1,sdata%nproc
      allocate(bufs(iproc)%val(sdata%mpi(iproc)%AB%nhalo))
      allocate(bufr(iproc)%val(sdata%mpi(iproc)%AB%nexcl))
   end do

   ! Prepare buffers to send
   do iproc=1,sdata%nproc
      bufs(iproc)%val = alpha(iproc)%valb(sdata%mpi(iproc)%AB%halo)
   end do

   ! Emulate mpi_alltoallv
   do iproc=1,sdata%nproc
      do jproc=1,sdata%nproc
         ! Communication iproc->jproc
         is_s = sdata%mpi(iproc)%AB%jhalodispl(jproc)+1
         is_e = sdata%mpi(iproc)%AB%jhalodispl(jproc)+sdata%mpi(iproc)%AB%jhalocounts(jproc)
         ir_s = sdata%mpi(jproc)%AB%jexcldispl(iproc)+1
         ir_e = sdata%mpi(jproc)%AB%jexcldispl(iproc)+sdata%mpi(jproc)%AB%jexclcounts(iproc)
         if (is_e-is_s/=ir_e-ir_s) call msgerror('wrong mpi_alltoallv indices')
         if (ir_e>=ir_s) bufr(jproc)%val(ir_s:ir_e) = bufs(iproc)%val(is_s:is_e)
      end do
   end do

   do iproc=1,sdata%nproc
      ! Copy zone B into zone A
      alpha(iproc)%vala = alpha(iproc)%valb(sdata%mpi(iproc)%isa_to_isb)

      ! Copy halo into zone B
      do i=1,sdata%mpi(iproc)%AB%nexcl
         alpha(iproc)%vala(sdata%mpi(iproc)%AB%excl(i)) = alpha(iproc)%vala(sdata%mpi(iproc)%AB%excl(i))+bufr(iproc)%val(i)
      end do
   end do
end if

! Release memory
do iproc=1,sdata%nproc
   deallocate(alpha(iproc)%valb)
end do

end subroutine alpha_com_BA

!----------------------------------------------------------------------
! Subroutine: alpha_com_CA
!> Purpose: communicate reduced field from zone C to zone A
!----------------------------------------------------------------------
subroutine alpha_com_CA(sdata,alpha)

implicit none

! Passed variables
type(sdatatype),intent(in) :: sdata                 !< Sampling data
type(alphatype),intent(inout) :: alpha(sdata%nproc) !< Subgrid variable

! Local variables
integer :: iproc,jproc,ir_s,ir_e,is_s,is_e,i
type(buftype) :: bufs(sdata%nproc),bufr(sdata%nproc)

! Allocation
do iproc=1,sdata%nproc
   allocate(alpha(iproc)%vala(sdata%mpi(iproc)%nsa))
end do

if (sdata%nproc==1) then
   ! Copy data
   alpha(1)%vala = alpha(1)%valc
elseif (sdata%nproc>1) then
   ! Allocation
   do iproc=1,sdata%nproc
      allocate(bufs(iproc)%val(sdata%mpi(iproc)%AC%nhalo))
      allocate(bufr(iproc)%val(sdata%mpi(iproc)%AC%nexcl))
   end do

   ! Prepare buffers to send
   do iproc=1,sdata%nproc
      bufs(iproc)%val = alpha(iproc)%valc(sdata%mpi(iproc)%AC%halo)
   end do

   ! Emulate mpi_alltoallv
   do iproc=1,sdata%nproc
      do jproc=1,sdata%nproc
         ! Communication iproc->jproc
         is_s = sdata%mpi(iproc)%AC%jhalodispl(jproc)+1
         is_e = sdata%mpi(iproc)%AC%jhalodispl(jproc)+sdata%mpi(iproc)%AC%jhalocounts(jproc)
         ir_s = sdata%mpi(jproc)%AC%jexcldispl(iproc)+1
         ir_e = sdata%mpi(jproc)%AC%jexcldispl(iproc)+sdata%mpi(jproc)%AC%jexclcounts(iproc)
         if (is_e-is_s/=ir_e-ir_s) call msgerror('wrong mpi_alltoallv indices')
         if (ir_e>=ir_s) bufr(jproc)%val(ir_s:ir_e) = bufs(iproc)%val(is_s:is_e)
      end do
   end do

   do iproc=1,sdata%nproc
      ! Copy zone C into zone A
      alpha(iproc)%vala = alpha(iproc)%valc(sdata%mpi(iproc)%isa_to_isc)

      ! Copy halo into zone C
      do i=1,sdata%mpi(iproc)%AC%nexcl
         alpha(iproc)%vala(sdata%mpi(iproc)%AC%excl(i)) = alpha(iproc)%vala(sdata%mpi(iproc)%AC%excl(i))+bufr(iproc)%val(i)
      end do
   end do
end if

! Release memory
do iproc=1,sdata%nproc
   deallocate(alpha(iproc)%valc)
end do

end subroutine alpha_com_CA

!----------------------------------------------------------------------
! Subroutine: alpha_copy_AB
!> Purpose: copy reduced field from zone A to zone B
!----------------------------------------------------------------------
subroutine alpha_copy_AB(sdatampi,alpha)

implicit none

! Passed variables
type(sdatampitype),intent(in) :: sdatampi !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Allocation
allocate(alpha%valb(sdatampi%nsb))

! Initialize
alpha%valb = 0.0

! Copy zone A into zone B
alpha%valb(sdatampi%isa_to_isb) = alpha%vala

! Release memory
deallocate(alpha%vala)

end subroutine alpha_copy_AB

!----------------------------------------------------------------------
! Subroutine: alpha_copy_AC
!> Purpose: copy reduced field from zone A to zone C
!----------------------------------------------------------------------
subroutine alpha_copy_AC(sdatampi,alpha)

implicit none

! Passed variables
type(sdatampitype),intent(in) :: sdatampi !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Allocation
allocate(alpha%valc(sdatampi%nsc))

! Initialize
alpha%valc = 0.0

! Copy zone A into zone C
alpha%valc(sdatampi%isa_to_isc) = alpha%vala

! Release memory
deallocate(alpha%vala)

end subroutine alpha_copy_AC

!----------------------------------------------------------------------
! Subroutine: alpha_copy_BC
!> Purpose: copy reduced field from zone B to zone C
!----------------------------------------------------------------------
subroutine alpha_copy_BC(sdatampi,alpha)

implicit none

! Passed variables
type(sdatampitype),intent(in) :: sdatampi !< Sampling data
type(alphatype),intent(inout) :: alpha    !< Subgrid variable

! Allocation
allocate(alpha%valc(sdatampi%nsc))

! Initialization
alpha%valc = 0.0

! Copy zone B into zone C
alpha%valc(sdatampi%isb_to_isc) = alpha%valb

! Release memory
deallocate(alpha%valb)

end subroutine alpha_copy_BC

end module module_apply_com
