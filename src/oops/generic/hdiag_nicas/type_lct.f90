!----------------------------------------------------------------------
! Module: type_lct
!> Purpose: LCT data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_lct

use model_interface, only: model_write
use tools_const, only: req,reqkm
use hdiag_tools, only: diag_filter,diag_interpolation
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr,isallnotmsr
use type_hdata, only: hdatatype
use type_mpl, only: mpl

implicit none

! LCT data derived type
type lcttype
   integer :: npack                         !< Pack buffer size

   ! LCT structure
   integer :: nscales                       !< Number of LCT scales
   integer,allocatable :: ncomp(:)          !< Number of LCT components
   real(kind_real),allocatable :: H(:)      !< LCT components
   real(kind_real),allocatable :: coef(:)   !< LCT coefficients

   ! LCT fit
   real(kind_real),allocatable :: raw(:,:)  !< Raw correlations
   real(kind_real),allocatable :: norm(:,:) !< Norm to take nsub into account
   real(kind_real),allocatable :: fit(:,:)  !< Fitted correlations
end type lcttype

logical,parameter :: write_cor = .true.     !< Write raw and fitted correlations

interface lct_alloc
  module procedure lct_alloc_basic
  module procedure lct_alloc
end interface

private
public :: lcttype
public :: lct_alloc,lct_dealloc,lct_pack,lct_unpack,lct_write

contains

!----------------------------------------------------------------------
! Subroutine: lct_alloc_basic
!> Purpose: lct object basic allocation
!----------------------------------------------------------------------
subroutine lct_alloc_basic(hdata,lct)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata !< HDIAG data
type(lcttype),intent(inout) :: lct  !< LCT

! Local variables
integer :: iscales

! Associate
associate(nam=>hdata%nam)

! Number of scales and components
lct%nscales = nam%lct_nscales
allocate(lct%ncomp(lct%nscales))
do iscales=1,lct%nscales
   if (nam%lct_diag(iscales)) then
      lct%ncomp(iscales) = 3
   else
      lct%ncomp(iscales) = 4
   end if
end do

! Allocation
allocate(lct%H(sum(lct%ncomp)))
allocate(lct%coef(lct%nscales))

! Initialization
call msr(lct%H)
call msr(lct%coef)

! End associate
end associate

end subroutine lct_alloc_basic

!----------------------------------------------------------------------
! Subroutine: lct_alloc
!> Purpose: lct object allocation
!----------------------------------------------------------------------
subroutine lct_alloc(hdata,ib,lct)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata !< HDIAG data
integer,intent(in) :: ib            !< Block index
type(lcttype),intent(inout) :: lct  !< LCT

! Associate
associate(nam=>hdata%nam,bpar=>hdata%bpar)

! Basic allocation
call lct_alloc_basic(hdata,lct)

! Allocation
allocate(lct%raw(nam%nc3,bpar%nl0r(ib)))
allocate(lct%norm(nam%nc3,bpar%nl0r(ib)))
allocate(lct%fit(nam%nc3,bpar%nl0r(ib)))

! Initialization
lct%npack = sum(lct%ncomp)+lct%nscales+nam%nc3*bpar%nl0r(ib)
lct%raw = 0.0
lct%norm = 0.0
call msr(lct%fit)

! End associate
end associate

end subroutine lct_alloc

!----------------------------------------------------------------------
! Subroutine: lct_dealloc
!> Purpose: lct object deallocation
!----------------------------------------------------------------------
subroutine lct_dealloc(lct)

implicit none

! Passed variables
type(lcttype),intent(inout) :: lct !< LCT

! Release memory
if (allocated(lct%H)) deallocate(lct%H)
if (allocated(lct%coef)) deallocate(lct%coef)
if (allocated(lct%raw)) deallocate(lct%raw)
if (allocated(lct%norm)) deallocate(lct%norm)
if (allocated(lct%fit)) deallocate(lct%fit)

end subroutine lct_dealloc

!----------------------------------------------------------------------
! Subroutine: lct_pack
!> Purpose: LCT packing
!----------------------------------------------------------------------
subroutine lct_pack(hdata,ib,lct,buf)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata           !< HDIAG data
integer,intent(in) :: ib                      !< Block index
type(lcttype),intent(in) :: lct               !< LCT
real(kind_real),intent(out) :: buf(lct%npack) !< Buffer

! Local variables
integer :: offset

! Associate
associate(nam=>hdata%nam,bpar=>hdata%bpar)

! Pack
offset = 0
buf(offset+1:offset+sum(lct%ncomp)) = pack(lct%H,.true.)
offset = offset+sum(lct%ncomp)
buf(offset+1:offset+lct%nscales) = lct%coef
offset = offset+lct%nscales
buf(offset+1:offset+nam%nc3*bpar%nl0r(ib)) = pack(lct%fit,.true.)

! End associate
end associate

end subroutine lct_pack

!----------------------------------------------------------------------
! Subroutine: lct_unpack
!> Purpose: LCT unpacking
!----------------------------------------------------------------------
subroutine lct_unpack(hdata,ib,lct,buf)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata          !< HDIAG data
integer,intent(in) :: ib                     !< Block index
type(lcttype),intent(inout) :: lct           !< LCT
real(kind_real),intent(in) :: buf(lct%npack) !< Buffer

! Local variables
integer :: offset
logical,allocatable :: mask_unpack(:,:)

! Associate
associate(nam=>hdata%nam,bpar=>hdata%bpar)

! Allocation
allocate(mask_unpack(nam%nc3,bpar%nl0r(ib)))
mask_unpack = .true.

! Unpack
offset = 0
lct%H = buf(offset+1:offset+sum(lct%ncomp))
offset = offset+sum(lct%ncomp)
lct%coef = buf(offset+1:offset+lct%nscales)
offset = offset+lct%nscales
lct%fit = unpack(buf(offset+1:offset+nam%nc3*bpar%nl0r(ib)),mask_unpack,lct%fit)

! End associate
end associate

end subroutine lct_unpack

!----------------------------------------------------------------------
! Subroutine: lct_write
!> Purpose: interpolate and write LCT
!----------------------------------------------------------------------
subroutine lct_write(hdata,lct)

implicit none

! Passed variables
type(hdatatype),intent(inout) :: hdata                                      !< HDIAG data
type(lcttype),intent(in) :: lct(hdata%nam%nc1,hdata%geom%nl0,hdata%bpar%nb) !< LCT array

! Local variables
integer :: ib,iv,il0,jl0r,jl0,ic1,jc3,icomp,ic0,iscales,offset
real(kind_real) :: fac,det,rmse,rmse_count
real(kind_real),allocatable :: fld_c1(:,:,:),fld(:,:,:)
logical :: valid
character(len=1) :: iscaleschar

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

do ib=1,bpar%nb
   write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

   ! Initialization
   offset = 0

   do iscales=1,lct(1,1,ib)%nscales
      ! Allocation
      allocate(fld_c1(hdata%nam%nc1,hdata%geom%nl0,lct(1,1,ib)%ncomp(iscales)+1))
      allocate(fld(hdata%geom%nc0,hdata%geom%nl0,lct(1,1,ib)%ncomp(iscales)+2))

      do il0=1,geom%nl0
         write(mpl%unit,'(a10,a,i3,a)') '','Level ',nam%levs(il0),': '

         ! Initialization
         fac = 1.0
         valid = .false.

         do while (.not.valid)
            ! Copy LCT
            call msr(fld_c1(:,il0,:))
            do ic1=1,nam%nc1
               fld_c1(ic1,il0,1:lct(1,1,ib)%ncomp(iscales)) = lct(ic1,il0,ib)%H(offset+1:offset+lct(1,1,ib)%ncomp(iscales))
               fld_c1(ic1,il0,lct(1,1,ib)%ncomp(iscales)+1) = lct(ic1,il0,ib)%coef(iscales)
            end do

            ! Check invalid points
            valid = .true.
            do ic1=1,nam%nc1
               if (geom%mask(hdata%c1_to_c0(ic1),il0).and.(.not.isallnotmsr(fld_c1(ic1,il0,:)))) valid = .false.
            end do

            if (.not.valid) then
               ! Filter LCT
               write(mpl%unit,'(a13,a,f9.2,a)') '','Filter LCT with radius ',fac*nam%diag_rhflt*reqkm,' km'
               do icomp=1,lct(1,1,ib)%ncomp(iscales)+1
                  call diag_filter(hdata,il0,'median',fac*nam%diag_rhflt,fld_c1(:,il0,icomp))
                  call diag_filter(hdata,il0,'average',fac*nam%diag_rhflt,fld_c1(:,il0,icomp))
               end do

               ! Update fac (increase smoothing)
               fac = 2.0*fac
            end if
         end do
      end do

      ! Interpolate LCT
      write(mpl%unit,'(a10,a)') '','Interpolate LCT'
      do icomp=1,lct(1,1,ib)%ncomp(iscales)+1
         call diag_interpolation(hdata,fld_c1(:,:,icomp),fld(:,:,icomp))
      end do

      ! Compute horizontal length-scale
      do il0=1,geom%nl0
         do ic0=1,geom%nc0
            if (geom%mask(ic0,il0)) then
               ! Compute determinant
               if (lct(1,1,ib)%ncomp(iscales)==3) then
                  det = fld(ic0,il0,1)*fld(ic0,il0,2)
               else
                  det = fld(ic0,il0,1)*fld(ic0,il0,2)-fld(ic0,il0,4)**2
               end if

               ! Length-scale = determinant^{1/4}
               if (det>0.0) fld(ic0,il0,lct(1,1,ib)%ncomp(iscales)+2) = 1.0/sqrt(sqrt(det))
            end if
         end do
      end do

      if (mpl%main) then
         ! Write LCT
         write(mpl%unit,'(a10,a)') '','Write LCT'
         iv = bpar%b_to_v2(ib)
         write(iscaleschar,'(i1)') iscales
         call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc',trim(nam%varname(iv))//'_H11_'//iscaleschar, &
       & fld(:,:,1)/req**2)
         call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc',trim(nam%varname(iv))//'_H22_'//iscaleschar, &
       & fld(:,:,2)/req**2)
         call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc',trim(nam%varname(iv))//'_H33_'//iscaleschar, &
       & fld(:,:,3))
         if (lct(1,1,ib)%ncomp(iscales)==4) call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc', &
       & trim(nam%varname(iv))//'_Hc12_'//iscaleschar,fld(:,:,4))
         call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc',trim(nam%varname(iv))//'_coef_'//iscaleschar, &
       & fld(:,:,lct(1,1,ib)%ncomp(iscales)+1))
         call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc',trim(nam%varname(iv))//'_Lh_'//iscaleschar, &
       & fld(:,:,lct(1,1,ib)%ncomp(iscales)+2)*reqkm)
      end if

      ! Release memory
      deallocate(fld_c1)
      deallocate(fld)
   end do

   ! Compute RMSE
   rmse = 0.0
   rmse_count = 0.0
   do il0=1,geom%nl0
      do ic1=1,nam%nc1
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               if (hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) then
                  if (isnotmsr(lct(ic1,il0,ib)%fit(jc3,jl0))) then
                     rmse = rmse+(lct(ic1,il0,ib)%fit(jc3,jl0)-lct(ic1,il0,ib)%raw(jc3,jl0))**2
                     rmse_count = rmse_count+1.0
                  end if
               end if
            end do
         end do
      end do
   end do
   if (rmse_count>0.0) rmse = sqrt(rmse/rmse_count)
   write(mpl%unit,'(a7,a,e15.8,a,i8,a)') '','LCT diag RMSE: ',rmse,' for ',int(rmse_count),' diagnostic points'

   if (mpl%main.and.write_cor) then
      ! Allocation
      allocate(fld(hdata%geom%nc0,hdata%geom%nl0,2))

      ! Select level
      il0 = 1

      ! Write raw LCT
      write(mpl%unit,'(a7,a)') '','Write LCT diag'
      call msr(fld)
      do ic1=1,nam%nc1
         ! Check diagnostic area
         valid = .true.
         do jl0r=1,bpar%nl0r(ib)
            jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
            do jc3=1,nam%nc3
               if (valid.and.hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) &
            &  valid = valid.and.(.not.isnotmsr(fld(hdata%c1c3_to_c0(ic1,jc3),jl0,1)))
            end do
         end do
         if (valid) then
            do jl0r=1,bpar%nl0r(ib)
               jl0 = bpar%l0rl0b_to_l0(jl0r,il0,ib)
               do jc3=1,nam%nc3
                  if (hdata%c1l0_log(ic1,il0).and.hdata%c1c3l0_log(ic1,jc3,jl0)) then
                     fld(hdata%c1c3_to_c0(ic1,jc3),jl0,1) = lct(ic1,il0,ib)%raw(jc3,jl0)
                     fld(hdata%c1c3_to_c0(ic1,jc3),jl0,2) = lct(ic1,il0,ib)%fit(jc3,jl0)
                  end if
               end do
            end do
         end if
      end do
      iv = bpar%b_to_v2(ib)
      call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc',trim(nam%varname(iv))//'_raw',fld(:,:,1))
      call model_write(nam,geom,trim(nam%prefix)//'_lct_gridded.nc',trim(nam%varname(iv))//'_fit',fld(:,:,2))

      ! Release memory
      deallocate(fld)
   end if
end do

! End associate
end associate

end subroutine lct_write

end module type_lct

