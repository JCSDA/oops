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
use module_diag_tools, only: diag_filter,diag_interpolation
use tools_const, only: req,reqkm
use tools_kinds, only: kind_real
use tools_interp, only: compute_grid_interp_bilin
use tools_missing, only: msr,isnotmsr,isallnotmsr
use type_hdata, only: hdatatype
use type_mpl, only: mpl

implicit none

! LCT data derived type
type lcttype
   integer :: npack                         !< Pack buffer size

   ! LCT structure
   integer :: nscales                       !< Number of LCT scales
   integer :: ncomp                         !< Number of LCT components
   real(kind_real),allocatable :: H(:,:)    !< LCT components
   real(kind_real),allocatable :: coef(:)   !< LCT coefficients

   ! LCT fit
   real(kind_real),allocatable :: raw(:,:)  !< Raw correlations
   real(kind_real),allocatable :: norm(:,:) !< Norm to take nsub into account
   real(kind_real),allocatable :: fit(:,:)  !< Fitted correlations
end type lcttype

logical,parameter :: write_cor = .true. !< Write raw and fitted correlations

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
type(lcttype),intent(inout) :: lct  !< lct

! Associate
associate(nam=>hdata%nam)

! Number of scales and components
lct%nscales = nam%lct_nscales
if (nam%lct_diag) then
   lct%ncomp = 3
else
   lct%ncomp = 4
end if

! Allocation
allocate(lct%H(lct%nscales,lct%ncomp))
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
type(lcttype),intent(inout) :: lct  !< lct

! Associate
associate(nam=>hdata%nam,bpar=>hdata%bpar)

! Basic allocation
call lct_alloc_basic(hdata,lct)

! Allocation
allocate(lct%raw(nam%nc,bpar%nl0(ib)))
allocate(lct%norm(nam%nc,bpar%nl0(ib)))
allocate(lct%fit(nam%nc,bpar%nl0(ib)))

! Initialization
lct%npack = lct%nscales*lct%ncomp+lct%nscales+nam%nc*bpar%nl0(ib)
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
deallocate(lct%H)
deallocate(lct%coef)
deallocate(lct%raw)
deallocate(lct%norm)
deallocate(lct%fit)

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
buf(offset+1:offset+lct%nscales*lct%ncomp) = pack(lct%H,.true.)
offset = offset+lct%nscales*lct%ncomp
buf(offset+1:offset+lct%nscales) = lct%coef
offset = offset+lct%nscales
buf(offset+1:offset+nam%nc*bpar%nl0(ib)) = pack(lct%fit,.true.)

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
logical,allocatable :: mask_unpack_H(:,:),mask_unpack_fit(:,:)

! Associate
associate(nam=>hdata%nam,bpar=>hdata%bpar)

! Allocation
allocate(mask_unpack_H(lct%nscales,lct%ncomp))
allocate(mask_unpack_fit(nam%nc,bpar%nl0(ib)))
mask_unpack_H = .true.
mask_unpack_fit = .true.

! Unpack
offset = 0
lct%H = unpack(buf(offset+1:offset+lct%nscales*lct%ncomp),mask_unpack_H,lct%H)
offset = offset+lct%nscales*lct%ncomp
lct%coef = buf(offset+1:offset+lct%nscales)
offset = offset+lct%nscales
lct%fit = unpack(buf(offset+1:offset+nam%nc*bpar%nl0(ib)),mask_unpack_fit,lct%fit)

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
integer :: ib,iv,il0,jl0,il0r,ic1,ic,icomp,ic0,iscales
real(kind_real) :: fac,det,rmse,rmse_count
real(kind_real),allocatable :: fld_nc1(:,:,:),fld(:,:,:)
logical :: valid
character(len=1) :: iscaleschar

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

do ib=1,bpar%nb
   ! Allocation
   allocate(fld_nc1(hdata%nam%nc1,hdata%geom%nl0,lct(1,1,ib)%ncomp+1))
   allocate(fld(hdata%geom%nc0,hdata%geom%nl0,lct(1,1,ib)%ncomp+2))

   do iscales=1,lct(1,1,ib)%nscales
      ! Initialization
      fac = 1.0
      call msr(fld_nc1)
   
      do while (.not.isallnotmsr(fld_nc1))
         ! Copy LCT
         call msr(fld_nc1)
         do jl0=1,geom%nl0
            do ic1=1,nam%nc1
               fld_nc1(ic1,jl0,1:lct(1,1,ib)%ncomp) = lct(ic1,jl0,ib)%H(iscales,:)
               fld_nc1(ic1,jl0,lct(1,1,ib)%ncomp+1) = lct(ic1,jl0,ib)%coef(iscales)
            end do
         end do
   
         ! Filter LCT
         write(mpl%unit,'(a7,a,f9.2,a)') '','Filter LCT with radius ',fac*nam%diag_rhflt*reqkm,' km'
         hdata%nc2 = nam%nc1
         do icomp=1,lct(1,1,ib)%ncomp+1
            call diag_filter(hdata,'median',fac*nam%diag_rhflt,fld_nc1(:,:,icomp))
            call diag_filter(hdata,'average',fac*nam%diag_rhflt,fld_nc1(:,:,icomp))
         end do
   
         ! Update fac (increase smoothing)
         fac = 2.0*fac
      end do
   
      ! Interpolate LCT
      write(mpl%unit,'(a7,a)') '','Interpolate LCT'
      do icomp=1,lct(1,1,ib)%ncomp+1
         call diag_interpolation(hdata,fld_nc1(:,:,icomp),fld(:,:,icomp))
      end do
   
      ! Compute horizontal length-scale
      do il0=1,geom%nl0
         do ic0=1,geom%nc0
            if (geom%mask(ic0,il0)) then
               ! Compute determinant
               if (nam%lct_diag) then
                  det = fld(ic0,il0,1)*fld(ic0,il0,2)
               else
                  det = fld(ic0,il0,1)*fld(ic0,il0,2)-fld(ic0,il0,4)**2
               end if
    
               ! Length-scale = determinant^{1/4}
               if (det>0.0) fld(ic0,il0,lct(1,1,ib)%ncomp+2) = 1.0/sqrt(sqrt(det))
            end if
         end do
      end do

      if (mpl%main) then
         ! Write LCT
         write(mpl%unit,'(a7,a)') '','Write LCT'
         iv = bpar%ib_to_iv(ib)
         write(iscaleschar,'(i1)') iscales
         call model_write(nam,geom,trim(nam%prefix)//'_lct.nc',trim(nam%varname(iv))//'_H11_'//iscaleschar, &
       & fld(:,:,1)/req**2)
         call model_write(nam,geom,trim(nam%prefix)//'_lct.nc',trim(nam%varname(iv))//'_H22_'//iscaleschar, &
       & fld(:,:,2)/req**2)
         call model_write(nam,geom,trim(nam%prefix)//'_lct.nc',trim(nam%varname(iv))//'_H33_'//iscaleschar, &
       & fld(:,:,3))
         if (.not.nam%lct_diag) call model_write(nam,geom,trim(nam%prefix)//'_lct.nc', &
       & trim(nam%varname(iv))//'_Hc12_'//iscaleschar,fld(:,:,4))
         call model_write(nam,geom,trim(nam%prefix)//'_lct.nc',trim(nam%varname(iv))//'_coef_'//iscaleschar, &
       & fld(:,:,lct(1,1,ib)%ncomp+1))
         call model_write(nam,geom,trim(nam%prefix)//'_lct.nc',trim(nam%varname(iv))//'_Lh_'//iscaleschar, &
       & fld(:,:,lct(1,1,ib)%ncomp+2)*reqkm)
      end if
   end do

   ! Compute RMSE
   rmse = 0.0
   rmse_count = 0.0
   do jl0=1,geom%nl0
      do ic1=1,nam%nc1
         do il0r=1,bpar%nl0(ib)
            il0 = bpar%il0rjl0ib_to_il0(il0r,jl0,ib)
            do ic=1,nam%nc
               if (hdata%ic1il0_log(ic1,jl0).and.hdata%ic1icil0_log(ic1,ic,il0)) then
                  if (isnotmsr(lct(ic1,jl0,ib)%fit(ic,il0))) then
                     rmse = rmse+(lct(ic1,jl0,ib)%fit(ic,il0)-lct(ic1,jl0,ib)%raw(ic,il0))**2
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
      ! Select level
      jl0 = 1
   
      ! Write raw LCT
      write(mpl%unit,'(a7,a)') '','Write LCT diag'
      call msr(fld)
      do ic1=1,nam%nc1
         ! Check diagnostic area
         valid = .true.
         do il0r=1,bpar%nl0(ib)
            il0 = bpar%il0rjl0ib_to_il0(il0r,jl0,ib)
            do ic=1,nam%nc
               if (valid.and.hdata%ic1il0_log(ic1,jl0).and.hdata%ic1icil0_log(ic1,ic,il0)) &
            &  valid = valid.and.(.not.isnotmsr(fld(hdata%ic1icil0_to_ic0(ic1,ic,il0),il0,1)))
            end do
         end do
         if (valid) then
            do il0r=1,bpar%nl0(ib)
               il0 = bpar%il0rjl0ib_to_il0(il0r,jl0,ib)
               do ic=1,nam%nc
                  if (hdata%ic1il0_log(ic1,jl0).and.hdata%ic1icil0_log(ic1,ic,il0)) then
                     fld(hdata%ic1icil0_to_ic0(ic1,ic,il0),il0,1) = lct(ic1,jl0,ib)%raw(ic,il0)
                     fld(hdata%ic1icil0_to_ic0(ic1,ic,il0),il0,2) = lct(ic1,jl0,ib)%fit(ic,il0)
                  end if
               end do
            end do
         end if
      end do
      call model_write(nam,geom,trim(nam%prefix)//'_lct.nc',trim(nam%varname(iv))//'_raw',fld(:,:,1))
      call model_write(nam,geom,trim(nam%prefix)//'_lct.nc',trim(nam%varname(iv))//'_fit',fld(:,:,2))
   end if
   
   ! Release memory
   deallocate(fld_nc1)
   deallocate(fld)
end do

! End associate
end associate

end subroutine lct_write

end module type_lct

