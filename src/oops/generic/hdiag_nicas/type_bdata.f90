!----------------------------------------------------------------------
! Module: type_bdata
!> Purpose: sample data derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_bdata

use model_interface, only: model_write
use netcdf
use hdiag_tools, only: diag_filter,diag_interpolation,diag_com_lg
use tools_display, only: msgwarning,msgerror,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_curve, only: curvetype
use type_geom, only: geomtype
use type_hdata, only: hdatatype
use type_mpl, only: mpl,mpl_bcast
use type_nam, only: namtype

implicit none

! B data derived type
type bdatatype
   ! Block name
   character(len=1024) :: cname                 !< Block name

   ! Namelist
   type(namtype),pointer :: nam                 !< Namelist

   ! Geometry
   type(geomtype),pointer :: geom               !< Geometry

   ! Data
   real(kind_real),allocatable :: coef_ens(:,:) !< Ensemble coefficient
   real(kind_real),allocatable :: coef_sta(:,:) !< Static coefficient 
   real(kind_real),allocatable :: rh0(:,:)      !< Fit support radius
   real(kind_real),allocatable :: rv0(:,:)      !< Fit support radius
   real(kind_real),allocatable :: rh0s(:,:)     !< Fit support radius  for sampling
   real(kind_real),allocatable :: rv0s(:,:)     !< Fit support radius, for sampling
   real(kind_real) :: wgt                       !< Block weight

   ! Transforms
   real(kind_real),allocatable :: trans(:,:)    !< Direct transform
   real(kind_real),allocatable :: transinv(:,:) !< Inverse transform
end type bdatatype

interface diag_to_bdata
  module procedure diag_to_bdata
  module procedure diag_c2a_to_bdata
end interface

private
public :: bdatatype
public :: bdata_alloc,bdata_dealloc,diag_to_bdata,bdata_read,bdata_write

contains

!----------------------------------------------------------------------
! Subroutine: bdata_alloc
!> Purpose: bdata object allocation
!----------------------------------------------------------------------
subroutine bdata_alloc(bdata,auto_block)

implicit none

! Passed variables
type(bdatatype),intent(inout) :: bdata !< Sampling data
logical,intent(in) :: auto_block       !< Autocovariance block key

! Associate
associate(nam=>bdata%nam,geom=>bdata%geom)

! Allocation
allocate(bdata%coef_ens(geom%nc0,geom%nl0))
allocate(bdata%coef_sta(geom%nc0,geom%nl0))
allocate(bdata%rh0(geom%nc0,geom%nl0))
allocate(bdata%rv0(geom%nc0,geom%nl0))
if (trim(nam%strategy)=='specific_multivariate') then
   allocate(bdata%rh0s(geom%nc0,geom%nl0))
   allocate(bdata%rv0s(geom%nc0,geom%nl0))
end if
if (nam%transform.and.auto_block) then
   allocate(bdata%trans(geom%nl0,geom%nl0))
   allocate(bdata%transinv(geom%nl0,geom%nl0))
end if

! Initialization
call msr(bdata%coef_ens)
call msr(bdata%coef_sta)
call msr(bdata%rh0)
call msr(bdata%rv0)
if (trim(nam%strategy)=='specific_multivariate') then
   call msr(bdata%rh0s)
   call msr(bdata%rv0s)
end if
call msr(bdata%wgt)
if (nam%transform.and.auto_block) then
   call msr(bdata%trans)
   call msr(bdata%transinv)
end if

! End associate
end associate

end subroutine bdata_alloc

!----------------------------------------------------------------------
! Subroutine: bdata_dealloc
!> Purpose: bdata object deallocation
!----------------------------------------------------------------------
subroutine bdata_dealloc(bdata,auto_block)

implicit none

! Passed variables
type(bdatatype),intent(inout) :: bdata !< Sampling data
logical,intent(in) :: auto_block       !< Autocovariance block key

! Associate
associate(nam=>bdata%nam)

! Release memory
deallocate(bdata%coef_ens)
deallocate(bdata%coef_sta)
deallocate(bdata%rh0)
deallocate(bdata%rv0)
if (trim(nam%strategy)=='specific_multivariate') then
   deallocate(bdata%rh0s)
   deallocate(bdata%rv0s)
end if
if (nam%transform.and.auto_block) then
   deallocate(bdata%trans)
   deallocate(bdata%transinv)
end if

! End associate
end associate

end subroutine bdata_dealloc

!----------------------------------------------------------------------
! Subroutine: diag_to_bdata
!> Purpose: copy diagnostics into bdata object
!----------------------------------------------------------------------
subroutine diag_to_bdata(hdata,diag,bdata)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                     !< HDIAG data
type(curvetype),intent(in) :: diag(hdata%bpar%nb+1)     !< Diagnostics
type(bdatatype),intent(inout) :: bdata(hdata%bpar%nb+1) !< B data

! Local variables
integer :: ib,il0
real(kind_real) :: rh0s(hdata%geom%nl0),rv0s(hdata%geom%nl0)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Common sampling support radii
if (trim(nam%strategy)=='specific_multivariate') then
   ! Initialization
   rh0s = huge(1.0)
   rv0s = huge(1.0)

   ! Get minimum
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         do il0=1,geom%nl0
            rh0s(il0) = min(rh0s(il0),diag(ib)%fit_rh(il0))
            rv0s(il0) = min(rv0s(il0),diag(ib)%fit_rv(il0))
         end do
      end if
   end do

   ! Copy 
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         do il0=1,geom%nl0
            bdata(ib)%rh0s(:,il0) = rh0s(il0)
            bdata(ib)%rv0s(:,il0) = rv0s(il0)
         end do
      end if
   end do
end if

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      if (bpar%nicas_block(ib)) then
         do il0=1,geom%nl0
            bdata(ib)%coef_ens(:,il0) = diag(ib)%raw_coef_ens(il0)
            bdata(ib)%rh0(:,il0) = diag(ib)%fit_rh(il0)
            bdata(ib)%rv0(:,il0) = diag(ib)%fit_rv(il0)
            select case (trim(nam%method))
            case ('cor','loc')
               bdata(ib)%coef_sta(:,il0) = 0.0
            case ('hyb-avg','hyb-rnd')
               bdata(ib)%coef_sta(:,il0) = diag(ib)%raw_coef_sta
            case ('dual-ens')
               call msgerror('dual-ens not ready yet for B data')
            end select
         end do
      end if
      if (isanynotmsr(diag(ib)%raw_coef_ens)) then
         bdata(ib)%wgt = sum(diag(ib)%raw_coef_ens,mask=isnotmsr(diag(ib)%raw_coef_ens)) &
                       & /float(count(isnotmsr(diag(ib)%raw_coef_ens)))
      else
         call msgerror('missing weight for global B data')
      end if
   end if
end do

! End associate
end associate

end subroutine diag_to_bdata

!----------------------------------------------------------------------
! Subroutine: diag_c2a_to_bdata
!> Purpose: copy local diagnostics into bdata object
!----------------------------------------------------------------------
subroutine diag_c2a_to_bdata(hdata,diag,bdata)

implicit none

! Passed variables
type(hdatatype),intent(in) :: hdata                            !< HDIAG data
type(curvetype),intent(in) :: diag(hdata%nc2a,hdata%bpar%nb+1) !< Diagnostic curves
type(bdatatype),intent(inout) :: bdata(hdata%bpar%nb+1)        !< B data

! Local variables
integer :: ib,i,ic2a,il0,ic0
real(kind_real) :: fld_c0(hdata%geom%nc0,hdata%geom%nl0)
real(kind_real),allocatable :: fld_c2(:,:)
real(kind_real),allocatable :: rh0s(:,:),rv0s(:,:)

! Associate
associate(nam=>hdata%nam,geom=>hdata%geom,bpar=>hdata%bpar)

! Common sampling support radii (no filtering)
if (trim(nam%strategy)=='specific_multivariate') then
   ! Allocation
   allocate(rh0s(geom%nc0,geom%nl0))
   allocate(rv0s(geom%nc0,geom%nl0))
   
   ! Initialization
   rh0s = huge(1.0)
   rv0s = huge(1.0)

   ! Get minimum
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         do i=1,2
            ! Allocation
            allocate(fld_c2(hdata%nc2a,geom%nl0))

            ! Copy data
            do ic2a=1,hdata%nc2a
               if (i==1) then
                  fld_c2(ic2a,:) = diag(ic2a,ib)%fit_rh
               elseif (i==2) then
                  fld_c2(ic2a,:) = diag(ic2a,ib)%fit_rv
               end if
            end do

            ! Local to global
            call diag_com_lg(hdata,fld_c2)
            if (.not.mpl%main) allocate(fld_c2(hdata%nc2,geom%nl0))
            call mpl_bcast(fld_c2,mpl%ioproc)

            ! Median filter
            do il0=1,geom%nl0
               call diag_filter(hdata,il0,'median',nam%diag_rhflt,fld_c2(:,il0))
            end do

            ! Interpolate
            call diag_interpolation(hdata,fld_c2,fld_c0)

            ! Get minimum
            do il0=1,geom%nl0
               do ic0=1,geom%nc0
                  if (geom%mask(ic0,il0)) then
                     if (i==1) then
                        rh0s(ic0,il0) = min(rh0s(ic0,il0),fld_c0(ic0,il0))
                     elseif (i==2) then
                        rv0s(ic0,il0) = min(rv0s(ic0,il0),fld_c0(ic0,il0))
                     end if
                  end if
               end do
            end do

            ! Release memory
            deallocate(fld_c2)
         end do
      end if
   end do

   ! Copy 
   do ib=1,bpar%nb+1
      if (bpar%B_block(ib).and.bpar%nicas_block(ib)) then
         bdata(ib)%rh0s = rh0s
         bdata(ib)%rv0s = rv0s
      end if
   end do
end if

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      if (bpar%nicas_block(ib)) then
         do i=1,4
            ! Allocation
            allocate(fld_c2(hdata%nc2a,geom%nl0))

            ! Copy data
            do ic2a=1,hdata%nc2a
               if (i==1) then
                  fld_c2(ic2a,:) = diag(ic2a,ib)%raw_coef_ens
               elseif (i==2) then
                  select case (trim(nam%method))
                  case ('cor','loc')
                     fld_c2(ic2a,:) = 0.0
                  case ('hyb-avg','hyb-rnd')
                     fld_c2(ic2a,:) = diag(ic2a,ib)%raw_coef_sta
                  case ('dual-ens')
                     call msgerror('dual-ens not ready yet for B data')
                  end select
               elseif (i==3) then
                  fld_c2(ic2a,:) = diag(ic2a,ib)%fit_rh
               elseif (i==2) then
                  fld_c2(ic2a,:) = diag(ic2a,ib)%fit_rv
               end if
            end do

            ! Local to global
            call diag_com_lg(hdata,fld_c2)
            if (.not.mpl%main) allocate(fld_c2(hdata%nc2,geom%nl0))
            call mpl_bcast(fld_c2,mpl%ioproc)

            ! Median filter
            do il0=1,geom%nl0
               call diag_filter(hdata,il0,'median',nam%diag_rhflt,fld_c2(:,il0))
            end do

            ! Interpolate
            if (i==1) then
               call diag_interpolation(hdata,fld_c2,bdata(ib)%coef_ens)
               if (isanynotmsr(fld_c2)) then
                  bdata(ib)%wgt = sum(fld_c2,mask=isnotmsr(fld_c2))/float(count(isnotmsr(fld_c2)))
               else
                 call msgerror('missing weight for local B data')
               end if
            elseif (i==2) then
               call diag_interpolation(hdata,fld_c2,bdata(ib)%coef_sta)
            elseif (i==3) then
               call diag_interpolation(hdata,fld_c2,bdata(ib)%rh0)
            elseif (i==4) then
               call diag_interpolation(hdata,fld_c2,bdata(ib)%rv0)
            end if

            ! Release memory
            deallocate(fld_c2)
         end do
      else
         ! Allocation
         allocate(fld_c2(hdata%nc2a,geom%nl0))

         ! Copy data
         do ic2a=1,hdata%nc2a
            fld_c2(ic2a,:) = diag(ic2a,ib)%raw_coef_ens
         end do

         ! Local to global
         call diag_com_lg(hdata,fld_c2)
         call mpl_bcast(fld_c2,mpl%ioproc)

         ! Median filter
         do il0=1,geom%nl0
            call diag_filter(hdata,il0,'median',nam%diag_rhflt,fld_c2(:,il0))
         end do

         ! Copy data
         if (isanynotmsr(fld_c2)) then
            bdata(ib)%wgt = sum(fld_c2,mask=isnotmsr(fld_c2))/float(count(isnotmsr(fld_c2)))
         else
            call msgerror('missing weight for local B data')
         end if

         ! Release memory
         deallocate(fld_c2)
      end if
   end if
end do

! End associate
end associate

end subroutine diag_c2a_to_bdata

!----------------------------------------------------------------------
! Subroutine: bdata_read
!> Purpose: read bdata object
!----------------------------------------------------------------------
subroutine bdata_read(bdata,auto_block,nicas_block)

implicit none

! Passed variables
type(bdatatype),intent(inout) :: bdata !< B data
logical,intent(in) :: auto_block       !< Autocovariance block key
logical,intent(in) :: nicas_block      !< NICAS block key

! Local variables
integer :: nc0_test,nl0_1_test,nl0_2_test,il0
integer :: info,ncid,nc0_id,nl0_1_id,nl0_2_id
integer :: coef_ens_id,coef_sta_id,rh0_id,rv0_id,rh0s_id,rv0s_id,trans_id,transinv_id
character(len=1024) :: subr = 'bdata_read'

! Associate
associate(nam=>bdata%nam,geom=>bdata%geom)

! Open file
info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_'//trim(bdata%cname)//'.nc',nf90_nowrite,ncid)
if (info==nf90_noerr) then
   if (nicas_block) then
      ! Check dimensions
      call ncerr(subr,nf90_inq_dimid(ncid,'nc0',nc0_id))
      call ncerr(subr,nf90_inquire_dimension(ncid,nc0_id,len=nc0_test))
      call ncerr(subr,nf90_inq_dimid(ncid,'nl0_1',nl0_1_id))
      call ncerr(subr,nf90_inquire_dimension(ncid,nl0_1_id,len=nl0_1_test))
      if (nam%transform.and.auto_block) then
         call ncerr(subr,nf90_inq_dimid(ncid,'nl0_2',nl0_2_id))
         call ncerr(subr,nf90_inquire_dimension(ncid,nl0_2_id,len=nl0_2_test))
      end if
      if ((geom%nc0/=nc0_test).or.(geom%nl0/=nl0_1_test)) call msgerror('wrong dimension when reading B')
      if (nam%transform.and.auto_block) then
         if (geom%nl0/=nl0_2_test) call msgerror('wrong dimension when reading B')
      end if

      ! Get arrays ID
      call ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))
      call ncerr(subr,nf90_inq_varid(ncid,'coef_sta',coef_sta_id))
      call ncerr(subr,nf90_inq_varid(ncid,'rh0',rh0_id))
      call ncerr(subr,nf90_inq_varid(ncid,'rv0',rv0_id))
      if (trim(nam%strategy)=='specific_multivariate') then
         call ncerr(subr,nf90_inq_varid(ncid,'rh0s',rh0s_id))
         call ncerr(subr,nf90_inq_varid(ncid,'rv0s',rv0s_id))
      end if
      if (nam%transform.and.auto_block) then
         call ncerr(subr,nf90_inq_varid(ncid,'trans',trans_id))
         call ncerr(subr,nf90_inq_varid(ncid,'transinv',transinv_id))
      end if

      ! Read arrays
      call ncerr(subr,nf90_get_var(ncid,coef_ens_id,bdata%coef_ens))
      call ncerr(subr,nf90_get_var(ncid,coef_sta_id,bdata%coef_sta))
      call ncerr(subr,nf90_get_var(ncid,rh0_id,bdata%rh0))
      call ncerr(subr,nf90_get_var(ncid,rv0_id,bdata%rv0))
      if (trim(nam%strategy)=='specific_multivariate') then
         call ncerr(subr,nf90_get_var(ncid,rh0s_id,bdata%rh0s))
         call ncerr(subr,nf90_get_var(ncid,rv0s_id,bdata%rv0s))
      end if
      if (nam%transform.and.auto_block) then
         call ncerr(subr,nf90_get_var(ncid,trans_id,bdata%trans))
         call ncerr(subr,nf90_get_var(ncid,transinv_id,bdata%transinv))
      end if
   end if

   ! Get main weight
   call ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',bdata%wgt))

   ! Close file
   call ncerr(subr,nf90_close(ncid))
else
   call msgwarning('cannot find B data to read, use namelist values')
   if (nicas_block) then
      bdata%coef_ens = 1.0
      bdata%coef_sta = 0.0
      do il0=1,geom%nl0
         bdata%rh0(:,il0) = nam%rh(il0)
         bdata%rv0(:,il0) = nam%rv(il0)
         if (trim(nam%strategy)=='specific_multivariate') then
            bdata%rh0s(:,il0) = nam%rh(il0)
            bdata%rv0s(:,il0) = nam%rv(il0)
         end if
      end do
      if (nam%transform.and.auto_block) then
         bdata%trans = 0.0
         do il0=1,geom%nl0
            bdata%trans(il0,il0) = 1.0
         end do
         bdata%transinv = bdata%trans
      end if
   end if
   bdata%wgt = 1.0
end if

! Check
if (any((bdata%rh0<0.0).and.isnotmsr(bdata%rh0))) call msgerror('rh0 should be positive')
if (any((bdata%rv0<0.0).and.isnotmsr(bdata%rv0))) call msgerror('rv0 should be positive')
if (trim(nam%strategy)=='specific_multivariate') then
   if (any((bdata%rh0s<0.0).and.isnotmsr(bdata%rh0s))) call msgerror('rh0s should be positive')
   if (any((bdata%rv0s<0.0).and.isnotmsr(bdata%rv0s))) call msgerror('rv0s should be positive')
end if

! End associate
end associate

end subroutine bdata_read

!----------------------------------------------------------------------
! Subroutine: bdata_write
!> Purpose: write bdata object
!----------------------------------------------------------------------
subroutine bdata_write(bdata,auto_block,nicas_block)

implicit none

! Passed variables
type(bdatatype),intent(in) :: bdata !< B data
logical,intent(in) :: auto_block    !< Autocovariance block key
logical,intent(in) :: nicas_block   !< NICAS block key

! Local variables
integer :: ncid,nc0_id,nl0_1_id,nl0_2_id
integer :: coef_ens_id,coef_sta_id,rh0_id,rv0_id,rh0s_id,rv0s_id,trans_id,transinv_id
character(len=1024) :: subr = 'bdata_write'

! Associate
associate(nam=>bdata%nam,geom=>bdata%geom)

! Processor verification
if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

! Create file
call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_'//trim(bdata%cname)//'.nc', &
 & or(nf90_clobber,nf90_64bit_offset),ncid))

if (nicas_block) then
   ! Define dimensions
   call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
   call ncerr(subr,nf90_def_dim(ncid,'nl0_1',geom%nl0,nl0_1_id))
   if (nam%transform.and.auto_block) call ncerr(subr,nf90_def_dim(ncid,'nl0_2',geom%nl0,nl0_2_id))

   ! Define arrays
   call ncerr(subr,nf90_def_var(ncid,'coef_ens',ncfloat,(/nc0_id,nl0_1_id/),coef_ens_id))
   call ncerr(subr,nf90_put_att(ncid,coef_ens_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'coef_sta',ncfloat,(/nc0_id,nl0_1_id/),coef_sta_id))
   call ncerr(subr,nf90_put_att(ncid,coef_sta_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'rh0',ncfloat,(/nc0_id,nl0_1_id/),rh0_id))
   call ncerr(subr,nf90_put_att(ncid,rh0_id,'_FillValue',msvalr))
   call ncerr(subr,nf90_def_var(ncid,'rv0',ncfloat,(/nc0_id,nl0_1_id/),rv0_id))
   call ncerr(subr,nf90_put_att(ncid,rv0_id,'_FillValue',msvalr))
   if (trim(nam%strategy)=='specific_multivariate') then
      call ncerr(subr,nf90_def_var(ncid,'rh0s',ncfloat,(/nc0_id,nl0_1_id/),rh0s_id))
      call ncerr(subr,nf90_put_att(ncid,rh0s_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_def_var(ncid,'rv0s',ncfloat,(/nc0_id,nl0_1_id/),rv0s_id))
      call ncerr(subr,nf90_put_att(ncid,rv0s_id,'_FillValue',msvalr))
   end if

   if (nam%transform.and.auto_block) then
      call ncerr(subr,nf90_def_var(ncid,'trans',ncfloat,(/nl0_1_id,nl0_2_id/),trans_id))
      call ncerr(subr,nf90_put_att(ncid,trans_id,'_FillValue',msvalr))
      call ncerr(subr,nf90_def_var(ncid,'transinv',ncfloat,(/nl0_1_id,nl0_2_id/),transinv_id))
      call ncerr(subr,nf90_put_att(ncid,transinv_id,'_FillValue',msvalr))
   end if
end if

! Write main weight
call ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',bdata%wgt))

! End definition mode
call ncerr(subr,nf90_enddef(ncid))

! Write arrays
if (nicas_block) then
   call ncerr(subr,nf90_put_var(ncid,coef_ens_id,bdata%coef_ens))
   call ncerr(subr,nf90_put_var(ncid,coef_sta_id,bdata%coef_sta))
   call ncerr(subr,nf90_put_var(ncid,rh0_id,bdata%rh0))
   call ncerr(subr,nf90_put_var(ncid,rv0_id,bdata%rv0))
   if (trim(nam%strategy)=='specific_multivariate') then
      call ncerr(subr,nf90_put_var(ncid,rh0s_id,bdata%rh0s))
      call ncerr(subr,nf90_put_var(ncid,rv0s_id,bdata%rv0s))
   end if
   if (nam%transform.and.auto_block) then
      call ncerr(subr,nf90_put_var(ncid,trans_id,bdata%trans))
      call ncerr(subr,nf90_put_var(ncid,transinv_id,bdata%transinv))
   end if
end if

! Close file
call ncerr(subr,nf90_close(ncid))

! Write gridded data (for visualisation)
if (nicas_block) then
   call model_write(nam,geom,trim(nam%prefix)//'_gridded_'//trim(bdata%cname)//'.nc','coef_ens',bdata%coef_ens)
   call model_write(nam,geom,trim(nam%prefix)//'_gridded_'//trim(bdata%cname)//'.nc','rh0',bdata%rh0)
   call model_write(nam,geom,trim(nam%prefix)//'_gridded_'//trim(bdata%cname)//'.nc','rv0',bdata%rv0)
   if (trim(nam%strategy)=='specific_multivariate') then
      call model_write(nam,geom,trim(nam%prefix)//'_gridded_'//trim(bdata%cname)//'.nc','rh0s',bdata%rh0s)
      call model_write(nam,geom,trim(nam%prefix)//'_gridded_'//trim(bdata%cname)//'.nc','rv0s',bdata%rv0s)
   end if
   call model_write(nam,geom,trim(nam%prefix)//'_gridded_'//trim(bdata%cname)//'.nc','coef_sta',bdata%coef_sta)
end if

! End associate
end associate

end subroutine bdata_write

end module type_bdata
