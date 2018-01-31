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

use netcdf
use hdiag_tools, only: diag_filter,diag_interpolation,diag_com_lg
use tools_const, only: rad2deg
use tools_display, only: msgwarning,msgerror,prog_init,prog_print
use tools_kinds, only: kind_real
use tools_missing, only: msvalr,msr,isnotmsr,isallnotmsr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_bpar, only: bpartype
use type_curve, only: curvetype
use type_geom, only: geomtype
use type_hdata, only: hdatatype
use type_mpl, only: mpl,mpl_bcast
use type_nam, only: namtype,namncwrite

implicit none

! B data derived type
type bdatatype
   ! Block name
   character(len=1024) :: cname                     !< Block name

   ! Data
   real(kind_real),allocatable :: coef_ens(:,:)     !< Ensemble coefficient
   real(kind_real),allocatable :: coef_sta(:,:)     !< Static coefficient
   real(kind_real),allocatable :: rh0(:,:)          !< Fit support radius
   real(kind_real),allocatable :: rv0(:,:)          !< Fit support radius
   real(kind_real),allocatable :: rh0s(:,:)         !< Fit support radius  for sampling
   real(kind_real),allocatable :: rv0s(:,:)         !< Fit support radius, for sampling
   real(kind_real) :: wgt                           !< Block weight
   real(kind_real),allocatable :: lon_c0_flt(:,:,:) !< Displaced longitude
   real(kind_real),allocatable :: lat_c0_flt(:,:,:) !< Displaced latitude
   real(kind_real),allocatable :: trans(:,:)        !< Direct transform
   real(kind_real),allocatable :: transinv(:,:)     !< Inverse transform
end type bdatatype

interface diag_to_bdata
  module procedure diag_to_bdata
  module procedure diag_c2a_to_bdata
end interface

private
public :: bdatatype
public :: bdata_alloc,bdata_copy,diag_to_bdata,bdata_read,bdata_write

contains

!----------------------------------------------------------------------
! Subroutine: bdata_alloc
!> Purpose: bdata object allocation
!----------------------------------------------------------------------
subroutine bdata_alloc(nam,geom,bpar,bdata)

implicit none

! Passed variables
type(namtype),target,intent(in) :: nam              !< Namelist
type(geomtype),target,intent(in) :: geom            !< Geometry
type(bpartype),intent(in) :: bpar                   !< Block parameters
type(bdatatype),allocatable,intent(out) :: bdata(:) !< B data

! Local variables
integer :: ib

! Allocation
allocate(bdata(bpar%nb+1))

do ib=1,bpar%nb+1
   ! Set name
   bdata(ib)%cname = 'bdata_'//trim(bpar%blockname(ib))

   if (bpar%diag_block(ib)) then
      ! Allocation
      allocate(bdata(ib)%coef_ens(geom%nc0,geom%nl0))
      allocate(bdata(ib)%coef_sta(geom%nc0,geom%nl0))
      allocate(bdata(ib)%rh0(geom%nc0,geom%nl0))
      allocate(bdata(ib)%rv0(geom%nc0,geom%nl0))
      if (trim(nam%strategy)=='specific_multivariate') then
         allocate(bdata(ib)%rh0s(geom%nc0,geom%nl0))
         allocate(bdata(ib)%rv0s(geom%nc0,geom%nl0))
      end if
      if (nam%transform.and.bpar%auto_block(ib)) then
         allocate(bdata(ib)%trans(geom%nl0,geom%nl0))
         allocate(bdata(ib)%transinv(geom%nl0,geom%nl0))
      end if

      ! Initialization
      call msr(bdata(ib)%coef_ens)
      call msr(bdata(ib)%coef_sta)
      call msr(bdata(ib)%rh0)
      call msr(bdata(ib)%rv0)
      if (trim(nam%strategy)=='specific_multivariate') then
         call msr(bdata(ib)%rh0s)
         call msr(bdata(ib)%rv0s)
      end if
      call msr(bdata(ib)%wgt)
      if (nam%transform.and.bpar%auto_block(ib)) then
         call msr(bdata(ib)%trans)
         call msr(bdata(ib)%transinv)
      end if
   end if

   if ((ib==bpar%nb+1).and.nam%displ_diag) then
      ! Allocation
      allocate(bdata(ib)%lon_c0_flt(geom%nc0,geom%nl0,2:nam%nts))
      allocate(bdata(ib)%lat_c0_flt(geom%nc0,geom%nl0,2:nam%nts))

      ! Initialization
      if (nam%displ_diag) then
         call msr(bdata(ib)%lon_c0_flt)
         call msr(bdata(ib)%lat_c0_flt)
      end if
   end if
end do

end subroutine bdata_alloc

!----------------------------------------------------------------------
! Subroutine: bdata_copy
!> Purpose: bdata object copy
!----------------------------------------------------------------------
subroutine bdata_copy(bpar,bdata_in,bdata_out)

implicit none

! Passed variables
type(bpartype),intent(in) :: bpar                     !< Block parameter
type(bdatatype),intent(in) :: bdata_in(bpar%nb+1)     !< Input B data
type(bdatatype),intent(inout) :: bdata_out(bpar%nb+1) !< Output B data

! Local variables
integer :: ib

! Copy data
do ib=1,bpar%nb+1
   if (allocated(bdata_in(ib)%coef_ens)) bdata_out(ib)%coef_ens = bdata_in(ib)%coef_ens
   if (allocated(bdata_in(ib)%coef_sta)) bdata_out(ib)%coef_sta = bdata_in(ib)%coef_sta
   if (allocated(bdata_in(ib)%rh0)) bdata_out(ib)%rh0 = bdata_in(ib)%rh0
   if (allocated(bdata_in(ib)%rv0)) bdata_out(ib)%rv0 = bdata_in(ib)%rv0
   if (allocated(bdata_in(ib)%rh0s)) bdata_out(ib)%rh0s = bdata_in(ib)%rh0s
   if (allocated(bdata_in(ib)%rv0s)) bdata_out(ib)%rv0s = bdata_in(ib)%rv0s
   if (allocated(bdata_in(ib)%lon_c0_flt)) bdata_out(ib)%lon_c0_flt = bdata_in(ib)%lon_c0_flt
   if (allocated(bdata_in(ib)%lat_c0_flt)) bdata_out(ib)%lat_c0_flt = bdata_in(ib)%lat_c0_flt
   if (allocated(bdata_in(ib)%trans)) bdata_out(ib)%trans = bdata_in(ib)%trans
   if (allocated(bdata_in(ib)%transinv)) bdata_out(ib)%transinv = bdata_in(ib)%transinv
end do

end subroutine bdata_copy

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
type(curvetype),intent(in) :: diag(hdata%nc2a,hdata%bpar%nb+1) !< Diagnostics
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
subroutine bdata_read(nam,geom,bpar,bdata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                   !< Namelist
type(geomtype),intent(in) :: geom                 !< Geometry
type(bpartype),intent(in) :: bpar                 !< Block parameters
type(bdatatype),intent(inout) :: bdata(bpar%nb+1) !< B data

! Local variables
integer :: ib,il0,its
integer :: nc0_test,nl0_1_test,nts_test,nl0_2_test
integer :: info,ncid,nc0_id,nl0_1_id,nts_id,nl0_2_id
integer :: coef_ens_id,coef_sta_id,rh0_id,rv0_id,rh0s_id,rv0s_id,lon_c0_flt_id,lat_c0_flt_id,trans_id,transinv_id
character(len=1024) :: subr = 'bdata_read'

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      ! Open file
      info = nf90_open(trim(nam%datadir)//'/'//trim(nam%prefix)//'_'//trim(bdata(ib)%cname)//'.nc',nf90_nowrite,ncid)
      if (info==nf90_noerr) then
         ! Check dimensions
         call ncerr(subr,nf90_inq_dimid(ncid,'nc0',nc0_id))
         call ncerr(subr,nf90_inquire_dimension(ncid,nc0_id,len=nc0_test))
         call ncerr(subr,nf90_inq_dimid(ncid,'nl0_1',nl0_1_id))
         call ncerr(subr,nf90_inquire_dimension(ncid,nl0_1_id,len=nl0_1_test))
         if ((geom%nc0/=nc0_test).or.(geom%nl0/=nl0_1_test)) call msgerror('wrong dimension when reading B')
         if ((ib==bpar%nb+1).and.nam%displ_diag) then
            call ncerr(subr,nf90_inq_dimid(ncid,'nts',nts_id))
            call ncerr(subr,nf90_inquire_dimension(ncid,nts_id,len=nts_test))
            if (nam%nts-1/=nts_test) call msgerror('wrong dimension when reading B')
         end if
         if (nam%transform.and.bpar%auto_block(ib)) then
            call ncerr(subr,nf90_inq_dimid(ncid,'nl0_2',nl0_2_id))
            call ncerr(subr,nf90_inquire_dimension(ncid,nl0_2_id,len=nl0_2_test))
            if (geom%nl0/=nl0_2_test) call msgerror('wrong dimension when reading B')
         end if

         ! Get arrays ID
         if (bpar%nicas_block(ib)) then
            call ncerr(subr,nf90_inq_varid(ncid,'coef_ens',coef_ens_id))
            call ncerr(subr,nf90_inq_varid(ncid,'coef_sta',coef_sta_id))
            call ncerr(subr,nf90_inq_varid(ncid,'rh0',rh0_id))
            call ncerr(subr,nf90_inq_varid(ncid,'rv0',rv0_id))
            if (trim(nam%strategy)=='specific_multivariate') then
               call ncerr(subr,nf90_inq_varid(ncid,'rh0s',rh0s_id))
               call ncerr(subr,nf90_inq_varid(ncid,'rv0s',rv0s_id))
            end if
         end if
         if ((ib==bpar%nb+1).and.nam%displ_diag) then
            call ncerr(subr,nf90_inq_varid(ncid,'lon_c0_flt',lon_c0_flt_id))
            call ncerr(subr,nf90_inq_varid(ncid,'lat_c0_flt',lat_c0_flt_id))
         end if
         if (nam%transform.and.bpar%auto_block(ib)) then
            call ncerr(subr,nf90_inq_varid(ncid,'trans',trans_id))
            call ncerr(subr,nf90_inq_varid(ncid,'transinv',transinv_id))
         end if

         ! Read arrays
         if (bpar%nicas_block(ib)) then
            call ncerr(subr,nf90_get_var(ncid,coef_ens_id,bdata(ib)%coef_ens))
            call ncerr(subr,nf90_get_var(ncid,coef_sta_id,bdata(ib)%coef_sta))
            call ncerr(subr,nf90_get_var(ncid,rh0_id,bdata(ib)%rh0))
            call ncerr(subr,nf90_get_var(ncid,rv0_id,bdata(ib)%rv0))
            if (trim(nam%strategy)=='specific_multivariate') then
               call ncerr(subr,nf90_get_var(ncid,rh0s_id,bdata(ib)%rh0s))
               call ncerr(subr,nf90_get_var(ncid,rv0s_id,bdata(ib)%rv0s))
            end if
         end if
         if ((ib==bpar%nb+1).and.nam%displ_diag) then
            call ncerr(subr,nf90_get_var(ncid,lon_c0_flt_id,bdata(ib)%lon_c0_flt))
            call ncerr(subr,nf90_get_var(ncid,lat_c0_flt_id,bdata(ib)%lat_c0_flt))
         end if
         if (nam%transform.and.bpar%auto_block(ib)) then
            call ncerr(subr,nf90_get_var(ncid,trans_id,bdata(ib)%trans))
            call ncerr(subr,nf90_get_var(ncid,transinv_id,bdata(ib)%transinv))
         end if

         ! Get main weight
         call ncerr(subr,nf90_get_att(ncid,nf90_global,'wgt',bdata(ib)%wgt))

         ! Close file
         call ncerr(subr,nf90_close(ncid))
      else
         ! Use namelist/default values
         call msgwarning('cannot find B data to read, use namelist values')
         if (bpar%nicas_block(ib)) then
            bdata(ib)%coef_ens = 1.0
            bdata(ib)%coef_sta = 0.0
            do il0=1,geom%nl0
               bdata(ib)%rh0(:,il0) = nam%rh(il0)
               bdata(ib)%rv0(:,il0) = nam%rv(il0)
               if (trim(nam%strategy)=='specific_multivariate') then
                  bdata(ib)%rh0s(:,il0) = nam%rh(il0)
                  bdata(ib)%rv0s(:,il0) = nam%rv(il0)
               end if
            end do
         end if
         if ((ib==bpar%nb+1).and.nam%displ_diag) then
            do its=2,nam%nts
               do il0=1,geom%nl0
                  bdata(ib)%lon_c0_flt(:,il0,its) = geom%lon
                  bdata(ib)%lat_c0_flt(:,il0,its) = geom%lat
               end do
            end do
         end if
         if (nam%transform.and.bpar%auto_block(ib)) then
            bdata(ib)%trans = 0.0
            do il0=1,geom%nl0
               bdata(ib)%trans(il0,il0) = 1.0
            end do
            bdata(ib)%transinv = bdata(ib)%trans
         end if
         bdata(ib)%wgt = 1.0
      end if

      ! Check
      if (bpar%nicas_block(ib)) then
         if (any((bdata(ib)%rh0<0.0).and.isnotmsr(bdata(ib)%rh0))) call msgerror('rh0 should be positive')
         if (any((bdata(ib)%rv0<0.0).and.isnotmsr(bdata(ib)%rv0))) call msgerror('rv0 should be positive')
         if (trim(nam%strategy)=='specific_multivariate') then
            if (any((bdata(ib)%rh0s<0.0).and.isnotmsr(bdata(ib)%rh0s))) call msgerror('rh0s should be positive')
            if (any((bdata(ib)%rv0s<0.0).and.isnotmsr(bdata(ib)%rv0s))) call msgerror('rv0s should be positive')
         end if
      end if
   end if
end do

end subroutine bdata_read

!----------------------------------------------------------------------
! Subroutine: bdata_write
!> Purpose: write bdata object
!----------------------------------------------------------------------
subroutine bdata_write(nam,geom,bpar,bdata)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                !< Namelist
type(geomtype),intent(in) :: geom              !< Geometry
type(bpartype),intent(in) :: bpar              !< Block parameters
type(bdatatype),intent(in) :: bdata(bpar%nb+1) !< B data

! Local variables
integer :: ib
integer :: ncid,nc0_id,nl0_1_id,nl0_2_id,nts_id
integer :: lon_id,lat_id,coef_ens_id,coef_sta_id,rh0_id,rv0_id,rh0s_id,rv0s_id,trans_id,transinv_id,lon_c0_flt_id,lat_c0_flt_id
character(len=1024) :: subr = 'bdata_write'

do ib=1,bpar%nb+1
   if (bpar%B_block(ib)) then
      ! Processor verification
      if (.not.mpl%main) call msgerror('only I/O proc should enter '//trim(subr))

      ! Create file
      call ncerr(subr,nf90_create(trim(nam%datadir)//'/'//trim(nam%prefix)//'_'//trim(bdata(ib)%cname)//'.nc', &
       & or(nf90_clobber,nf90_64bit_offset),ncid))
      call namncwrite(nam,ncid)

      ! Define dimensions
      call ncerr(subr,nf90_def_dim(ncid,'nc0',geom%nc0,nc0_id))
      call ncerr(subr,nf90_def_dim(ncid,'nl0_1',geom%nl0,nl0_1_id))
      if (nam%transform.and.bpar%auto_block(ib)) call ncerr(subr,nf90_def_dim(ncid,'nl0_2',geom%nl0,nl0_2_id))
      if ((ib==bpar%nb+1).and.nam%displ_diag) call ncerr(subr,nf90_def_dim(ncid,'nts',nam%nts-1,nts_id))

      ! Define arrays
      if (bpar%nicas_block(ib)) then
         call ncerr(subr,nf90_def_var(ncid,'lon',ncfloat,(/nc0_id/),lon_id))
         call ncerr(subr,nf90_put_att(ncid,lon_id,'_FillValue',msvalr))
         call ncerr(subr,nf90_def_var(ncid,'lat',ncfloat,(/nc0_id/),lat_id))
         call ncerr(subr,nf90_put_att(ncid,lat_id,'_FillValue',msvalr))
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
      end if
      if ((ib==bpar%nb+1).and.nam%displ_diag) then
         call ncerr(subr,nf90_def_var(ncid,'lon_c0_flt',ncfloat,(/nc0_id,nl0_1_id,nts_id/),lon_c0_flt_id))
         call ncerr(subr,nf90_put_att(ncid,lon_c0_flt_id,'_FillValue',msvalr))
         call ncerr(subr,nf90_def_var(ncid,'lat_c0_flt',ncfloat,(/nc0_id,nl0_1_id,nts_id/),lat_c0_flt_id))
         call ncerr(subr,nf90_put_att(ncid,lat_c0_flt_id,'_FillValue',msvalr))
      end if
      if (nam%transform.and.bpar%auto_block(ib)) then
         call ncerr(subr,nf90_def_var(ncid,'trans',ncfloat,(/nl0_1_id,nl0_2_id/),trans_id))
         call ncerr(subr,nf90_put_att(ncid,trans_id,'_FillValue',msvalr))
         call ncerr(subr,nf90_def_var(ncid,'transinv',ncfloat,(/nl0_1_id,nl0_2_id/),transinv_id))
         call ncerr(subr,nf90_put_att(ncid,transinv_id,'_FillValue',msvalr))
      end if

      ! Write main weight
      call ncerr(subr,nf90_put_att(ncid,nf90_global,'wgt',bdata(ib)%wgt))

      ! End definition mode
      call ncerr(subr,nf90_enddef(ncid))

      ! Write arrays
      if (bpar%nicas_block(ib)) then
         call ncerr(subr,nf90_put_var(ncid,lon_id,geom%lon*rad2deg))
         call ncerr(subr,nf90_put_var(ncid,lat_id,geom%lat*rad2deg))
         call ncerr(subr,nf90_put_var(ncid,coef_ens_id,bdata(ib)%coef_ens))
         call ncerr(subr,nf90_put_var(ncid,coef_sta_id,bdata(ib)%coef_sta))
         call ncerr(subr,nf90_put_var(ncid,rh0_id,bdata(ib)%rh0))
         call ncerr(subr,nf90_put_var(ncid,rv0_id,bdata(ib)%rv0))
         if (trim(nam%strategy)=='specific_multivariate') then
            call ncerr(subr,nf90_put_var(ncid,rh0s_id,bdata(ib)%rh0s))
            call ncerr(subr,nf90_put_var(ncid,rv0s_id,bdata(ib)%rv0s))
         end if
      end if
      if ((ib==bpar%nb+1).and.nam%displ_diag) then
         call ncerr(subr,nf90_put_var(ncid,lon_c0_flt_id,bdata(ib)%lon_c0_flt))
         call ncerr(subr,nf90_put_var(ncid,lat_c0_flt_id,bdata(ib)%lat_c0_flt))
      end if
      if (nam%transform.and.bpar%auto_block(ib)) then
         call ncerr(subr,nf90_put_var(ncid,trans_id,bdata(ib)%trans))
         call ncerr(subr,nf90_put_var(ncid,transinv_id,bdata(ib)%transinv))
      end if

      ! Close file
      call ncerr(subr,nf90_close(ncid))
   end if
end do

end subroutine bdata_write

end module type_bdata
