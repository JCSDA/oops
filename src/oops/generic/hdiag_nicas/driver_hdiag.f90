!----------------------------------------------------------------------
! Module: driver_hdiag
!> Purpose: hybrid_diag driver
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module driver_hdiag

use model_interface, only: model_write
use module_average, only: compute_avg,compute_avg_lr,compute_avg_asy,compute_bwgtsq,compute_bwavg
use module_displacement, only: compute_displacement
use module_dualens, only: compute_dualens
use module_fit, only: compute_fit
use module_hybridization, only: compute_hybridization
use module_localization, only: compute_localization
use module_moments, only: compute_moments
use module_sampling, only: setup_sampling
use netcdf
use tools_const, only: reqkm
use tools_display, only: vunitchar,prog_init,prog_print,msgerror,msgwarning,aqua,aqua,peach,peach,purple,purple,black
use tools_kinds, only: kind_real
use tools_missing, only: msvali,msvalr,msr,isnotmsi,isnotmsr,isanynotmsr
use tools_nc, only: ncerr,ncfloat
use type_avg, only: avgtype,avg_dealloc
use type_bdata, only: bdatatype,bdata_alloc,diag_to_bdata,bdata_read,bdata_write
use type_bpar, only: bpartype
use type_curve, only: curvetype,curve_alloc,curve_dealloc,curve_normalization,curve_write,curve_write_all,curve_write_local
use type_displ, only: displtype,displ_alloc,displ_dealloc,displ_write
use type_geom, only: geomtype
use type_hdata, only: hdatatype
use type_mom, only: momtype,mom_dealloc
use type_mpl, only: mpl
use type_nam, only: namtype

implicit none

private
public :: run_hdiag

contains

!----------------------------------------------------------------------
! Subroutine: run_hdiag
!> Purpose: hybrid_diag
!----------------------------------------------------------------------
subroutine run_hdiag(nam,geom,bpar,bdata,ens1)

implicit none

! Passed variables
type(namtype),target,intent(inout) :: nam                                                  !< Namelist
type(geomtype),target,intent(in) :: geom                                                   !< Geometry
type(bpartype),target,intent(in) :: bpar                                                   !< Block parameters
type(bdatatype),allocatable,intent(inout) :: bdata(:)                                      !< B data
real(kind_real),intent(in),optional :: ens1(geom%nc0a,geom%nl0,nam%nv,nam%nts,nam%ens1_ne) !< Ensemble 1

! Local variables
integer :: ib,jl0,ic2,ildw
character(len=7) :: lonchar,latchar
character(len=1024) :: filename
type(avgtype) :: avg_1(bpar%nb+1),avg_2(bpar%nb+1)
type(avgtype),allocatable :: avg_1_nc2(:,:)
type(curvetype) :: cor_1(bpar%nb+1),cor_2(bpar%nb+1)
type(curvetype),allocatable :: cor_1_nc2(:,:)
type(curvetype) :: loc_1(bpar%nb+1),loc_2(bpar%nb+1)
type(curvetype) :: loc_3(bpar%nb+1),loc_4(bpar%nb+1)
type(curvetype),allocatable :: loc_1_nc2(:,:)
type(displtype) :: displ
type(hdatatype) :: hdata
type(momtype) :: mom_1(bpar%nb),mom_2(bpar%nb)

! Allocate B data
allocate(bdata(bpar%nb+1))
do ib=1,bpar%nb+1
   if (bpar%diag_block(ib)) then
      bdata(ib)%nam => nam
      bdata(ib)%geom => geom
      bdata(ib)%cname = 'bdata_'//trim(bpar%blockname(ib))
      call bdata_alloc(bdata(ib))
   end if
end do

if (nam%new_hdiag) then
   ! Set pointers
   hdata%nam => nam
   hdata%geom => geom
   hdata%bpar => bpar

   ! Setup sampling
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,i5,a)') '--- Setup sampling (nc1 = ',nam%nc1,')'

   call setup_sampling(hdata)

   if (nam%displ_diag) then
      ! Compute displacement diagnostic
      write(*,'(a)') '-------------------------------------------------------------------'
      write(*,'(a)') '--- Compute displacement diagnostic'

      if (present(ens1)) then
         call compute_displacement(hdata,displ,ens1)
      else
        call compute_displacement(hdata,displ)
      end if
   end if

   ! Compute sample moments
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute sample moments'

   ! Compute ensemble 1 sample moments
   write(mpl%unit,'(a7,a)') '','Ensemble 1:'
   if (present(ens1)) then
      call compute_moments(hdata,'ens1',displ,mom_1,ens1)
   else
      call compute_moments(hdata,'ens1',displ,mom_1)
   end if

   if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) then
      ! Compute randomized sample moments
      write(mpl%unit,'(a7,a)') '','Ensemble 2:'
      call compute_moments(hdata,'ens2',displ,mom_2)
   end if

   ! Compute statistics
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a)') '--- Compute statistics'

   ! Allocation
   if (nam%local_diag) then
      allocate(avg_1_nc2(hdata%nc2,bpar%nb+1))
      allocate(cor_1_nc2(hdata%nc2,bpar%nb+1))
   end if

   do ib=1,bpar%nb
      if (bpar%diag_block(ib)) then
         write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))
         write(mpl%unit,'(a7,a)') '','Compute global statistics'

         ! Compute global statistics
         call compute_avg(hdata,ib,mom_1(ib),avg_1(ib))

         ! Compute global asymptotic statistics
         call compute_avg_asy(hdata,ib,nam%ne,avg_1(ib))

         if (nam%local_diag) then
            write(mpl%unit,'(a7,a)') '','Compute local statistics'

            ! Compute local statistics
            call compute_avg(hdata,ib,mom_1(ib),avg_1_nc2(:,ib))

            ! Compute local asymptotic statistics
            call compute_avg_asy(hdata,ib,nam%ne,avg_1_nc2(:,ib))
         end if

         ! Print results
         do jl0=1,geom%nl0
            write(mpl%unit,'(a10,a,i3,a4,a21,a,e9.2,a,a,a,f8.2,a)') '','Level: ',nam%levs(jl0),' ~> ', &
             & 'raw cov. / cor. (1): ',trim(peach),avg_1(ib)%m11(1,bpar%il0rz(jl0,ib),jl0),trim(black),' / ', &
             & trim(peach),avg_1(ib)%cor(1,bpar%il0rz(jl0,ib),jl0),trim(black)
         end do

         if (trim(nam%method)=='hyb-avg') then
            ! Static covariance = ensemble covariance
            avg_1(ib)%m11sta = avg_1(ib)%m11*avg_1(ib)%m11
            avg_1(ib)%stasq = avg_1(ib)%m11**2
         elseif (trim(nam%method)=='hyb-rnd') then
            ! Compute randomized averaged statistics
            call compute_avg(hdata,ib,mom_2(ib),avg_2(ib))

            ! Static covariance = randomized covariance
            avg_1(ib)%m11sta = avg_1(ib)%m11*avg_2(ib)%m11
            avg_1(ib)%stasq = avg_2(ib)%m11**2
         elseif (trim(nam%method)=='dual-ens') then
            ! Compute low-resolution averaged statistics
            call compute_avg(hdata,ib,mom_2(ib),avg_2(ib))
            call compute_avg_asy(hdata,ib,nam%ens2_ne,avg_2(ib))

            ! LR covariance/HR covariance product average
            call compute_avg_lr(hdata,ib,mom_1(ib),mom_2(ib),avg_1(ib),avg_2(ib))
         end if

         ! Print results
         select case (trim(nam%method))
         case ('hyb-avg','hyb-rnd','dual-ens')
            do jl0=1,geom%nl0
               write(mpl%unit,'(a10,a,i3,a4,a21,a,e9.2,a,a,a,f8.2,a)') '','Level: ',nam%levs(jl0),' ~> ', &
                & 'raw cov. / cor. (2): ',trim(peach),avg_2(ib)%m11(1,bpar%il0rz(jl0,ib),jl0),trim(black),' / ', &
                & trim(peach),avg_2(ib)%cor(1,bpar%il0rz(jl0,ib),jl0),trim(black)
            end do
         end select
      end if
   end do

   if (bpar%diag_block(bpar%nb+1)) then
      ! Define block weights from avg_1
      call compute_bwgtsq(hdata,avg_1(1:bpar%nb))

      ! Compute global block averages
      call compute_bwavg(hdata,avg_1)
      if ((trim(nam%method)=='hyb-rnd').or.(trim(nam%method)=='dual-ens')) call compute_bwavg(hdata,avg_2)
      if (nam%local_diag) then
         ! Compute local block averages
         call compute_bwavg(hdata,avg_1_nc2)
      end if
   end if

   ! Copy correlation
   do ib=1,bpar%nb+1
      if (bpar%diag_block(ib)) then
         ! Allocation
         call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_cor',cor_1(ib))
         select case (trim(nam%method))
         case ('hyb-avg','hyb-rnd')
            call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_cor_hyb',cor_2(ib))
         case ('dual-ens')
            call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_cor_lr',cor_2(ib))
         end select

         ! Copy
         cor_1(ib)%raw = avg_1(ib)%cor
         call curve_normalization(hdata,ib,cor_1(ib))
         select case (trim(nam%method))
         case ('hyb-avg')
            cor_2(ib)%raw = avg_1(ib)%cor
            call curve_normalization(hdata,ib,cor_2(ib))
         case ('hyb-rnd','dual-ens')
            cor_2(ib)%raw = avg_2(ib)%cor
            call curve_normalization(hdata,ib,cor_2(ib))
         end select

         if (nam%local_diag) then
            do ic2=1,hdata%nc2
               ! Allocation
               call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_cor_nc2',cor_1_nc2(ic2,ib))

               ! Copy
               cor_1_nc2(ic2,ib)%raw = avg_1_nc2(ic2,ib)%cor
               cor_1_nc2(ic2,ib)%raw_coef_ens = 1.0
            end do
         end if
      end if
   end do

   if (any(bpar%fit_block)) then
      ! Compute first ensemble correlation fit
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute first ensemble correlation fit'

      do ib=1,bpar%nb+1
         if (bpar%fit_block(ib)) then
            write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

            ! Compute global fit
            write(mpl%unit,'(a7,a)') '','Compute global fit'
            call compute_fit(hdata,ib,cor_1(ib))

            if (nam%local_diag) then
               ! Compute local fit
               write(mpl%unit,'(a7,a)',advance='no') '','Compute local fit:'
               call compute_fit(hdata,ib,cor_1_nc2(:,ib))
            end if

            ! Print results
             do jl0=1,geom%nl0
                write(mpl%unit,'(a10,a,i3,a,f8.2,a,f8.2,a)') '','Level: ',nam%levs(jl0), &
              & ' ~> cor. support radii: '//trim(aqua),cor_1(ib)%fit_rh(jl0)*reqkm,trim(black)//' km  / ' &
              & //trim(aqua),cor_1(ib)%fit_rv(jl0),trim(black)//' '//trim(vunitchar)
             end do
         end if
      end do

      select case (trim(nam%method))
      case ('hyb-avg','hyb-rnd','dual-ens')
         ! Compute second ensemble correlation fit
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a)') '--- Compute second ensemble correlation fit'

         do ib=1,bpar%nb+1
            if (bpar%fit_block(ib)) then
               write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

               ! Compute global fit
               write(mpl%unit,'(a7,a)') '','Compute global fit'
               call compute_fit(hdata,ib,cor_2(ib))

               ! Print results
               do jl0=1,geom%nl0
                   write(mpl%unit,'(a10,a,i3,a,f8.2,a,f8.2,a)') '','Level: ',nam%levs(jl0), &
                 & ' ~> cor. support radii: '//trim(aqua),cor_2(ib)%fit_rh(jl0)*reqkm,trim(black)//' km  / ' &
                 & //trim(aqua),cor_2(ib)%fit_rv(jl0),trim(black)//' '//trim(vunitchar)
                end do
            end if
         end do
      end select
   end if

   select case (trim(nam%method))
   case ('loc','hyb-avg','hyb-rnd','dual-ens')
      ! Compute localization diagnostic and fit
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute localization diagnostic and fit'

      ! Allocation
      if (nam%local_diag) allocate(loc_1_nc2(hdata%nc2,bpar%nb+1))

      do ib=1,bpar%nb+1
         if (bpar%diag_block(ib)) then
            write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

            ! Allocation
            call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_loc',loc_1(ib))
            if (nam%local_diag) then
               do ic2=1,hdata%nc2
                  call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_loc_nc2',loc_1_nc2(ic2,ib))
               end do
            end if

            ! Compute global localization
            write(mpl%unit,'(a7,a)') '','Compute global localization'
            call compute_localization(hdata,ib,avg_1(ib),loc_1(ib))

            if (nam%local_diag) then
               ! Compute local localization
               write(mpl%unit,'(a7,a)',advance='no') '','Compute local localization:'
               call compute_localization(hdata,ib,avg_1_nc2(:,ib),loc_1_nc2(:,ib))
            end if

            ! Print results
            do jl0=1,geom%nl0
               write(mpl%unit,'(a10,a,i3,a4,a21,a,f8.2,a)') '','Level: ',nam%levs(jl0),' ~> ', &
             & 'diagonal value: ',trim(peach),loc_1(ib)%raw_coef_ens(jl0),trim(black)
               if (bpar%fit_block(ib)) then
                  write(mpl%unit,'(a45,a,f8.2,a,f8.2,a)') 'loc. support radii: ',trim(aqua),loc_1(ib)%fit_rh(jl0)*reqkm, &
                & trim(black)//' km  / '//trim(aqua),loc_1(ib)%fit_rv(jl0),trim(black)//' '//trim(vunitchar)
               end if
            end do
         end if
      end do
   end select

   select case (trim(nam%method))
   case ('hyb-avg','hyb-rnd')
      ! Compute static hybridization diagnostic and fit
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute static hybridization diagnostic and fit'

      do ib=1,bpar%nb+1
         if (bpar%diag_block(ib)) then
            write(mpl%unit,'(a7,a,a)') '','Block: ',trim(bpar%blockname(ib))

            ! Allocation
            call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_loc_hyb',loc_2(ib))

            ! Compute global static hybridization
            write(mpl%unit,'(a7,a)') '','Compute global static hybridization'
            call compute_hybridization(hdata,ib,avg_1(ib),loc_2(ib))

            ! Print results
            do jl0=1,geom%nl0
               write(mpl%unit,'(a10,a,i3,a4,a21,a,f8.2,a)') '','Level: ',nam%levs(jl0),' ~> ', &
             & 'diagonal value: ',trim(peach),loc_2(ib)%raw_coef_ens(jl0),trim(black)
               if (bpar%fit_block(ib)) then
                  write(mpl%unit,'(a45,a,f8.2,a,f8.2,a)') 'loc. support radii: ',trim(aqua),loc_2(ib)%fit_rh(jl0)*reqkm, &
                & trim(black)//' km  / '//trim(aqua),loc_2(ib)%fit_rv(jl0),trim(black)//' '//trim(vunitchar)
               end if
            end do
            write(mpl%unit,'(a10,a,f8.2,a)') '','Raw static coeff.: ',trim(purple),loc_2(ib)%raw_coef_sta,trim(black)
         end if
      end do
   end select

   if (trim(nam%method)=='dual-ens') then
      ! Compute low-resolution localization diagnostic and fit
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute low-resolution localization diagnostic and fit'

      do ib=1,bpar%nb+1
         if (bpar%diag_block(ib)) then
            ! Allocation
            call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_loc_lr',loc_2(ib))

            ! Compute global low-resolution localization
            write(mpl%unit,'(a7,a)') '','Compute global low-resolution localization'
            call compute_localization(hdata,ib,avg_2(ib),loc_2(ib))

            ! Print results
            do jl0=1,geom%nl0
               write(mpl%unit,'(a10,a,i3,a4,a21,a,f8.2,a)') '','Level: ',nam%levs(jl0),' ~> ', &
             & 'diagonal value: ',trim(peach),loc_2(ib)%raw_coef_ens(jl0),trim(black)
               if (bpar%fit_block(ib)) then
                  write(mpl%unit,'(a45,a,f8.2,a,f8.2,a)') 'loc. support radii: ',trim(aqua),loc_2(ib)%fit_rh(jl0)*reqkm, &
                & trim(black)//' km  / '//trim(aqua),loc_2(ib)%fit_rv(jl0),trim(black)//' '//trim(vunitchar)
               end if
            end do
         end if
      end do

      ! Compute dual-ensemble hybridization diagnostic and fit
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Compute dual-ensemble hybridization diagnostic and fit'

      do ib=1,bpar%nb+1
         if (bpar%diag_block(ib)) then
            ! Allocation
            call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_loc_deh',loc_3(ib))
            call curve_alloc(hdata,ib,trim(bpar%blockname(ib))//'_loc_deh_lr',loc_4(ib))

            ! Compute global dual-ensemble hybridization
            write(mpl%unit,'(a7,a)') '','Compute global dual-ensemble hybridization'
            call compute_dualens(hdata,ib,avg_1(ib),avg_2(ib),loc_3(ib),loc_4(ib))

            ! Print results
            do jl0=1,geom%nl0
               write(mpl%unit,'(a10,a,i3,a4,a21,a,f8.2,a)') '','Level: ',nam%levs(jl0),' ~> ', &
             & 'diagonal value (HR): ',trim(peach),loc_3(ib)%raw_coef_ens(jl0),trim(black)
               write(mpl%unit,'(a45,a,f8.2,a)') 'diagonal value (LR): ',trim(peach),loc_4(ib)%raw_coef_ens(jl0),trim(black)
               if (bpar%fit_block(ib)) then
                  write(mpl%unit,'(a45,a,f8.2,a,f8.2,a)') 'loc. support radii (HR): ',trim(aqua),loc_3(ib)%fit_rh(jl0)*reqkm, &
                & trim(black)//' km  / '//trim(aqua),loc_2(ib)%fit_rv(jl0),trim(black)//' '//trim(vunitchar)
                  write(mpl%unit,'(a45,a,f8.2,a,f8.2,a)') 'loc. support radii (LR): ',trim(aqua),loc_4(ib)%fit_rh(jl0)*reqkm, &
                & trim(black)//' km  / '//trim(aqua),loc_2(ib)%fit_rv(jl0),trim(black)//' '//trim(vunitchar)
               end if
            end do
         end if
      end do
   end if

   if (trim(nam%fit_type)/='none') then
      ! Copy diagnostics into B data
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Copy diagnostics into B data'

      select case (trim(nam%method))
      case ('cor')
         do ib=1,bpar%nb+1
            if (bpar%diag_block(ib)) then
               if (nam%local_diag) then
                  call diag_to_bdata(hdata,ib,cor_1_nc2(:,ib),bdata(ib))
               else
                  call diag_to_bdata(hdata,ib,cor_1(ib),bdata(ib))
               end if
            end if
         end do
      case ('loc')
         do ib=1,bpar%nb+1
            if (bpar%diag_block(ib)) then
               if (nam%local_diag) then
                  call diag_to_bdata(hdata,ib,loc_1_nc2(:,ib),bdata(ib))
               else
                  call diag_to_bdata(hdata,ib,loc_1(ib),bdata(ib))
               end if
            end if
         end do
      case default
         call msgerror('bdata not implemented yet for this method')
      end select

      if (mpl%main) then
         ! Write B data
         write(mpl%unit,'(a)') '-------------------------------------------------------------------'
         write(mpl%unit,'(a,i5,a)') '--- Write B data'

         do ib=1,bpar%nb+1
            if (bpar%diag_block(ib)) call bdata_write(bdata(ib))
         end do
      end if
   end if

   if (mpl%main) then
      ! Write data
      write(mpl%unit,'(a)') '-------------------------------------------------------------------'
      write(mpl%unit,'(a)') '--- Write data'

      ! Displacement
      if (nam%displ_diag) call displ_write(hdata,trim(nam%prefix)//'_displ_diag.nc',displ)

      ! Full variances
      if (nam%full_var) then
         filename = trim(nam%prefix)//'_full_var.nc'
         do ib=1,bpar%nb
            if (bpar%diag_block(ib)) call model_write(nam,geom,filename,trim(bpar%blockname(ib))//'_var', &
          & sum(mom_1(ib)%m2full,dim=3)/float(mom_1(ib)%nsub))
         end do
      end if

      ! Diagnostics
      call curve_write_all(hdata,trim(nam%prefix)//'_diag.nc',cor_1,cor_2,loc_1,loc_2,loc_3,loc_4)

      if (nam%local_diag) then
         ! Fit support radii maps
         if (any(bpar%fit_block)) then
            call curve_write_local(hdata,trim(nam%prefix)//'_local_diag_cor.nc',cor_1_nc2)

            select case (trim(nam%method))
            case ('loc','hyb-avg','hyb-rnd','dual-ens')
               call curve_write_local(hdata,trim(nam%prefix)//'_local_diag_loc.nc',loc_1_nc2)
            end select
         end if

         ! Local diagnostics
         do ildw=1,nam%nldwv
            if (isnotmsi(hdata%nn_ldwv_index(ildw))) then
               write(lonchar,'(f7.2)') nam%lon_ldwv(ildw)
               write(latchar,'(f7.2)') nam%lat_ldwv(ildw)

               call curve_write_all(hdata,trim(nam%prefix)//'_diag_'//trim(adjustl(lonchar))//'-'//trim(adjustl(latchar))//'.nc', &
             & cor_1_nc2(hdata%nn_ldwv_index(ildw),:),cor_2,loc_1_nc2(hdata%nn_ldwv_index(ildw),:),loc_2,loc_3,loc_4)
            else
               call msgwarning('missing local profile')
            end if
         end do
      end if


      ! Write local diagnostic profiles
      do ildw=1,nam%nldwv
         if (isnotmsi(hdata%nn_ldwv_index(ildw))) then
            ! Diagnostic profile
            write(lonchar,'(f7.2)') nam%lon_ldwv(ildw)
            write(latchar,'(f7.2)') nam%lat_ldwv(ildw)
            call curve_write_all(hdata,trim(nam%prefix)//'_diag_'//trim(adjustl(lonchar))//'-'//trim(adjustl(latchar))//'.nc', &
          & cor_1_nc2(:,hdata%nn_ldwv_index(ildw)),cor_2,loc_1_nc2(:,hdata%nn_ldwv_index(ildw)),loc_2,loc_3,loc_4)
         else
            call msgwarning('missing local profile')
         end if
      end do
   end if
elseif (nam%new_param) then
   ! Read B data
   write(mpl%unit,'(a)') '-------------------------------------------------------------------'
   write(mpl%unit,'(a,i5,a)') '--- Read B data'

   do ib=1,bpar%nb+1
      if (bpar%diag_block(ib)) call bdata_read(bdata(ib))
   end do
end if

end subroutine run_hdiag

end module driver_hdiag
