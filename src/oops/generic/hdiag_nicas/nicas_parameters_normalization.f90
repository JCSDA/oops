!----------------------------------------------------------------------
! Module: nicas_parameters_normalization.f90
!> Purpose: normalization routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_parameters_normalization

use omp_lib
use tools_display, only: msgerror,prog_init,prog_print
use tools_kinds,only: kind_real
use tools_missing, only: msr,isnotmsi,msi
use tools_qsort, only: qsort
use type_geom, only: fld_com_lg
use type_mpl, only: mpl
use type_ndata, only: ndatatype

implicit none

private
public :: compute_normalization

contains

!----------------------------------------------------------------------
! Subroutine: compute_normalization
!> Purpose: compute normalization
!----------------------------------------------------------------------
subroutine compute_normalization(ndata)

implicit none

! Passed variables
type(ndatatype),intent(inout) :: ndata !< NICAS data

! Local variables
integer :: il0i,i_s,ic1,ic1b,jc1b,is,js,isc,jsb,jsc,ic0,ic0a,il0,il1,ih,jv,nlr,ilr,jlr,ic,isc_add,progint
integer,allocatable :: ineh(:,:),inev(:),ines(:,:),inec(:),order(:),isc_list(:)
integer,allocatable :: h_col(:,:,:),v_col(:,:),s_col(:,:,:),c_ind(:,:)
real(kind_real) :: S_add
real(kind_real),allocatable :: h_S(:,:,:),v_S(:,:),s_S(:,:,:),c_S(:,:)
real(kind_real),allocatable :: S_list(:),S_list_tmp(:)
logical :: conv
logical,allocatable :: done(:)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

! Compute horizontal interpolation inverse mapping
allocate(ineh(geom%nc0a,geom%nl0i))
ineh = 0
do il0i=1,geom%nl0i
   do i_s=1,ndata%h(il0i)%n_s
      ic0a = ndata%h(il0i)%row(i_s)
      ineh(ic0a,il0i) = ineh(ic0a,il0i)+1
   end do
end do
allocate(h_col(maxval(ineh),geom%nc0a,geom%nl0i))
allocate(h_S(maxval(ineh),geom%nc0a,geom%nl0i))
ineh = 0
do il0i=1,geom%nl0i
   do i_s=1,ndata%h(il0i)%n_s
      ic0a = ndata%h(il0i)%row(i_s)
      ineh(ic0a,il0i) = ineh(ic0a,il0i)+1
      h_col(ineh(ic0a,il0i),ic0a,il0i) = ndata%h(il0i)%col(i_s)
      h_S(ineh(ic0a,il0i),ic0a,il0i) = ndata%h(il0i)%S(i_s)
   end do
end do

! Compute vertical interpolation inverse mapping
allocate(inev(geom%nl0))
inev = 0
do i_s=1,ndata%v%n_s
   il0 = ndata%v%row(i_s)
   inev(il0) = inev(il0)+1
end do
allocate(v_col(maxval(inev),geom%nl0))
allocate(v_S(maxval(inev),geom%nl0))
inev = 0
do i_s=1,ndata%v%n_s
   il0 = ndata%v%row(i_s)
   inev(il0) = inev(il0)+1
   v_col(inev(il0),il0) = ndata%v%col(i_s)
   v_S(inev(il0),il0) = ndata%v%S(i_s)
end do

! Compute subsampling interpolation inverse mapping
allocate(ines(ndata%nc1b,ndata%nl1))
ines = 0
do il1=1,ndata%nl1
   do i_s=1,ndata%s(il1)%n_s
      ic1b = ndata%s(il1)%row(i_s)
      ines(ic1b,il1) = ines(ic1b,il1)+1
   end do
end do
allocate(s_col(maxval(ines),ndata%nc1b,ndata%nl1))
allocate(s_S(maxval(ines),ndata%nc1b,ndata%nl1))
ines = 0
do il1=1,ndata%nl1
   do i_s=1,ndata%s(il1)%n_s
      ic1b = ndata%s(il1)%row(i_s)
      ines(ic1b,il1) = ines(ic1b,il1)+1
      jc1b = ndata%s(il1)%col(i_s)
      jsb = ndata%c1bl1_to_sb(jc1b,il1)
      jsc = ndata%sb_to_sc_nor(jsb)
      s_col(ines(ic1b,il1),ic1b,il1) = jsc
      s_S(ines(ic1b,il1),ic1b,il1) = ndata%s(il1)%S(i_s)
   end do
end do

! Compute convolution inverse mapping
allocate(inec(ndata%nsc_nor))
inec = 0
do i_s=1,ndata%c_nor%n_s
   isc = ndata%c_nor%col(i_s)
   jsc = ndata%c_nor%row(i_s)
   is = ndata%sc_nor_to_s(isc)
   js = ndata%sc_nor_to_s(jsc)
   inec(isc) = inec(isc)+1
   inec(jsc) = inec(jsc)+1
end do
allocate(c_ind(maxval(inec),ndata%nsc_nor))
allocate(c_S(maxval(inec),ndata%nsc_nor))
call msi(c_ind)
call msr(c_S)
inec = 0
do i_s=1,ndata%c_nor%n_s
   isc = ndata%c_nor%col(i_s)
   jsc = ndata%c_nor%row(i_s)
   is = ndata%sc_nor_to_s(isc)
   js = ndata%sc_nor_to_s(jsc)
   inec(isc) = inec(isc)+1
   c_ind(inec(isc),isc) = jsc
   c_S(inec(isc),isc) = ndata%c_nor%S(i_s)
   inec(jsc) = inec(jsc)+1
   c_ind(inec(jsc),jsc) = isc
   c_S(inec(jsc),jsc) = ndata%c_nor%S(i_s)
end do

! Re-order indices
do isc=1,ndata%nsc_nor
   allocate(order(inec(isc)))
   call qsort(inec(isc),c_ind(1:inec(isc),isc),order)
   c_S(1:inec(isc),isc) = c_S(order(1:inec(isc)),isc)
   deallocate(order)
end do

! Allocation
allocate(ndata%norm(geom%nc0a,geom%nl0))
allocate(done(geom%nc0a))
call msr(ndata%norm)

! Compute normalization weights
do il0=1,geom%nl0
   il0i = min(il0,geom%nl0i)
   write(mpl%unit,'(a10,a,i3,a)',advance='no') '','Level ',nam%levs(il0),': '
   call prog_init(progint,done)

   !$omp parallel do schedule(static) private(ic0a,ic0,isc_list,order,S_list,S_list_tmp,nlr) &
   !$omp&                             private(isc_add,S_add,ih,ic1b,ic1,jv,il1,is,ilr,conv,ic,jlr,isc,jsc)
   do ic0a=1,geom%nc0a
      ! Index
      ic0 = geom%c0a_to_c0(ic0a)

      if (geom%mask(ic0,il0)) then
         ! Allocation
         allocate(isc_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         allocate(order(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         allocate(S_list(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         if (nam%lsqrt) then
            allocate(S_list_tmp(ndata%nsc_nor))
         else
            allocate(S_list_tmp(ineh(ic0a,il0i)*inev(il0)*maxval(ines)))
         end if

         ! Initialization
         isc_list = 0
         S_list = 0.0

         ! Adjoint interpolation
         nlr = 0
         do ih=1,ineh(ic0a,il0i)
            ic1b = h_col(ih,ic0a,il0i)
            ic1 = ndata%c1b_to_c1(ic1b)
            do jv=1,inev(il0)
               il1 = v_col(jv,il0)
               if ((ndata%l1_to_l0(il1)>=ndata%vbot(ic1)).and.(ndata%l1_to_l0(il1)<=ndata%vtop(ic1))) then
                  do is=1,ines(ic1b,il1)
                     isc_add = s_col(is,ic1b,il1)
                     S_add = h_S(ih,ic0a,il0i)*v_S(jv,il0)*s_S(is,ic1b,il1)
                     if (nlr==0) then
                        ilr = 1
                        nlr = 1
                     else
                        do ilr=1,nlr
                           if (isc_add==isc_list(ilr)) exit
                        end do
                        if (ilr==nlr+1) nlr = nlr+1
                     end if
                     isc_list(ilr) = isc_add
                     S_list(ilr) = S_list(ilr)+S_add
                  end do
               end if
            end do
         end do

         if (nam%lsqrt) then
            ! Initialization
            S_list_tmp = 0.0
            do ilr=1,nlr
               isc = isc_list(ilr)
               S_list_tmp(isc) = S_list(ilr)
            end do

            ! Convolution
            do ilr=1,nlr
               isc = isc_list(ilr)
               do ic=1,inec(isc)
                  jsc = c_ind(ic,isc)
                  S_list_tmp(jsc) = S_list_tmp(jsc)+c_S(ic,isc)*S_list(ilr)
               end do
            end do

            ! Sum of squared values
            ndata%norm(ic0a,il0) = sum(S_list_tmp**2)

            ! Normalization factor
            ndata%norm(ic0a,il0) = 1.0/sqrt(ndata%norm(ic0a,il0))
         else
            ! Sort arrays
            call qsort(nlr,isc_list(1:nlr),order(1:nlr))
            S_list(1:nlr) = S_list(order(1:nlr))

            ! Convolution
            S_list_tmp(1:nlr) = S_list(1:nlr)
            do ilr=1,nlr
               isc = isc_list(ilr)
               conv = .true.
               ic = 1
               do jlr=ilr+1,nlr
                  jsc = isc_list(jlr)
                  if (ic<=inec(isc)) then
                     do while (c_ind(ic,isc)/=jsc)
                        ic = ic+1
                        if (ic>inec(isc)) then
                           conv = .false.
                           exit
                        end if
                     end do
                  else
                     conv = .false.
                  end if
                  if (conv) then
                     S_list(ilr) = S_list(ilr)+c_S(ic,isc)*S_list_tmp(jlr)
                     S_list(jlr) = S_list(jlr)+c_S(ic,isc)*S_list_tmp(ilr)
                  end if
               end do
            end do

            ! Normalization factor
            if (nlr>0) ndata%norm(ic0a,il0) = 1.0/sqrt(sum(S_list(1:nlr)*S_list_tmp(1:nlr)))
         end if

         ! Print progression
         done(ic0a) = .true.
         call prog_print(progint,done)

         ! Release memory
         deallocate(isc_list)
         deallocate(order)
         deallocate(S_list)
         deallocate(S_list_tmp)
      end if
   end do
   !$omp end parallel do
   write(mpl%unit,'(a)') '100%'
end do

! Release memory
deallocate(done)

! End associate
end associate

end subroutine compute_normalization

end module nicas_parameters_normalization
