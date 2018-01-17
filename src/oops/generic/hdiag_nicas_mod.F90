! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic HDIAG_NICAS localization

module hdiag_nicas_mod

use iso_c_binding
use kinds
use config_mod
use unstructured_grid_mod
use driver_hdiag, only: run_hdiag
use driver_nicas, only: run_nicas
use driver_test, only: run_test
use model_oops, only: model_oops_coord
use module_apply_localization, only: apply_localization
use tools_const, only: req
use tools_display, only: listing_setup,msgerror
use tools_missing, only: msi,msr
use type_bdata, only: bdatatype
use type_bpar, only: bpartype,bpar_alloc
use type_geom, only: geomtype,compute_grid_mesh
use type_mpl, only: mpl,mpl_start,mpl_barrier
use type_nam, only: namtype,namcheck
use type_ndata, only: ndataloctype
use type_randgen, only: create_randgen

implicit none
private
public hdiag_nicas, create_hdiag_nicas, delete_hdiag_nicas, hdiag_nicas_multiply

! ------------------------------------------------------------------------------

!>  Derived type containing the data

type hdiag_nicas
  type(namtype) :: nam
  type(geomtype) :: geom
  type(bpartype) :: bpar
  type(bdatatype),allocatable :: bdata(:)
  type(ndataloctype),allocatable :: ndataloc(:)
end type hdiag_nicas

! ------------------------------------------------------------------------------

#define LISTED_TYPE hdiag_nicas

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"
!#include "linkedList.intf.h"

!> Global registry
type(registry_t) :: hdiag_nicas_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "util/linkedList_c.f"
!#include "linkedList.h"

! ------------------------------------------------------------------------------
!  C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_hdiag_nicas_c(key, c_conf, cnv, cnc0a, clats, clons, careas, cnlev, cvunit, cimask, cglbind, cnts, & 
 & cens1_ne, cens1) bind(c, name='create_hdiag_nicas_f90')
implicit none
integer(c_int), intent(inout) :: key
type(c_ptr), intent(in) :: c_conf
integer(c_int), intent(in) :: cnv, cnc0a, cnlev, cens1_ne, cnts
real(c_double), intent(in) :: clats(cnc0a), clons(cnc0a), careas(cnc0a), cvunit(cnlev), cens1(cens1_ne*cnts*cnv*cnlev*cnc0a)
integer(c_int), intent(in) :: cimask(cnlev*cnc0a), cglbind(cnc0a)
type(hdiag_nicas), pointer :: self
integer :: nv,nc0a,nlev,nts,ens1_ne
real(kind=kind_real) :: lats(cnc0a), lons(cnc0a), areas(cnc0a), vunit(cnlev), ens1(cens1_ne*cnts*cnv*cnlev*cnc0a)
integer :: imask(cnlev*cnc0a), glbind(cnc0a)
call hdiag_nicas_registry%init()
call hdiag_nicas_registry%add(key)
call hdiag_nicas_registry%get(key,self)
nv=cnv
lats=clats
lons=clons
areas=careas
vunit=cvunit
imask=cimask
glbind=cglbind
nts=cnts
ens1_ne=cens1_ne
ens1=cens1
call create_hdiag_nicas(self, c_conf, nv, lats, lons, areas, vunit, imask, glbind, nts, ens1_ne, ens1)
end subroutine create_hdiag_nicas_c

! ------------------------------------------------------------------------------

subroutine delete_hdiag_nicas_c(key) bind(c, name='delete_hdiag_nicas_f90')
implicit none
integer(c_int), intent(inout) :: key
type(hdiag_nicas), pointer :: self
call hdiag_nicas_registry%get(key,self)
call delete_hdiag_nicas(self)
call hdiag_nicas_registry%remove(key)
end subroutine delete_hdiag_nicas_c

! ------------------------------------------------------------------------------

subroutine hdiag_nicas_multiply_c(key, idx) bind(c, name='hdiag_nicas_multiply_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx
type(hdiag_nicas), pointer :: self
type(unstructured_grid), pointer :: udx
call hdiag_nicas_registry%get(key,self)
call unstructured_grid_registry%get(idx, udx)
call hdiag_nicas_multiply(self, udx)
end subroutine hdiag_nicas_multiply_c

! ------------------------------------------------------------------------------
!  End C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_hdiag_nicas(self, c_conf, nv, lats, lons, areas, vunit, imask, glbind, nts, ens1_ne, ens1_vec)

implicit none
type(hdiag_nicas), intent(inout) :: self
type(c_ptr), intent(in) :: c_conf
integer, intent(in) :: nv
real(kind=kind_real), intent(in) :: lats(:)
real(kind=kind_real), intent(in) :: lons(:)
real(kind=kind_real), intent(in) :: areas(:)
real(kind=kind_real), intent(in) :: vunit(:)
integer, intent(in) :: imask(:)
integer, intent(in) ::  glbind(:)
integer, intent(in) :: nts
integer, intent(in) :: ens1_ne
real(kind=kind_real), intent(in) :: ens1_vec(:)
integer :: offset,ie,its,ic0a,iv,il0,ib,il
real(kind=kind_real), allocatable :: ens1(:,:,:,:,:)
logical,allocatable :: mask_unpack(:,:)

! Initialize mpl
call mpl_start

! Get sizes 
self%geom%nc0a = size(lats)
self%geom%nlev = size(vunit)

! Read JSON
call hdiag_nicas_read_conf(c_conf,self%nam)

! Force other namelist variables
self%nam%datadir = '.'
self%nam%model = 'oops'
self%nam%colorlog = .false.
self%nam%nv = nv
self%nam%nts = nts
self%nam%ens1_ne = ens1_ne
self%nam%ens1_ne_offset = 0
self%nam%ens1_nsub = 1
do iv=1,self%nam%nv
   write(self%nam%varname(iv),'(a,i2.2)') 'var_',iv
   self%nam%addvar2d(iv) = ''
end do
self%nam%nl = self%geom%nlev
do il=1,self%nam%nl
   self%nam%levs(il) = il
end do

! Setup display
call listing_setup(self%nam%colorlog,self%nam%logpres)

! Check namelist parameters
call namcheck(self%nam)

! Write parallel setup
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a,i3,a,i2,a)') '--- Parallelization with ',mpl%nproc,' MPI tasks and ',mpl%nthread,' OpenMP threads'

! Initialize random number generator
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Initialize random number generator'
call create_randgen(self%nam)

! Initialize coordinates
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Initialize geometry'
call model_oops_coord(self%nam,self%geom,lats,lons,areas,vunit,imask,glbind)

write(mpl%unit,*) 'Number of levels:',self%nam%nl,self%geom%nl0,self%geom%nlev

! Initialize block parameters
call bpar_alloc(self%nam,self%geom,self%bpar)

! Compute grid mesh
write(mpl%unit,'(a)') '-------------------------------------------------------------------'
write(mpl%unit,'(a)') '--- Compute grid mesh'
call compute_grid_mesh(self%nam,self%geom)

! Transform ensemble data from vector to array
allocate(ens1(self%geom%nc0a,self%geom%nl0,self%nam%nv,self%nam%nts,self%nam%ens1_ne))
allocate(mask_unpack(self%geom%nl0,self%nam%nv))
call msr(ens1)
offset = 0
do ie=1,self%nam%ens1_ne
   do its=1,self%nam%nts
      do ic0a=1,self%geom%nc0a
         ens1(ic0a,:,:,its,ie) = unpack(ens1_vec(offset+1:offset+self%geom%nl0*self%nam%nv),mask_unpack,ens1(ic0a,:,:,its,ie))
         offset = offset+self%geom%nl0*self%nam%nv
      end do
   end do
end do

! Call hybrid_diag driver
call run_hdiag(self%nam,self%geom,self%bpar,self%bdata,ens1)

! Call NICAS driver
call run_nicas(self%nam,self%geom,self%bpar,self%bdata,self%ndataloc)

! Call test driver
call run_test(self%nam,self%geom,self%bpar,self%bdata,self%ndataloc)

! Close listing files
call mpl_barrier
write(mpl%unit,*) 'HDIAG_NICAS setup done'
if ((mpl%main.and..not.self%nam%colorlog).or..not.mpl%main) close(unit=mpl%unit)

end subroutine create_hdiag_nicas

!-------------------------------------------------------------------------------

subroutine hdiag_nicas_read_conf(c_conf,nam)
implicit none
type(c_ptr), intent(in) :: c_conf
type(namtype), intent(out) :: nam
integer :: il,its,ildwh,ildwv,idir
character(len=3) :: ilchar,itschar,ildwhchar,ildwvchar,idirchar

! Default initialization

! general_param default
nam%prefix = ''
nam%default_seed = .false.
nam%load_ensemble = .true.

! driver_param default
nam%method = ''
nam%strategy = ''
nam%new_hdiag = .false.
nam%new_param = .false.
nam%new_mpi = .false.
nam%check_adjoints = .false.
nam%check_pos_def = .false.
nam%check_mpi = .false.
nam%check_sqrt = .false.
nam%check_dirac = .false.
nam%check_perf = .false.
nam%check_hdiag = .false.
nam%new_lct = .false.
nam%new_obsop = .false.

! model_param default
nam%logpres = .false.
nam%transform = .false.

! sampling_param default
nam%sam_write = .true.
nam%sam_read = .false.
nam%mask_type = ''
call msr(nam%mask_th)
nam%mask_check = .false.
call msi(nam%nc1)
call msi(nam%ntry)
call msi(nam%nrep)
call msi(nam%nc)
call msr(nam%dc)
call msi(nam%nl0r)

! diag_param default
call msi(nam%ne)
nam%gau_approx = .false.
nam%local_diag = .false.
call msr(nam%local_rad)
nam%displ_diag = .false.
call msr(nam%displ_rad)
call msi(nam%displ_niter)
call msr(nam%displ_rhflt)
call msr(nam%displ_tol)

! fit_param default
nam%fit_type = ''
nam%fit_wgt = .false.
nam%lhomh = .false.
nam%lhomv = .false.
call msr(nam%rvflt)
nam%lct_diag = .false.

! output_param default
call msi(nam%nldwh)
call msi(nam%il_ldwh)
call msi(nam%ic_ldwh)
call msi(nam%nldwv)
call msr(nam%lon_ldwv)
call msr(nam%lat_ldwv)
nam%flt_type = ''
call msr(nam%diag_rhflt)

! nicas_param default
nam%lsqrt = .false.
call msr(nam%resol)
nam%network = .false.
call msi(nam%mpicom)
call msi(nam%ndir)
call msr(nam%londir)
call msr(nam%latdir)
call msi(nam%levdir)
call msi(nam%ivdir)
call msi(nam%itsdir)

! Setup from configuration

! general_param
nam%prefix = config_get_string(c_conf,1024,"prefix")
nam%default_seed = integer_to_logical(config_get_int(c_conf,"default_seed"))

! driver_param default
nam%method = config_get_string(c_conf,1024,"method")
nam%strategy = config_get_string(c_conf,1024,"strategy")
nam%new_hdiag = integer_to_logical(config_get_int(c_conf,"new_hdiag"))
nam%new_param = integer_to_logical(config_get_int(c_conf,"new_param"))
nam%new_mpi = integer_to_logical(config_get_int(c_conf,"new_mpi"))
nam%check_adjoints = integer_to_logical(config_get_int(c_conf,"check_adjoints"))
nam%check_pos_def = integer_to_logical(config_get_int(c_conf,"check_pos_def"))
nam%check_mpi = integer_to_logical(config_get_int(c_conf,"check_mpi"))
nam%check_sqrt = integer_to_logical(config_get_int(c_conf,"check_sqrt"))
nam%check_dirac = integer_to_logical(config_get_int(c_conf,"check_dirac"))
nam%check_hdiag = integer_to_logical(config_get_int(c_conf,"check_hdiag"))
nam%new_lct = integer_to_logical(config_get_int(c_conf,"new_lct"))
nam%new_obsop = integer_to_logical(config_get_int(c_conf,"new_obsop"))

! model_param default
nam%logpres = integer_to_logical(config_get_int(c_conf,"logpres"))
nam%transform = integer_to_logical(config_get_int(c_conf,"transform"))

! sampling_param default
nam%mask_type = config_get_string(c_conf,1024,"mask_type")
nam%mask_th = config_get_real(c_conf,"mask_th")
nam%mask_check = integer_to_logical(config_get_int(c_conf,"mask_check"))
nam%nc1 = config_get_int(c_conf,"nc1")
nam%ntry = config_get_int(c_conf,"ntry")
nam%nrep = config_get_int(c_conf,"nrep")
nam%nc = config_get_int(c_conf,"nc")
nam%dc = config_get_real(c_conf,"dc")/req
nam%nl0r = config_get_int(c_conf,"nl0r")

! diag_param default
nam%ne = config_get_int(c_conf,"ne")
nam%gau_approx = integer_to_logical(config_get_int(c_conf,"gau_approx"))
nam%full_var = integer_to_logical(config_get_int(c_conf,"full_var"))
nam%local_diag = integer_to_logical(config_get_int(c_conf,"local_diag"))
nam%local_rad = config_get_real(c_conf,"local_rad")/req
nam%displ_diag = integer_to_logical(config_get_int(c_conf,"displ_diag"))
nam%displ_rad = config_get_real(c_conf,"displ_rad")/req
nam%displ_niter = config_get_int(c_conf,"displ_niter")
nam%displ_rhflt = config_get_real(c_conf,"displ_rhflt")/req
nam%displ_tol = config_get_real(c_conf,"displ_tol")

! fit_param default
nam%fit_type = config_get_string(c_conf,1024,"fit_type")
nam%fit_wgt = integer_to_logical(config_get_int(c_conf,"fit_wgt"))
nam%lhomh = integer_to_logical(config_get_int(c_conf,"lhomh"))
nam%lhomv = integer_to_logical(config_get_int(c_conf,"lhomv"))
nam%rvflt = config_get_real(c_conf,"rvflt")
nam%lct_diag = integer_to_logical(config_get_int(c_conf,"lct_diag"))
nam%lct_nscales = config_get_int(c_conf,"lct_nscales")

! output_param default
nam%nldwh = config_get_int(c_conf,"nldwh")
do ildwh=1,nam%nldwh
   write(ildwvchar,'(i3)') ildwh
   nam%il_ldwh(ildwh) = config_get_int(c_conf,"il_ldwh("//trim(adjustl(ildwhchar))//")")
   nam%ic_ldwh(ildwh) = config_get_int(c_conf,"ic_ldwh("//trim(adjustl(ildwhchar))//")")
end do
nam%nldwv = config_get_int(c_conf,"nldwv")
do ildwv=1,nam%nldwv
   write(ildwvchar,'(i3)') ildwv
   nam%lon_ldwv(ildwv) = config_get_real(c_conf,"lon_ldwv("//trim(adjustl(ildwvchar))//")")
   nam%lat_ldwv(ildwv) = config_get_real(c_conf,"lat_ldwv("//trim(adjustl(ildwvchar))//")")
end do
nam%flt_type = config_get_string(c_conf,1024,"flt_type")
nam%diag_rhflt = config_get_real(c_conf,"diag_rhflt")/req

! nicas_param default
nam%lsqrt = integer_to_logical(config_get_int(c_conf,"lsqrt"))
nam%resol = config_get_real(c_conf,"resol")
nam%network = integer_to_logical(config_get_int(c_conf,"network"))
nam%mpicom = config_get_int(c_conf,"mpicom")
nam%ndir = config_get_int(c_conf,"ndir")
do idir=1,nam%ndir
   write(idirchar,'(i3)') idir
   nam%londir(idir) = config_get_real(c_conf,"londir("//trim(adjustl(idirchar))//")")
   nam%latdir(idir) = config_get_real(c_conf,"latdir("//trim(adjustl(idirchar))//")")
   nam%levdir(idir) = config_get_int(c_conf,"levdir("//trim(adjustl(idirchar))//")")
   nam%ivdir(idir) = config_get_int(c_conf,"ivdir("//trim(adjustl(idirchar))//")")
   nam%itsdir(idir) = config_get_int(c_conf,"itsdir("//trim(adjustl(idirchar))//")")
end do

end subroutine hdiag_nicas_read_conf

!-------------------------------------------------------------------------------

logical function integer_to_logical(i)
implicit none
integer,intent(in) :: i

if (i==0) then
   integer_to_logical = .false.
elseif (i==1) then
   integer_to_logical = .true.
else
   call abor1_ftn('wrong integer in integer_to_logical')
end if

end function integer_to_logical

!-------------------------------------------------------------------------------

subroutine delete_hdiag_nicas(self)
implicit none
type(hdiag_nicas), intent(inout) :: self

end subroutine delete_hdiag_nicas

!-------------------------------------------------------------------------------

subroutine hdiag_nicas_multiply(self,dx)
implicit none
type(hdiag_nicas), intent(inout) :: self
type(unstructured_grid), intent(inout) :: dx
integer :: iv,ic0a,il0
real(kind_real) :: fld(self%geom%nc0a,self%geom%nl0,self%nam%nv,self%nam%nts)
type(column_element), pointer :: current

! Multiply with HDIAG_NICAS
write(mpl%unit,*) "HDIAG_NICAS multiply"

! Initialization
ic0a = 1
current => dx%head

do while (associated(current))
   ! Copy field
   fld(ic0a,:,:,1) = current%column%fld

   ! Update index and pointer
   ic0a = ic0a+1
   current => current%next
end do

! Apply localization
call apply_localization(self%nam,self%geom,self%bpar,self%ndataloc,fld)

! Initialization
ic0a = 1
current => dx%head

do while (associated(current))
   ! Copy field
   current%column%fld = fld(ic0a,:,:,1)

   ! Update index and pointer
   ic0a = ic0a+1
   current => current%next
end do

end subroutine hdiag_nicas_multiply

!-------------------------------------------------------------------------------

end module hdiag_nicas_mod
