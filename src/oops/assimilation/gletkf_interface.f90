!
! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gletkf_interface

  use iso_c_binding
  use letkf
  implicit none

  private

contains
  !------------------------------------------------------------------------------
  !> 
  subroutine letkf_core_c(nobsl, chxens, chxens_orig, cdep,&
                          cwts_ensmean, cwts_ensperts,&
                          crdiaginv_loc, nanals, neigv,&
                          getkf_inflation, denkf, getkf, mult_infl) &
                          bind(c, name='letkf_core_f90')
    ! Passed variables
    integer(c_int) :: nobsl, nanals, neigv
    integer(c_int) :: getkf_inflation, denkf, getkf
    real(c_float) :: chxens(nanals,nobsl), chxens_orig(nanals,nobsl), &
                      cdep(nobsl), cwts_ensmean(nanals), &
                      cwts_ensperts(nanals,nanals), crdiaginv_loc(nobsl)
    real(c_float) :: mult_infl

    ! getkf_inflation, denkf, getkf
    ! are passed as integer but cast to logical here
    call letkf_core(nobsl,chxens,chxens_orig,cdep,&
                      cwts_ensmean,cwts_ensperts,&
                      crdiaginv_loc,nanals,neigv,&
                      getkf_inflation==1,denkf==1,getkf==1, mult_infl)

  end subroutine letkf_core_c
end module gletkf_interface 

