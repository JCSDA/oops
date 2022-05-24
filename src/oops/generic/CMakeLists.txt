# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

oops_target_sources(
  AnalyticInitBase.h
  GeoVaLsWriter.h
  IdentityLinearModel.h
  IdentityModel.h
  LinearModelBase.h
  LocalizationBase.h
  ModelBase.h
  ObsErrorBase.h
  ObsErrorDiag.cc
  ObsErrorDiag.h
  ObsFilterBase.h
  ObsFilterParametersBase.h
  PseudoModel.h
  PseudoModelState4D.h
  UnstructuredInterpolator.h
  VerticalLocEV.h
  fft_gpoint2spectral_f.F90
  fft_init_f.F90
  fft_interface_f.h
  fft_multiple_f.F90
  fft_pass_proc_f.F90
  fft_spectral2gpoint_f.F90
  gc99.cc
  gc99.h
  gc99_mod.F90
  instantiateLinearModelFactory.h
  instantiateModelFactory.h
  instantiateObsErrorFactory.h
  soar.cc
  soar.h
  )
