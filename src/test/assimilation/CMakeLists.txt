# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

oops_target_sources(
  # Many of the .cc files in this directory are test executables, and are compiled
  # in the ecbuild_add_test commands in the higher-level CMakeLists
  FullGMRES.h
  SolveMatrixEquation.h
  SpectralLMP.h
  TestVector3D.h
  TriDiagSolve.h
  Vector3D.cc
  Vector3D.h
  rotmat.h
  )
