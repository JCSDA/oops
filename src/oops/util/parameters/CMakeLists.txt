# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

oops_target_sources(
  ConfigurationParameter.cc
  ConfigurationParameter.h
  EmptyParameters.h
  GenericParameters.h
  HasDiracParameters_.h
  HasParameters_.h
  HasReadParameters_.h
  HasWriteParameters_.h
  IgnoreOtherParameters.cc
  IgnoreOtherParameters.h
  NumericConstraints.h
  ObjectJsonSchema.cc
  ObjectJsonSchema.h
  OptionalParameter.cc
  OptionalParameter.h
  OptionalPolymorphicParameter.h
  Parameter.h
  ParameterBase.cc
  ParameterBase.h
  ParameterConstraint.h
  ParameterTraits.cc
  ParameterTraits.h
  ParameterTraitsAnyOf.h
  ParameterTraitsScalarOrMap.h
  Parameters.cc
  Parameters.h
  ParametersOrConfiguration.h
  PolymorphicParameter.h
  PolymorphicParameterTraits.h
  PropertyJsonSchema.cc
  PropertyJsonSchema.h
  RequiredParameter.h
  RequiredPolymorphicParameter.h
  )
