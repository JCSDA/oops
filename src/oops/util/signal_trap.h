#pragma once
/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

namespace util {

/// \brief set up floating point signal trapping
/// \details This routine will set up trapping two floating point exceptions: divide
/// by zero and invalid operation (eg, sqrt(-1.0)). Note that the apple intel emulator
/// Rosetta disallows floating point exception trapping so this is not supported.
/// However, macOS native arm64 and native intel architectures are supported.
/// Note the abortFlag is an int instead of a bool since it can be set from the environment
/// variable OOPS_ABORTFPE (which is set to 1 by default).
/// \param abortFlag when == 1, abort the process when a floating point signal is caught
void trap_sigfpe(const int abortFlag);

}  // namespace util
