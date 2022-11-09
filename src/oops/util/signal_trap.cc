/*
 * (C) Copyright 2019-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#if defined  (__APPLE__) && defined (__aarch64__)
void trap_sigfpe(const int abortflg)
{
    // do nothing
}
#else

#ifdef __APPLE__
#include <xmmintrin.h>  // Apple-specific signal handling
#else
#include <cfenv>       // feenableexcept
#endif

#include <cerrno>      // strerror(errno)
#include <csignal>     // sigaction, siginfo_t
#include <cstring>     // strerror
#include <iostream>    // cout, cerr

// The behavior of boost::stacktrace is controlled by preprocessor macros. See
// https://www.boost.org/doc/libs/1_76_0/doc/html/stacktrace/configuration_and_build.html
// for details. The correct libs and macros are set using CMake (in backtrace_deps.cmake).
#include <boost/stacktrace.hpp>  // boost stacktraces

#include "oops/util/abor1_cpp.h"  // ABORT macro
#include "oops/util/Logger.h"     // required for oops::Log

void trap_sigfpe(const int);                    // user function traps SIGFPE
void sigfpe_handler(int, siginfo_t *, void *);  // called when relevant SIGFPE occurs

// Global bool required here because sa.sigaction function does not provide it
static int do_abort;  // Whether to abort on SIGFPE

// #define LOGIT_STDOUT std::cout
// #define LOGIT_STDERR std::cerr
#define LOGIT_STDOUT oops::Log::info()
#define LOGIT_STDERR oops::Log::error()

// This is the user function to enable handling of SIGFPE
void trap_sigfpe(const int abortflg)
{
  struct sigaction sig_action = {};  // passed to sigaction (init to empty)

  do_abort = abortflg;
#ifdef __APPLE__
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DIV_ZERO);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_OVERFLOW);
#else
  int ret;
  if ((ret = feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)) != 0)
    LOGIT_STDERR << "Call to feenableexcept returned" << ret << std::endl;
#endif

#ifdef __PGI
  ret = fedisableexcept(FE_OVERFLOW);
#endif

  sig_action.sa_flags = SA_SIGINFO;          // handler specified in sa_sigaction
  sig_action.sa_sigaction = sigfpe_handler;  // function name
  sigemptyset(&sig_action.sa_mask);          // initialize mask
  sigaddset(&sig_action.sa_mask, SIGFPE);    // disable SIGFPE while another is being processed

  if (sigaction(SIGFPE, &sig_action, NULL) == 0)
    LOGIT_STDOUT << "Trapping of SIGFPE enabled" << std::endl;
  else
    LOGIT_STDERR << "Call to sigaction failed so trapping of SIGFPE disabled:" << strerror(errno)
                 << std::endl;
}

// This is the signal handler invoked when SIGFPE encountered
// Arguments "sig" and "ucontext" are required but unused here
void sigfpe_handler(int sig, siginfo_t *info, void *ucontext) {
  //  myrank = eckit::mpi::comm().rank();
  LOGIT_STDERR << "Caught SIGFPE: ";

  switch (info->si_code) {
  case FPE_INTDIV:
    LOGIT_STDERR << "integer divide by zero" << std::endl;
    break;
  case FPE_INTOVF:  // Cannot as yet get this one to trigger
    LOGIT_STDERR << "integer overflow" << std::endl;
    break;
  case FPE_FLTDIV:
    LOGIT_STDERR << "floating-point divide by zero" << std::endl;
    break;
  case FPE_FLTOVF:
    LOGIT_STDERR << "floating-point overflow" << std::endl;
    break;
  case FPE_FLTINV:
    LOGIT_STDERR << "floating-point invalid operation" << std::endl;
    break;
  case FPE_FLTSUB:  // Cannot as yet get this one to trigger
    LOGIT_STDERR << "subscript out of range" << std::endl;
    break;
  default:
    LOGIT_STDERR << "Arithmetic Exception" << std::endl;
    break;
  }

  if (do_abort) {
    LOGIT_STDERR << boost::stacktrace::stacktrace() << std::endl;
    ABORT("A SIGFPE was encountered. If available, see frame #2 in the stacktrace for details.");
  }
}

#endif
