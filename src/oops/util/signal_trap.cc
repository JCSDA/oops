/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#ifdef __APPLE__
#include <xmmintrin.h>  // Apple-specific signal handling
#else
#include <cfenv>       // feenableexcept
#endif

#include <execinfo.h>  // backtrace*
#include <cerrno>      // strerror(errno)
#include <csignal>     // sigaction, siginfo_t
#include <cstdlib>     // abort
#include <cstring>     // strerror
#include <iostream>    // cout, cerr
#include "oops/util/Logger.h"  // required for oops::Log

extern void trap_sigfpe(void);                         // user function traps SIGFPE
extern void sigfpe_handler(int, siginfo_t *, void *);  // called when relevant SIGFPE occurs

// #define LOGIT_STDOUT std::cout
// #define LOGIT_STDERR std::cerr
#define LOGIT_STDOUT oops::Log::info()
#define LOGIT_STDERR oops::Log::error()

// This is the user function to enable handling of SIGFPE
void trap_sigfpe(void)
{
  int ret;
  struct sigaction sig_action = {};  // passed to sigaction (init to empty)

#ifdef __APPLE__
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DIV_ZERO);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_OVERFLOW);
#else
  if ((ret = feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)) != 0)
    LOGIT_STDERR << "Call to feenableexcept returned" << ret << std::endl;
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
void sigfpe_handler(int sig, siginfo_t *info, void *ucontext) {
  static const int maxfuncs = 20;  // gather no more than 10 functions in the backtrace
  void *stack[maxfuncs];           // call stack
  char **stacknames;               // backtrace function names
  size_t nfuncs;                   // number of functions returned by backtrace
  size_t n;                        // iterator over nfuncs

  //  myrank = eckit::mpi::comm().rank();
  LOGIT_STDERR << "Caught SIGFPE: " << std::endl;

  switch (info->si_code) {
  case FPE_INTDIV:
    LOGIT_STDERR << "integer divide by zero" << std::endl;
    break;
  case FPE_INTOVF:  // Cannot as yet get this one to trigger
    LOGIT_STDERR << "integer overflow)" << std::endl;
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

  nfuncs = backtrace(stack, maxfuncs);                // number of functions in the backtrace
  stacknames = backtrace_symbols(&stack[0], nfuncs);  // generate the backtrace names from symbols
  // Loop through the callstack, printing the backtrace
  for (n = 0; n < nfuncs; ++n) {
    LOGIT_STDERR << stacknames[n] << std::endl;
  }
  abort();                                            // exit
}
