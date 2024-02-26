/*
 * (C) Copyright 2019-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#include "oops/util/signal_trap.h"

#include <cerrno>      // strerror(errno)
#include <cfenv>
#include <csignal>     // sigaction, siginfo_t
#include <cstring>     // strerror
#include <iostream>    // cout, cerr
#include <string>

// TODO(srherbener) 10/19/23: It was decided for now to only trap for invalid and
// divide by zero, and to remove trapping for overflow due to the many variations of
// handling floating point data on both hardware (eg, extended precision) and
// compilers (eg, floating point models).
//
// The invalid (eg, sqrt(-1)) and divide by zero exceptions are clearly defined and behave
// consistently across the variations mentioned above.
//
// The only floating point exception that was tested up to this point was divide by zero,
// so we are still adding a test for invalid operations.
//
// The trapping of overflow can be revisited in the future.

// Note that this routine is using the iostream channels (cout, cerr) directly instead
// of the oops logging channels (eg, oops::Log::info()). This was done so that the ctest
// could run test code directly (instead of as an OOPS application), create floating point
// exceptions that cause an abort, and then use the CMake WILL_FAIL test property to catch
// the abort.

#include "oops/util/abor1_cpp.h"
#include "oops/util/Stacktrace.h"

namespace util {

static bool do_abort;

//------------------------------------------------------------------------
// private functions
//------------------------------------------------------------------------
static void enable_floating_point_exceptions(const int abortFlag);
static void sigfpe_abort();

#if defined(__APPLE__)
#if defined(__aarch64__) || defined(__arm64__)

//------------------------------------------------------------------------
// macOS, arm architecture
//------------------------------------------------------------------------

static void fpe_signal_handler(const int sig);

const char trapSigOsArchName[] = "macOS, arm64";

//------------------------------------------------------------------------
void enable_floating_point_exceptions(const int abortFlag) {
    // Need a global variable to pass in the abort flag
    do_abort = (abortFlag > 0);

    // Set config to trap on divide by zero, invalid operation and overflow.
    // Note that the arm64 architecture sends a SIGILL signal (not a SIGFPE).
    // The fenv_t struct provides two data members:
    //  __fpcr floating point control register
    //  __fpsr floating point status register
    std::fenv_t env;
    std::fegetenv(&env);
    env.__fpcr = env.__fpcr |
        (__fpcr_trap_invalid | __fpcr_trap_divbyzero);
    std::fesetenv(&env);

    std::signal(SIGILL, fpe_signal_handler);
}

//------------------------------------------------------------------------
void fpe_signal_handler(const int sig)
{
    // Use the fenv_t struct to query the results
    fenv_t test_env;
    fegetenv(&test_env);
    std::cerr << std::endl << "Floating point exception detected: ";
    if (std::fetestexcept(FE_DIVBYZERO)) {
        std::cerr << "Divide by zero";
    } else if (std::fetestexcept(FE_INVALID)) {
        std::cerr << "Invalid operation";
    } else {
        std::cerr << "Arithmetic exception";
    }
    std::cerr << std::endl;

    if (do_abort) {
        sigfpe_abort();
    }
}

#else
//------------------------------------------------------------------------
// macOS, intel architecture
//------------------------------------------------------------------------
#include <xmmintrin.h>

static void fpe_signal_action(const int sig, siginfo_t *sip, void *scp);

const char trapSigOsArchName[] = "macOS, Intel";

//------------------------------------------------------------------------
void enable_floating_point_exceptions(const int abortFlag) {
    // Need a global variable to pass in the abort flag
    do_abort = (abortFlag > 0);

    // Set the main floating point control register using the fenv style. Note
    // that 0 means enable, 1 means disable.
    // The fenv_t struct provides three data members:
    //  __fpcr floating point control register
    //  __fpsr floating point status register
    //  __mxcsr SSE floating point control/status register (see macros below)
    fenv_t env;
    fegetenv(&env);
    env.__control = env.__control & ~(FE_INVALID | FE_DIVBYZERO);
    fesetenv(&env);

    // Set the mxcsr register (fpe for SSE instructions) directly via the provided
    // macros. Note that a value of 0 means enable, and a value of 1 means disable
    // In this case, the intel architecture issues a SIGFPE when a floating
    // point exception occurs
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DIV_ZERO);
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    // Need to utilize the sigaction style of exception handling
    struct sigaction act;
    act.sa_sigaction = fpe_signal_action;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_SIGINFO;
    if (sigaction(SIGFPE, &act, NULL) != 0) {
        std::cerr << "Call to sigaction failed so trapping of SIGFPE disabled:"
                  << strerror(errno) << std::endl;
    }
}

//------------------------------------------------------------------------
void fpe_signal_action(int sig, siginfo_t *sip, void *scp) {
    /* see signal.h for codes */
    if (sig == SIGFPE) {
        std::cerr << std::endl << "Floating point exception detected: ";
    }

    switch (sip->si_code) {
        case FPE_FLTINV:
            std::cerr << "Invalid operation";
            break;
        case FPE_FLTDIV:
            std::cerr << "Divide by zero";
            break;
        default:
            std::cerr << "Arithmetic exception";
    }
    std::cerr << std::endl;

    if (do_abort) {
        sigfpe_abort();
    }
}

#endif
#else
//------------------------------------------------------------------------
// Linux
//------------------------------------------------------------------------

static void fpe_signal_action(int sig, siginfo_t *info, void *ucontext);

const char trapSigOsArchName[] = "Linux";

//------------------------------------------------------------------------
void enable_floating_point_exceptions(const int abortFlag) {
    // Need a global variable to pass in the abort flag
    do_abort = (abortFlag > 0);

    int ret;
    if ((ret = feenableexcept (FE_DIVBYZERO | FE_INVALID)) != 0)
      std::cerr << "Call to feenableexcept returned" << ret << std::endl;

    struct sigaction sig_action = {};             // passed to sigaction (init to empty)
    sig_action.sa_flags = SA_SIGINFO;             // handler specified in sa_sigaction
    sig_action.sa_sigaction = fpe_signal_action;  // function name
    sigemptyset(&sig_action.sa_mask);             // initialize mask
    sigaddset(&sig_action.sa_mask, SIGFPE);       // disable SIGFPE while another is
                                                  // being processed

    if (sigaction(SIGFPE, &sig_action, NULL) != 0) {
        std::cerr << "Call to sigaction failed so trapping of SIGFPE disabled:"
                  << strerror(errno) << std::endl;
  }
}

//------------------------------------------------------------------------
// This is the signal handler invoked when SIGFPE encountered
// Arguments "sig" and "ucontext" are required but unused here
void fpe_signal_action(int sig, siginfo_t *info, void *ucontext) {
    /* see signal.h for codes */
    if (sig == SIGFPE) {
        std::cerr << std::endl << "Floating point exception detected: ";
    }

    switch (info->si_code) {
    case FPE_FLTDIV:
      std::cerr << "Divide by zero";
      break;
    case FPE_FLTINV:
      std::cerr << "Invalid operation";
      break;
    default:
      std::cerr << "Arithmetic exception";
      break;
    }
    std::cerr << std::endl;

    if (do_abort) {
        sigfpe_abort();
    }
}
#endif

//------------------------------------------------------------------------
// Shared functions
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void trap_sigfpe(const int abortFlag) {
    // Set up trap for:
    //    divide by zero
    //    invalid operation
    //    overflow
    std::cout << "INFO: util::trap_sigfpe: Enabling FPE signal trapping for "
              << trapSigOsArchName << std::endl;

    // call the OS/architecture specific code for enabling trapping and handling
    // of floating point exception signals.
    enable_floating_point_exceptions(abortFlag);
}

//------------------------------------------------------------------------
void sigfpe_abort() {
    std::cerr << std::endl << stacktrace_current() << std::endl;
    ABORT("Trapped a floating point exception");
}

}   // namespace util
