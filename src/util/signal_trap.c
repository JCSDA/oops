/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifdef _AIX43

#include <signal.h>
#include <fptrap.h>
#include <memory.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>

#define ALARM_CALL 2*60

/*
  ------------------------------------------------------------------------
  'signal_trap_dump' is a GLOBAL variable that is a copy of the value
  of the 'create_dump' parameter passed into the 'signal_trap' subroutine.
  ------------------------------------------------------------------------
*/

int signal_trap_dump = 0;

/*
  --------------------------------------------------------------------
  'signal_trap_incr' is a special variable that is used to ensure that
  if more than one thread of the process gets a signal more or less
  concurrently, only one of them produces a traceback. This makes it 
  easier to see what happened. This variable is of type 
  'static volatile sig_atomic_t (defined in <sys/signal.h>.
  The compiler gives a warning message that says:

                   Duplicate type qualifier "volatile" ignored.

  This warning message should be ignored. Do not be temptes to remove
  the 'volatile' from the definition.
  --------------------------------------------------------------------
*/

static volatile sig_atomic_t signal_trap_incr = 0;
/* Sami's version static sig_atomic_t signal_trap_incr = 0;*/

extern void signal_trace(int, siginfo_t *, void *);

/*======================================================*/
int signal_trap_(int32_t *create_dump, int32_t *signals)
/*======================================================*/
{
/*
   ----------------------------------------------------------
   The entry-point used if '-qextname' compiler flag is used.
   ----------------------------------------------------------
*/
   return signal_trap(create_dump, signals);
}

/*=====================================================*/
int signal_trap(int32_t *create_dump, int32_t *signals)
/*=====================================================*/
{
/************************************************************************
 NAME
       signal_trap - Signal handling with traceback for Fortran programs.

 SYNOPSIS
       int signal_trap (int32_t *create_dump, int32_t *signals)
       int signal_trap_(int32_t *create_dump, int32_t *signals) [for -qextname]

 DESCRIPTION
       The signal_trap subroutine installs a signal handler for various
       signals that could cause a Fortran program to abort. The handler
       prints out a traceback and optionally creates a core dump.

       "sigaction" is used rather than "signal", since with "signal" 
       there is a timing window that allows other threads within the
       process to run with the SIG_DFL handler and this could cause the
       process to abort before the traceback has completed.

       The "static volatile sig_atomic_t" variable is used to ensure that
       only 1 thread is updating it at any one time.

       Since the CPU does not normaly create a SIGFPE signal on a 
       floating-point error, the program initialises floating-point error
       trapping (for TRP_INVALID, TRP_DIV_BY_ZERO and TRP_OVERFLOW) if
       trapping for SIGFPE is requested. See man pages for "fp_trap" and
       "fp_enable".

       This subroutine works for serial code and also for multi-threaded
       codes, such as those produced using OpenMP.

 PARAMETERS
       create_dump - A pointer to an "int". If the integer has the value:
                     0        - no core dump is produced;
                     non-zero -  a core dump is produced.
       signals     - An array of signals to trap. This list is terminated
                     by an array element containing the value 0 (zero/NULL).
                     If the first element of the array is zero, then this
                     is taken to mean that the default set of signals,
                     stored in 'signals_to_trap' below, is to be used.

 RETURN VALUES
       -1 - an error occurred;
        0 - Floating-point error trapping was not set;
       >0 - The mode of floating-point trapping set
            the mode returned is usually 1 (FP_TRAP_SYNC).

 FORTRAN EXAMPLE
       PROGRAM IFS_MODEL

       INTEGER*4  CORE_DUMP_FLAG, IRETURN, SIGNALS(1), SIGNAL_TRAP

       CORE_DUMP_FLAG = 0
       SIGNALS(1)     = 0

       IRETURN = SIGNAL_TRAP(CORE_DUMP_FLAG, SIGNALS)

       IF (IRETURN .LT. 0) THEN
          PRINT *, 'ERROR'
       ELSEIF (IRETURN .EQ. 0) THEN
          PRINT *, 'FPE TRAPPING IS NOT SET'
       ELSE
          PRINT *, 'FPE TRAPPING MODE =', IRETURN
       ENDIF

       CALL SUB()

       END

       SUBROUTINE sub()

       REAL A

       A=sqrt(LOG(rand()))

       PRINT *, A

       RETURN

       END

 C EXAMPLE

   Save the file in signal_test.c and compile with

   xlc -q64 -g signal_test.c -o signal_test /usr/local/lib/libec.a  -lm -lxlf90

   then run

   ./signal_test

   Output should be something like:

FPE TRAPPING MODE = 1

  Signal received: SIGFPE - Floating-point exception
    Signal generated for floating-point exception:

  Traceback:
    Location 0x0000377c
    Offset 0x00000044 in procedure raise
    Offset 0x000000b0 in procedure fp_raise_xcp
    Offset 0x00000088 in procedure __sqrt_raise_xcp
    Offset 0x00000120 in procedure sqrt
    Offset 0x0000001c in procedure sub
    Offset 0x00000068 in procedure main
    --- End of call chain ---

#include <stdlib.h>
#include <signal.h>
#include <fptrap.h>
#include <memory.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>

extern int signal_trap(int32_t *, int32_t *);
void sub(void);

int
main(int argc, char *argv[])
{

   int32_t  *core_dump_flag = NULL;
   int32_t  *signals = NULL;
   int32_t  ireturn = 0;

   ireturn = signal_trap(core_dump_flag, signals);

   if (ireturn < 0) {
          printf("ERROR\n"); 
   } else if (ireturn == 0) {
          printf("FPE TRAPPING IS NOT SET\n");
   } else {
          printf("FPE TRAPPING MODE = %d\n",ireturn);
   }

   sub();

   printf("%f\n",ireturn);

}

void
sub(void) {

   float a;

   a=sqrt(log(rand()));

   printf("%f\n",a);

   return;

}

 
************************************************************************/

/*
  -------------------------------------------------------------------
  'signals_to_trap' is an array containing the signals that are to be
  caught and passed to the signal handler. You can add more signals,
  just make sure that the list is terminated by 'NULL'.
  -------------------------------------------------------------------
*/

signal_t signals_to_trap[] = {
		              SIGFPE,
                              SIGILL,
		              SIGBUS,
		              SIGSEGV,
		              SIGXCPU,
/* <<<<< INSERT NEW SIGNALS BETWEEN HERE >>>>>*/

/* <<<<<             AND HERE            >>>>>*/
		              NULL
		             };

   struct sigaction sa;
   signal_t         *sig, *sigs;
   int              ret = 0;
   char             errmsg[1024];

   signal_trap_dump = *create_dump;
   sigs = sig       = (*signals != NULL ? signals : signals_to_trap);

   memset(&sa, 0, sizeof(struct sigaction));
   sa.sa_flags      = SA_SIGINFO;
   sa.sa_sigaction  = signal_trace;
   sigemptyset(&sa.sa_mask);

/*
  -----------------------------------------------------------------------
  Add each signal to the "signal mask", as we don't want multiple signals
  interrupting the signal handler.
  -----------------------------------------------------------------------
*/

   while (*sig != NULL)
   {
      sigaddset(&sa.sa_mask, *sig++);
   }

   sig = sigs;
   while (*sig != NULL)
   {
      if (sigaction(*sig, &sa, NULL) == -1)
      {
         sprintf(errmsg, "Call to 'sigaction' in routine 'signal_trap' for signal number %d failed\n (line %d in file %s)\n",
	                 (int)*sig, __LINE__, __FILE__);
	 perror(errmsg);
	 exit(1);
      }

/*
  -----------------------------------------------------------------------
  If we are trapping Floating-Point Error, then set the processor in SYNC
  modes and enable TRP_INVALID, TRP_DIV_BY_ZERO and TRP_OVERFLOW.
  -----------------------------------------------------------------------
*/

      if (*sig++ == SIGFPE)
      {
         if (((ret = fp_trap(FP_TRAP_FASTMODE)) == FP_TRAP_UNIMPL) ||
                                           (ret == FP_TRAP_ERROR))
         {
            sprintf(errmsg, "Call to 'fp_trap' in signal_trap failed (return code = %d)\n (line %d in file %s)\n",
                            ret, __LINE__, __FILE__);
            perror(errmsg);
            exit(1);
         }

         fp_enable(TRP_INVALID | TRP_DIV_BY_ZERO | TRP_OVERFLOW);
      }
   }

   return ret;
}

/*===================================================================*/
void signal_trace(int signo, siginfo_t *sigcode, void *sigcontextptr)
/*===================================================================*/
{
/*
  -------------------------------------------------------------------------
  This is the actual signal-handler that is called by the Kernel when one
  of the trapped SIGNALs occurs. All but the first thread that cause it to
  be called are put to sleep for about 5 mins, giving the first thread time
  to print out the traceback and optionally dump the memory before exiting. 
  We set an alarm call for 2 minutes, to ensure that if we start looping,
  which can happen if the stack is overwritten, we will not loop forever.
  -------------------------------------------------------------------------
*/

   (void) alarm(ALARM_CALL);
   (void) signal(SIGALRM, SIG_DFL);

   if (signal_trap_incr++)
   {
      sleep(ALARM_CALL+10);
   }
   else
   {
      setlinebuf(stderr);
      setlinebuf(stdout);

      if (signal_trap_dump) 
      {
         xl__trcedump(signo, sigcode, sigcontextptr);
      }
      else
      {
         xl__trce(signo, sigcode, sigcontextptr);
      }
   }
}

#else
void signal_trap() { }
void signal_trap_() { }
#endif 
