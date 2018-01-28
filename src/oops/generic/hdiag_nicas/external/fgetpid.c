//----------------------------------------------------------------------
// fgetpid_
// A. C. Lorenz (U.K. Met-Office)
// Modified by Benjamin Menetrier for hdiag_nicas
//----------------------------------------------------------------------
#include <unistd.h>

void fgetpid_(int *id)
{
    *id = (int)getpid();
    return;
}
