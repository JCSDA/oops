//----------------------------------------------------------------------
// fgetpid_
// A. C. Lorenz (U.K. Met-Office)
// Modified by Benjamin Menetrier for BUMP
//----------------------------------------------------------------------
#include <unistd.h>

void fgetpid_(int *id)
{
    *id = (int)getpid();
    return;
}
