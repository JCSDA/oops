//---------------------------------------------------------------------
/// Purpose: random number generator class API
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "type_rng.h"
#include "type_rng.hpp"

using namespace std;

// Constructor
RNG* rng_create(unsigned long int default_seed){
    return new rng(default_seed);
}

// Destructor
void rng_delete(RNG* rng){
    delete rng;
}

// Reseed generator
void rng_reseed(RNG* rng, unsigned long int seed) {
    rng->rng_reseed(seed);
}

// Random integer generator
void rand_integer(RNG* rng, int binf, int bsup, int *ir) {
    rng->rand_integer(binf, bsup, ir);
}

// Random real generator
void rand_real(RNG* rng, double binf, double bsup, double *rr) {
    rng->rand_real(binf, bsup, rr);
}

// Sampling initialization
void initialize_sampling(RNG* rng, int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]) {
    rng->initialize_sampling( n, lon, lat, mask, rh, ntry, nrep, ns, ihor);
}
