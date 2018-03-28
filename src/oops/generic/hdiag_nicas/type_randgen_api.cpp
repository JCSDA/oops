//---------------------------------------------------------------------
/// Purpose: random number generator class API
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "type_randgen.h"
#include "type_randgen.hpp"

using namespace std;

// Constructor
RANDGEN* create_randgen(unsigned long int default_seed){
    return new randGen(default_seed);
}

// Destructor
void delete_randgen(RANDGEN* randgen){
    delete randgen;
}

// Reseed generator
void reseed_randgen(RANDGEN* randgen, unsigned long int seed) {
    randgen->reseed_randgen(seed);
}

// Get version
void get_version(const RANDGEN* randgen, int *version) {
    randgen->get_version(version);
}

// Random integer generator
void rand_integer(RANDGEN* randgen, int binf, int bsup, int *ir) {
    randgen->rand_integer(binf, bsup, ir);
}

// Random real generator
void rand_real(RANDGEN* randgen, double binf, double bsup, double *rr) {
    randgen->rand_real(binf, bsup, rr);
}

// Sampling initialization
void initialize_sampling(RANDGEN* randgen, int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]) {
    randgen->initialize_sampling( n, lon, lat, mask, rh, ntry, nrep, ns, ihor);
}
