//---------------------------------------------------------------------
/// Purpose: random numbers generator class header
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
extern "C"
{
    class rng;
    typedef rng RNG;

    // Constructor
    RNG* rng_create(unsigned long int default_seed);

    // Destructor
    void rng_delete(RNG* rng);

    // Reseed generator
    void rng_reseed(RNG* rng, unsigned long int seed);

    // Random integer generator
    void rand_integer(RNG* rng, int binf, int bsup, int *ir);

    // Random real generator
    void rand_real(RNG* rng, double binf, double bsup, double *rr);

    // Sampling initialization
    void initialize_sampling(RNG* rng, int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]);
}
