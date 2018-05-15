//---------------------------------------------------------------------
/// Purpose: random numbers generator class header
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------

class rng {
    public:
        // Constructor
        rng(int default_seed);

        // Destructor
        ~rng();

        // Reseed generator
        void rng_reseed(int default_seed);

        // Random integer generator
        void rand_integer(int binf, int bsup, int *ir);

        // Random real generator
        void rand_real(double binf, double bsup, double *rr);

        // Sampling initialization
        void initialize_sampling(int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]);
    private:
        // Linear congruential generator
        unsigned long int a_;
        unsigned long int c_;
        unsigned long int m_;
        unsigned long int seed_;
        double lcg();
};
