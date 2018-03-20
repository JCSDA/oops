//---------------------------------------------------------------------
/// Purpose: random numbers generator class header
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------

class rng {
    public:
        // Constructor
        rng(unsigned long int default_seed);

        // Destructor
        ~rng();

        // Reseed generator
        void rng_reseed(unsigned long int seed);

        // Random integer generator
        void rand_integer(int binf, int bsup, int *ir);

        // Random real generator
        void rand_real(double binf, double bsup, double *rr);

        // Sampling initialization
        void initialize_sampling(int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]);
    private:
        // Linear congruential generator
        unsigned long int a_=1103515245;
        unsigned long int c_=12345;
        unsigned long int m_=2147483648;
        unsigned long int seed_;
        double lcg();
};
