//---------------------------------------------------------------------
/// Purpose: random number generator class implementation
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "type_randgen.h"
#include "type_randgen.hpp"
#include "external/Cover_Tree.h"
#include "external/Cover_Tree_Point.h"
#include <ostream>
#include <iomanip>
#include <cmath>
#if __cplusplus > 199711L
#include <random>
#endif

using namespace std;

// Constructor
randGen::randGen(int default_seed) {
    // Initialize random number generator
    if (default_seed==0) {
#if __cplusplus > 199711L
        std::random_device rd;
        gen = new std::mt19937(rd());
        version_ = 1;
#else
        seed = clock();
        version_ = 0;
#endif
    }
    else {
        seed = -14051987;
        version_ = 0;
    }
}

// Destructor
randGen::~randGen(){}

// Random integer generator
void randGen::rand_integer(int binf, int bsup, int *ir) {
    if (version_==1) {
#if __cplusplus > 199711L
        // Initialize uniform distribution
        std::uniform_int_distribution<int> dis(binf,bsup);

        // Generate random integer
        *ir=dis(*gen);
#endif
    }
    else {
        // Generate random integer
        int range=bsup-binf+1;
        double r;
        r = ran3();
        r *= range;
        *ir=binf+(int)r;
    }
    return;
}

// Random real generator
void randGen::rand_real(double binf, double bsup, double *rr) {
    if (version_==1) {
#if __cplusplus > 199711L
        // Initialize uniform distribution
        std::uniform_real_distribution<double> dis(binf,bsup);

        // Generate random real
        *rr=dis(*gen);
#endif
    }
    else {
       // Generate random real
       *rr=binf+ran3()*(bsup-binf);
    }
    return;
}

// Sampling initialization
void randGen::initialize_sampling(int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]) {
    // Declaration
    int ir;

    // Copy mask (updated then)
    int mask_copy[n];
    for(int i=0;i<n;i++) {
        mask_copy[i]=mask[i];
    }

    // Initialize tree
    CoverTree<CoverTreePoint> cTree(1.0e10);
    int is=0;

    // Fill the tree
    while(is<ns) {
        // Find new point
        double distmax=0.0;
        int irmax=-1;
        int itry=0;
        while(itry<ntry) {
            // Generate random real
            rand_integer(0,n-1,&ir);
            if(mask_copy[ir]==1) {
                if(is>0) {
                    // Find nearest neighbor
                    vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(-999,lon[ir],lat[ir]),1));

                    // Check distance
                    double dist=neighbors[0].getDist()/sqrt(0.5*(pow(rh[neighbors[0].getIndex()],2.0)+pow(rh[ir],2.0)));
                    if(dist>distmax) {
                        distmax=dist;
                        irmax=ir;
                    }
                }
                else {
                    irmax=ir;
                }
            }
            itry++;
        }

        // Insert point
        if(irmax>0) {
            ihor[is]=irmax+1;
            cTree.insert(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]));
            mask_copy[irmax]=0;
            is++;
        }
    }

    // Get minimum distance
    double distmininit=HUGE_VAL;
    for(int is=0;is<ns;is++) {
        vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]),2));
        double dist=neighbors[1].getDist()/sqrt(0.5*(pow(rh[neighbors[1].getIndex()],2.0)+pow(rh[ihor[is]-1],2.0)));
        if(dist<distmininit) {
            distmininit=dist;
        }
    }

    if (nrep>0) {
        // Improve sampling with replacements
        int irep=0;
        while(irep<nrep) {
            // Get minimum distance
            double distmin=HUGE_VAL;
            int ismin=0;
            for(int is=0;is<ns;is++) {
                vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]),2));
                double dist=neighbors[1].getDist()/sqrt(0.5*(pow(rh[neighbors[1].getIndex()],2.0)+pow(rh[ihor[is]-1],2.0)));
                if(dist<distmin) {
                    distmin=dist;
                    ismin=is;
                }
            }

            // Remove point
            cTree.remove(CoverTreePoint(ismin,lon[ihor[ismin]-1],lat[ihor[ismin]-1]));

            // Find new point
            double distmax=0.0;
            int irmax=-1;
            int itry=0;
            while(itry<ntry) {
                // Generate random number
                rand_integer(0,n-1,&ir);
                if(mask_copy[ir]==1) {
                    // Find nearest neighbor
                    vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(-999,lon[ir],lat[ir]),1));

                    // Check distance
                    double dist=neighbors[0].getDist()/sqrt(0.5*(pow(rh[neighbors[0].getIndex()],2.0)+pow(rh[ir],2.0)));
                    if((dist>distmax) && (dist>distmininit)) {
                        distmax=dist;
                        irmax=ir;
                    }
                }
                itry++;
            }

            // Insert point
            if(irmax>0) {
                // Insert new point
                ihor[ismin]=irmax+1;
                cTree.insert(CoverTreePoint(ismin,lon[ihor[ismin]-1],lat[ihor[ismin]-1]));
                mask_copy[irmax]=0;
            }
            else {
                // Re-insert old point
                cTree.insert(CoverTreePoint(ismin,lon[ihor[ismin]-1],lat[ihor[ismin]-1]));
            }
            irep++;
        }

        // Find final minimum distance
        double distmin=HUGE_VAL;
        for(int is=0;is<ns;is++) {
            vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]),2));
            double dist=neighbors[1].getDist()/sqrt(0.5*(pow(rh[neighbors[1].getIndex()],2.0)+pow(rh[ihor[is]-1],2.0)));
            if(dist<distmin) {
                distmin=dist;
            }
        }
    }
    return;
}

// ran3 generator
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
double randGen::ran3() {
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;

    if (seed < 0 || iff == 0) {
        iff=1;
        mj=MSEED-(seed < 0 ? -seed : seed);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;k++)
            for (i=1;i<=55;i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext=0;
        inextp=31;
        seed=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return (double)mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
