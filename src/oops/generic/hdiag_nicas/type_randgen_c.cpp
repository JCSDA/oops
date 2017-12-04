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
#else
        srand(time(0));
#endif
    }
    else {
#if __cplusplus > 199711L
        gen = new std::mt19937;
#else
        srand(1);
#endif
    }
}

// Destructor
randGen::~randGen(){}

// Random integer generator
void randGen::rand_integer(int binf, int bsup, int *ir) const{
#if __cplusplus > 199711L
    // Initialize uniform distribution
    std::uniform_int_distribution<int> dis(binf,bsup);

    // Generate random integer
    *ir=dis(*gen);
#else
    // Generate random integer
    int range=bsup-binf+1;
    *ir=rand()%range+binf;
#endif
    return;
}

// Random real generator
void randGen::rand_real(double binf, double bsup, double *rr) const{
#if __cplusplus > 199711L
    // Initialize uniform distribution
    std::uniform_real_distribution<double> dis(binf,bsup);

    // Generate random real
    *rr=dis(*gen);
#else
    // Generate random real
    *rr=binf+(rand()/(RAND_MAX/(bsup-binf)));
#endif
    return;
}

// Sampling initialization
void randGen::initialize_sampling(int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]) const{
#if __cplusplus > 199711L
    // Initialize uniform distribution
    std::uniform_int_distribution<int> dis(0,n-1);
#endif

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
#if __cplusplus > 199711L
            ir=dis(*gen);
#else
            ir=rand()%n;
#endif
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
#if __cplusplus > 199711L
                ir=dis(*gen);
#else
                ir=rand()%n;
#endif
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

