// ================================================================
//
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
//
// ================================================================

// ================================================================
// 
// Authors:
// Created:
//
// ================================================================

#ifndef CLUSTER_SUM_H
#define CLUSTER_SUM_H

#include "Potential.h"

namespace zeno {

/// Computes cluster sum.
///

template <class T>
class ClusterSum {
 public:
    ClusterSum(std::vector<Particle<T> *> * particles);
    virtual ~ClusterSum();
    virtual double value() = 0;
    virtual int numValues() {
      return values.size();
    }
    virtual std::vector<double> getValues() {
      values[0] = value();
      return values;
    }

 protected:
    std::vector<Particle<T> *> * particles;
    std::vector<double> values;
};

///Sub class of ClusterSum to compute cluster sum for chains.
///

template <class T>
class ClusterSumChain : public ClusterSum<T> {
public:
    ClusterSumChain(std::vector<Particle<T> *> * particles, double diameter, double ringFac, double chainFac);
    ~ClusterSumChain();
    double value();

private:
    double diameter;
    double ringFac;
    double chainFac;
};

///Sub class of ClusterSum to compute cluster sum using Wheatley Recursion.
///
template <class T>
class ClusterSumWheatleyRecursion : public ClusterSum<T>{
public:
    ClusterSumWheatleyRecursion(std::vector<Particle<T> *> * particles, Potential<T> const * potential, double temperature, int nDerivatives);
    ~ClusterSumWheatleyRecursion();
    double value();
    virtual std::vector<double> getValues();

private:
    Potential<T> const * potential;
    std::vector<double> derivativeValues;
    double temperature;
    double preFac;
};


/// Constructs the class to compute cluster sum.
///

using namespace zeno;

template <class T>
ClusterSum<T>::
ClusterSum(std::vector<Particle<T> *> * particles):
particles(particles){
  values.resize(1);
}

template <class T>
ClusterSum<T>::
  ~ClusterSum() {
}

/// Constructs a sub class of ClusterSum to compute cluster sum for chains.
///
template <class T>
ClusterSumChain<T>::
ClusterSumChain(std::vector<Particle<T> *> * particles, double diameter,
                double ringFac, double chainFac):
        ClusterSum<T>(particles), diameter(diameter), ringFac(ringFac), chainFac(chainFac) {
    ClusterSum<T>::values.resize(1);
}

template <class T>
ClusterSumChain<T>::
~ClusterSumChain() {
}

/// Computes cluster sum for chains.
///
template <class T>
double
ClusterSumChain<T>::
value(){
    const int n = ClusterSum<T>::particles->size();
    const int nf1 = (1 << (n - 1));
    double fValues[n][n];
    for(int iMol1 = 0; iMol1 < n; ++iMol1){
        for(int iMol2 = iMol1 + 1; iMol2 < n; ++iMol2){
            Vector3<T> x = ClusterSum<T>::particles->at(iMol1)->getCenter();
            Vector3<T> y = ClusterSum<T>::particles->at(iMol2)->getCenter();
            Vector3<T> distCenterVec = x - y;
            T distCenterSqr = distCenterVec.getMagnitudeSqr();
            bool overlapped = distCenterSqr < diameter*diameter;
            fValues[iMol1][iMol2] = fValues[iMol2][iMol1] = overlapped ? 1 : 0;
        }
    }
    double nC[n-1][nf1];
    //nC[m][i] is the number of chains beginning at last vertex and ending at m, traversing all points in i
    //Start with all pairwise paths from last vertex to each other vertex
    for(int m = 0; m < (n - 1); ++m){
        nC[m][(1 << m)] = fValues[m][n - 1];
    }
    //All other paths
    for(int i = 1; i < nf1; ++i) {//i excludes the last bit, which is implicit in all partitions
        //the following two loops generate all pairs formed by each bit in i with each bit not in i
        //loop over bits not in i; start with full complement of i (i^(nf1-1)), and in each iteration
        //get lowest bit (im=(iC&-iC)) and strip it from complement (iC^=im) until complement is empty (iC=0)
        for(int iC = i ^ (nf1 - 1), im = (iC & -iC); iC > 0; iC ^= im, im = (iC & -iC)) {
            int m = log2(im);
            int iim = i|im;
            nC[m][iim] = 0;
            //loop over bits in i, in same manner as loop over complement
            for (int it = i, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
                int k = log2(ik);
                nC[m][iim] += fValues[m][k] * nC[k][i];
            }
        }//end for(iC)
    }//end for(i)
    double ringValue = 0;
    double chainValue = 0;
    if (ringFac != 0.0) {
        for (int m = 0; m < n - 1; m++) {
            ringValue += nC[m][nf1-1] * fValues[m][n-1];
        }
    }
    if (chainFac != 0.0) {
        //Sum chains in which last (n-1) vertex is not a leaf.
        //Consider all partitions, counting paths beginning in one partition and ending in its complement
        //Use same looping structure as employed above
        for (int iS = 1; iS < nf1; iS += 2) {//keep 1 in iS-partition to prevent double counting
            //loop over bits not in iS
            int iSComp = iS^(nf1-1);
            for (int iC = iSComp, im = (iC & -iC); iC > 0; iC ^= im, im = (iC & -iC)) {
                int m = log2(im);
                //loop over bits in iS
                for (int it = iS, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
                    int k = log2(ik);
                    chainValue += nC[m][iSComp] * nC[k][iS];
                }
            }
        }
        //Sum chains where last (n-1) vertex is a leaf
        for (int m = 0; m < n - 1; m++) {
            chainValue += nC[m][nf1-1];
        }
    }//end if(chainFrac)
    return chainFac*chainValue + ringFac*ringValue;
}

/// Constucts a sub class of ClusterSum to compute cluster sum using Wheatley Recursion.
///
template <class T>
ClusterSumWheatleyRecursion<T>::
ClusterSumWheatleyRecursion(std::vector<Particle<T> *> * particles, Potential<T> const * potential, double temperature, int nDerivatives):
ClusterSum<T>(particles), potential(potential), temperature(temperature) {
    const int n = ClusterSum<T>::particles->size();
    int factorial = 1;
    for (int m = 2; m <= n; ++m){
        factorial *= m;
    }
    preFac = -(n - 1.0)/factorial;
    ClusterSum<T>::values.resize(1+nDerivatives);
}

template <class T>
ClusterSumWheatleyRecursion<T>::
~ClusterSumWheatleyRecursion() {
}

/// Computes cluster sum using Wheatley Recursion.
///
template <class T>
double
ClusterSumWheatleyRecursion<T>::
value() {
    return getValues()[0];
}

/// Computes cluster sum using Wheatley Recursion.
///
template <class T>
std::vector<double>
ClusterSumWheatleyRecursion<T>::
getValues() {
    const int nDer = ClusterSum<T>::values.size() - 1;
    const int n = ClusterSum<T>::particles->size();
    const int nf = (1 << n);
    double fQ[nf][nDer+1], fC[nf][nDer+1];
    double fA[nf][nDer+1], fB[nf][nDer+1];
    double binomial[nDer+1][nDer+1];
    int factorial[nDer+1];
    factorial[1] = 1;
    for (int m=2; m<=nDer; m++) {
        factorial[m] = factorial[m-1]*m;
    }
    for (int m=0; m<=nDer; m++) {
        for (int l=0; l<=m; l++) {
            binomial[m][l] = factorial[m] / (factorial[l] * factorial[m-l]);
        }
    }
    for(int iMol1 = 0; iMol1 < n; ++iMol1){
        int i = 1 << iMol1;
        fQ[i][0] = 1.0;
        for(int iMol2 = iMol1 + 1; iMol2 < n; ++iMol2){
            double u = potential->energy2(ClusterSum<T>::particles->at(iMol1), ClusterSum<T>::particles->at(iMol2));
            double e = std::exp(-u/temperature);
            if (e < 1e-12) e = 0;
            int ij = i|(1<<iMol2);
            fQ[ij][0] = e;
            if (e == 0) {
              for (int m=1; m<=nDer; m++) {
                fQ[ij][m] = 0;
              }
              continue;
            }
            for (int m=1; m<=nDer; m++) {
              fQ[ij][m] = -fQ[ij][m-1]*u;
            }
        }
    }
    //generate all partitions and compute
    for (int i = 3; i < nf; ++i){
        int j = i & -i; //lowest bit in i
        if (i == j) continue; //1-point set
        int k = i & ~j;
        if (k == (k & -k)) {
            // 2-point set
            continue;
        }
        fQ[i][0] = fQ[k][0];
        if (fQ[i][0] == 0) {
            for (int m=1; m<=nDer; m++) {
                fQ[i][m] = 0;
            }
            continue;
        }
        for (int l = (j << 1); l < i; l = (l << 1)){
            if ( (l&i) == 0 ) continue;
            fQ[i][0] *= fQ[l|j][0];
        }
        if (fQ[i][0] == 0) {
            for (int m=1; m<=nDer; m++) {
                fQ[i][m] = 0;
            }
            continue;
        }
        double c = std::log(fQ[i][0])*temperature;
        for (int m=1; m<=nDer; m++) {
            fQ[i][m] = fQ[i][m-1] * c;
        }
    }
    //Compute the fC's
    for (int i = 1; i < nf; ++i){
        for (int m=0; m<=nDer; m++) {
            fC[i][m] = fQ[i][m];
        }
        int iLowBit = i & -i;
        int inc = iLowBit << 1;
        for (int j = iLowBit; j < i; j += inc){
            int jComp = i & ~j;
            while ((j|jComp) != i && j<i){
                int jHighBits = j ^ iLowBit;
                int jlow = jHighBits & -jHighBits;
                j += jlow;
                jComp = (i & ~j);
            }
            if (j==i) break;
            for (int m=0; m<=nDer; m++) {
                for (int l=0; l<=m; l++) {
                    fC[i][m] -= binomial[m][l]*fC[j][l] * fQ[jComp][m-l];
                }
            }
        }
    }
    //find fA1
    for (int i = 2; i < nf; i += 2){
        //all even sets don't contain 1
        for (int m=0; m<=nDer; m++) {
            fB[i][m] = fC[i][m];
        }
    }
    for (int m=0; m<=nDer; m++) {
        fA[1][m] = 0;
        fB[1][m] = fC[1][m];
    }
    for (int i = 3; i < nf; i += 2){
        //every set will contain 1
        for (int m=0; m<=nDer; m++) {
            fA[i][m] = 0;
            fB[i][m] = fC[i][m];
        }
        int ii = i - 1;//all bits in i but lowest
        int iLow2Bit = (ii & -ii);//next lowest bit
        int jBits = 1 | iLow2Bit;
        if (jBits == i) continue;
        int iii = ii ^ iLow2Bit; //i with 2 lowest bits off
        int jInc = (iii & -iii);//3rd lowest bit, alsso increment for j
        for (int j = jBits; j < i; j += jInc){//sum over partitions of i containing j Bits
            int jComp = (i & ~j);//subset of i complementing j
            while ((j|jComp) != i && j<i){
                int jHighBits = j ^ jBits;
                int jlow = jHighBits & -jHighBits;
                j += jlow;
                jComp = (i & ~j);
            }
            if (j == i) break;
            for (int m=0; m<=nDer; m++) {
                for (int l=0; l<=m; l++) {
                    fA[i][m] += binomial[m][l]*fB[j][l] * fC[jComp|1][m-l];
                }
            }
        }
        //remove from B graphs that contain articulation point 0.
        for (int m=0; m<=nDer; m++) {
            fB[i][m] -= fA[i][m];
        }
    }
    for (int v = 1; v < n; ++v){
        int vs1 = 1 << v;
        for (int i = vs1 + 1; i < nf; ++i){
            for (int m=0; m<=nDer; m++) {
                fA[i][m] = 0;
            }
            if ( (i & vs1) == 0 ) continue;
            int iLowBit = (i & -i);
            if ( iLowBit == i ) continue;
            int jBits;
            int ii = i ^ iLowBit;
            int iLow2Bit = (ii & -ii);
            if ( iLowBit!=vs1 && iLow2Bit!=vs1 ){
                //v is not in the lowest 2 bits
                jBits = iLowBit | vs1;
                //we can only increment by the 2nd lowest
                int jInc = iLow2Bit;
                for (int j = jBits; j < i; j += jInc){
                    if ( (j&jBits) != jBits ){
                        j |= vs1;
                        if (j==i) break;
                    }
                    int jComp = i & ~j;
                    while ((j|jComp) != i && j<i){
                        int jHighBits = j^jBits;
                        int jlow = jHighBits & -jHighBits;
                        j += jlow;
                        j |= vs1;
                        jComp = (i & ~j);
                    }
                    if (j==i) break;
                    for  (int m=0; m<=nDer; m++) {
                        for (int l=0; l<=m; l++) {
                            fA[i][m] += binomial[m][l]*fB[j][l] * (fB[jComp|vs1][m-l] + fA[jComp|vs1][m-l]);
                        }
                    }
                }
            }
            else{
                //lowest 2 bits contain v
                jBits = iLowBit | iLow2Bit;
                if (jBits == i) continue; // no bits left jComp
                int iii = ii ^ iLow2Bit;
                int jInc = ( iii & -iii);
                //at this point jBits has (lowest bit + v)
                for (int j = jBits; j < i; j += jInc){//sum over partitions of i
                    int jComp = i & ~j;
                    while ((j|jComp) != i && j<i){
                        int jHighBits = j^jBits;
                        int jlow = jHighBits & -jHighBits;
                        j += jlow;
                        jComp = (i & ~j);
                    }
                    if (j==i) break;
                    for (int m=0; m<=nDer; m++) {
                        for  (int l=0; l<=m; l++) {
                            fA[i][m] += binomial[m][l]*fB[j][l] * (fB[jComp|vs1][m-l] + fA[jComp|vs1][m-l]);
                        }
                    }
                }
            }
            for (int m=0; m<=nDer; m++) {
                fB[i][m] -= fA[i][m];
            }
        }
    }
    for (int m=0; m<=nDer; m++) {
        ClusterSum<T>::values[m] = preFac*fB[nf-1][m];
    }
    return ClusterSum<T>::values;
}

}
#endif
