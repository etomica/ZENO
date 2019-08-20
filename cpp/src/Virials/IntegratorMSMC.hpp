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

#include "IntegratorMSMC.h"


///  Constructs the class to perform a step.
///

template <class T,
  class RandomNumberGenerator>
IntegratorMSMC<T,
               RandomNumberGenerator>::
IntegratorMSMC(Parameters const * parameters,
                int threadNum,
                Timer const * totalTimer,
                RandomNumberGenerator * randomNumberGenerator,
                std::vector<Sphere<double> *> & boundingSpheres,
                std::vector<int> & numParticles,
                std::vector<MixedModel<T> *> & models,
                OverlapTester<T> const & overlapTester,
                std::vector<MCMove<T, RandomNumberGenerator>> & mcMoves,
                std::vector<double> & moveProbs,
                ClusterSum<T, RandomNumberGenerator> * clusterSum,
                MeterOverlap<T, RandomNumberGenerator> & meterOverlap) :
              parameters(parameters),
              threadNum(threadNum),
              totalTimer(totalTimer),
              randomNumberGenerator(randomNumberGenerator),
              boundingSpheres(boundingSpheres),
              numParticles(numParticles),
              overlapTester(overlapTester),
              mcMoves(mcMoves),
              moveProbs(moveProbs),
              clusterSum(clusterSum),
              meterOverlap(meterOverlap){
    for(int i = 0; i < numParticles.size(); ++i)
    {
        for (int j =0; j < numParticles[i]; ++j)
        {
            particles.push_back(new Particle<T> (models[i], boundingSpheres[i]));
        }
    }

}

template <class T,
  class RandomNumberGenerator>
IntegratorMSMC<T,
               RandomNumberGenerator>::
  ~IntegratorMSMC() {

}

/// Carries out a step. A step performed by the integrator consists of selecting a MCMove, performing the trial defined by MCMove and collecting data and statistics.
///
template <class T,
  class RandomNumberGenerator>
void
IntegratorMSMC<T,
               RandomNumberGenerator>::
doStep(){
    double random = randomNumberGenerator->getRandIn01();
    double cumProb = 0.0;
    int currentMove = -1;
    for(int i = 0; i < mcMoves.size(); ++i)
    {
        cumProb += moveProbs[i];
        if(random < cumProb)
        {
            currentMove = i;
        }
        mcMoves[i].doTrial();
    }
    meterOverlap.collectData();
}

/// Returns particles.
///
template <class T,
        class RandomNumberGenerator>
std::vector<Particle<T> *>
IntegratorMSMC<T,RandomNumberGenerator>::
getParticles(){
    return particles;
}

/// Returns random number generator.
///
template <class T,
        class RandomNumberGenerator>
RandomNumberGenerator *
IntegratorMSMC<T,RandomNumberGenerator>::
getRandomNumberGenerator(){
    return randomNumberGenerator;
}

/// Returns cluster sum.
///
template <class T,
        class RandomNumberGenerator>
ClusterSum<T, RandomNumberGenerator> *
IntegratorMSMC<T,RandomNumberGenerator>::
getClusterSum(){
    return clusterSum;
}

/// Returns random utilities.
///
template <class T,
        class RandomNumberGenerator>
RandomUtilities<T, RandomNumberGenerator> *
IntegratorMSMC<T,RandomNumberGenerator>::
getRandomUtilities(){
    return & randomUtilities;
}

