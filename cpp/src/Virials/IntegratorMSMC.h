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

#ifndef INTEGRATOR_MSMC_H
#define INTEGRATOR_MSMC_H

#include <vector>

#include "../Parameters.h"
#include "../Timer.h"
#include "../Geometry/Sphere.h"
#include "../Geometry/MixedModel.h"
#include "OverlapTester.h"
#include "MCMove.h"
#include "ClusterSum.h"
#include "MeterOverlap.h"
#include "Particle.h"
#include "RandomUtilities.h"

/// Performs a step. A step performed by the integrator consists of selecting a MCMove, performing the trial defined by MCMove and collecting data and statistics.
///

template <class T,
        class RandomNumberGenerator>
class MCMove;

template <class T,
        class RandomNumberGenerator>
class ClusterSum;

template <class T,
  class RandomNumberGenerator>
class IntegratorMSMC {
 public:
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
                MeterOverlap<T, RandomNumberGenerator> & meterOverlap);

  ~IntegratorMSMC();

  void doStep();
  std::vector<Particle<T> *> getParticles();
  RandomNumberGenerator * getRandomNumberGenerator();
  ClusterSum<T, RandomNumberGenerator> * getClusterSum();
  RandomUtilities<T, RandomNumberGenerator> * getRandomUtilities();
private:
  Parameters const * parameters;
  int threadNum;
  Timer const * totalTimer;
  RandomNumberGenerator * randomNumberGenerator;
  std::vector<Sphere<double> *> & boundingSpheres;
  std::vector<int> & numParticles;
  OverlapTester<T> const & overlapTester;
  std::vector<Particle<T> *> particles;
  std::vector<MCMove<T, RandomNumberGenerator>> & mcMoves;
  std::vector<double> & moveProbs;
  ClusterSum<T, RandomNumberGenerator> * clusterSum;
  RandomUtilities<T, RandomNumberGenerator> randomUtilities;
  MeterOverlap<T, RandomNumberGenerator> & meterOverlap;
};

#endif

