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

#include "ResultsVirial.h"

template <class T, class RandomNumberGenerator>
ResultsVirial<T, RandomNumberGenerator>::
    ResultsVirial(int numThreads)
    :refAverage(NULL),
    refOverlapAverage(NULL),
    targetAverage(NULL),
    targetOverlapAverage(NULL),
    refAverageReduced(0),
    refOverlapAverageReduced(0),
    targetAverageReduced(0),
    targetOverlapAverageReduced(0),
    numThreads(numThreads),
    reduced(true){
    refAverage = new Uncertain<double>[numThreads];
    refOverlapAverage = new Uncertain<double>[numThreads];
    targetAverage = new Uncertain<double>[numThreads];
    targetOverlapAverage = new Uncertain<double>[numThreads];

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverage[threadNum] = 0.0;
        refOverlapAverage[threadNum] = 0.0;
        targetAverage[threadNum] = 0.0;
        targetOverlapAverage[threadNum] = 0.0;
    }
}
template <class T, class RandomNumberGenerator>
ResultsVirial<T, RandomNumberGenerator>::
~ResultsVirial() {
    delete[] refAverage;
    delete[] refOverlapAverage;
    delete[] targetAverage;
    delete[] targetOverlapAverage;
}

template <class T, class RandomNumberGenerator>
void
ResultsVirial<T, RandomNumberGenerator>::
putData(int threadNum,
        MeterOverlap<T, RandomNumberGenerator> refMeter,
        MeterOverlap<T, RandomNumberGenerator> targetMeter){

    reduced = false;

    double ** stats = refMeter.getStatistics();
    refAverage[threadNum]=Uncertain<double>(stats[0][MeterOverlap<T, RandomNumberGenerator>::AVG_AVG], stats[0][MeterOverlap<T, RandomNumberGenerator>::AVG_ERR]);
    refOverlapAverage[threadNum]=Uncertain<double>(stats[1][MeterOverlap<T, RandomNumberGenerator>::AVG_AVG], stats[1][MeterOverlap<T, RandomNumberGenerator>::AVG_ERR]);
    stats = targetMeter.getStatistics();
    targetAverage[threadNum]=Uncertain<double>(stats[0][MeterOverlap<T, RandomNumberGenerator>::AVG_AVG], stats[0][MeterOverlap<T, RandomNumberGenerator>::AVG_ERR]);
    targetOverlapAverage[threadNum]=Uncertain<double>(stats[1][MeterOverlap<T, RandomNumberGenerator>::AVG_AVG], stats[1][MeterOverlap<T, RandomNumberGenerator>::AVG_ERR]);
}

template <class T, class RandomNumberGenerator>
void
ResultsVirial<T, RandomNumberGenerator>::
reduce() {
    if (reduced) {
        return;
    }
    refAverageReduced = 0.0;
    refOverlapAverageReduced = 0.0;
    targetAverageReduced = 0.0;
    targetOverlapAverageReduced = 0.0;

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverageReduced += refAverage[threadNum];
        refOverlapAverageReduced += refOverlapAverage[threadNum];
        targetAverageReduced += targetAverage[threadNum];
        targetOverlapAverageReduced += targetOverlapAverage[threadNum];
    }
    reduced = true;
}