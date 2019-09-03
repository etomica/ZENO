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

ResultsVirial::
    ResultsVirial(int numThreads, double refIntegral)
    :refAverage(NULL),
    refOverlapAverage(NULL),
    targetAverage(NULL),
    targetOverlapAverage(NULL),
    refAverageReduced(0),
    refOverlapAverageReduced(0),
    targetAverageReduced(0),
    targetOverlapAverageReduced(0),
    numThreads(numThreads),
    reduced(true),
    refNumSteps(NULL),
    targetNumSteps(NULL),
    refNumStepsReduced(0),
    targetNumStepsReduced(0),
    refIntegral(refIntegral){
    refAverage = new Uncertain<double>[numThreads];
    refOverlapAverage = new Uncertain<double>[numThreads];
    targetAverage = new Uncertain<double>[numThreads];
    targetOverlapAverage = new Uncertain<double>[numThreads];
    refNumSteps = new long long[numThreads];
    targetNumSteps = new long long[numThreads];

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverage[threadNum] = 0.0;
        refOverlapAverage[threadNum] = 0.0;
        targetAverage[threadNum] = 0.0;
        targetOverlapAverage[threadNum] = 0.0;
        refNumSteps[threadNum] = 0;
        targetNumSteps[threadNum] = 0;
    }
}

ResultsVirial::
~ResultsVirial() {
    delete[] refAverage;
    delete[] refOverlapAverage;
    delete[] targetAverage;
    delete[] targetOverlapAverage;
    delete[] refNumSteps;
    delete[] targetNumSteps;
}

void
ResultsVirial::
putData(int threadNum,
        MeterOverlap<double> * refMeter,
        MeterOverlap<double> * targetMeter){

    reduced = false;

    double ** stats = refMeter->getStatistics();
    refAverage[threadNum]=Uncertain<double>(stats[0][MeterOverlap<double>::AVG_AVG], pow(stats[0][MeterOverlap<double>::AVG_ERR],2));
    refOverlapAverage[threadNum]=Uncertain<double>(stats[1][MeterOverlap<double>::AVG_AVG], pow(stats[1][MeterOverlap<double>::AVG_ERR],2));
    refNumSteps[threadNum] = refMeter->getNumSamples();
    stats = targetMeter->getStatistics();
    targetAverage[threadNum]=Uncertain<double>(stats[0][MeterOverlap<double>::AVG_AVG], pow(stats[0][MeterOverlap<double>::AVG_ERR],2));
    targetOverlapAverage[threadNum]=Uncertain<double>(stats[1][MeterOverlap<double>::AVG_AVG], pow(stats[1][MeterOverlap<double>::AVG_ERR],2))/refMeter->getAlpha()[0];
    targetNumSteps[threadNum] = targetMeter->getNumSamples();
}

void
ResultsVirial::
reduce() {
    if (reduced) {
        return;
    }
    refAverageReduced = 0.0;
    refOverlapAverageReduced = 0.0;
    targetAverageReduced = 0.0;
    targetOverlapAverageReduced = 0.0;
    refNumStepsReduced = 0;
    targetNumStepsReduced = 0;

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverageReduced += refAverage[threadNum];
        refOverlapAverageReduced += refOverlapAverage[threadNum];
        refNumStepsReduced += refNumSteps[threadNum];
        targetAverageReduced += targetAverage[threadNum];
        targetOverlapAverageReduced += targetOverlapAverage[threadNum];
        targetNumStepsReduced += targetNumSteps[threadNum];
    }
    reduced = true;
}

long long
ResultsVirial::
getNumSteps() const{
    assert(reduced);
    return refNumStepsReduced + targetNumStepsReduced;
}

double
ResultsVirial::
getRefFrac() const{
    assert(reduced);
    return (double)refNumStepsReduced/(refNumStepsReduced + targetNumStepsReduced);
}

Uncertain<double>
ResultsVirial::
getRefAverageReduced() const {
    assert(reduced);
    return refAverageReduced;
}

Uncertain<double>
ResultsVirial::
getRefOverlapAverageReduced() const {
    assert(reduced);
    return refOverlapAverageReduced;
}

Uncertain<double>
ResultsVirial::
getTargetAverageReduced() const {
    assert(reduced);
    return targetAverageReduced;
}

Uncertain<double>
ResultsVirial::
getTargetOverlapAverageReduced() const {
    assert(reduced);
    return targetOverlapAverageReduced;
}

double
ResultsVirial::
getRefIntegral() const {
    return refIntegral;
}