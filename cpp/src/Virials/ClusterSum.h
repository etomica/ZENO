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

#include "OverlapTester.h"

/// Computes cluster sum.
///
template <class T,
        class RandomNumberGenerator>
class IntegratorMSMC;

template <class T,
        class RandomNumberGenerator>
class ClusterSum {
 public:
    ClusterSum(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, OverlapTester<T> const & overlapTester);
    virtual ~ClusterSum();
    virtual double value() = 0;

 protected:
    IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC;
    OverlapTester<T> const & overlapTester;
};

///Sub class of ClusterSum to compute cluster sum for chains.
///

template <class T,
        class RandomNumberGenerator>
class ClusterSumChain : public ClusterSum<T, RandomNumberGenerator> {
public:
    ClusterSumChain(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, OverlapTester<T> const & overlapTester, double ringFac, double chainFac);
    ~ClusterSumChain();
    double value();

private:
      double ringFac;
      double chainFac;
};

///Sub class of ClusterSum to compute cluster sum using Wheatley Recursion.
///
template <class T,
        class RandomNumberGenerator>
class ClusterSumWheatleyRecursion : public ClusterSum<T, RandomNumberGenerator>{
public:
    ClusterSumWheatleyRecursion(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, OverlapTester<T> const & overlapTester);
    ~ClusterSumWheatleyRecursion();
    double value();

private:
    double preFac;
};
#endif

