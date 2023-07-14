// ================================================================
//
// This software was developed by employees of the National Institute of
// Standards and Technology (NIST), an agency of the Federal Government.
// Pursuant to title 17 United States Code Section 105, works of NIST employees
// are not subject to copyright protection in the United States and are
// considered to be in the public domain. Permission to freely use, copy,
// modify, and distribute this software and its documentation without fee is
// hereby granted, provided that this notice and disclaimer of warranty appears
// in all copies.
//
// THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
// EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
// WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
// FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
// THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
// EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
// DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
// RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
// BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
// SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
// SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
// SERVICES PROVIDED HEREUNDER.
//
// Distributions of NIST software should also include copyright and licensing
// statements of any third-party software that are legally bundled with the
// code in compliance with the conditions of those licenses.
// 
// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Created: 2019-05-08
//
// ================================================================

#include "Zeno.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <thread>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "BoundingSphereGenerator.h"

#include "Walker/WalkerExterior.h"
#include "Walker/SamplerInterior.h"

#include "Potential.h"

#include "Virials/VirialAlpha.h"
#include "Virials/VirialProduction.h"

using namespace zeno;

Zeno::Zeno(std::vector<MixedModel<double>> * modelsToProcess, Potential<double> & potential)
  : mpiSize(1),
    mpiRank(0),
    model(),
    potential(potential),
    modelBoundingSphere(),
    resultsZeno(nullptr),
    resultsInterior(nullptr),
    resultsVirial(nullptr),
    initializeTimer(),
    preprocessTimer(),
    walkOnSpheresTimer(),
    walkOnSpheresReductionTimer(),
    interiorSamplingTimer(),
    interiorSamplingReductionTimer(),
    virialTimer(),
    virialReductionTimer(),
    totalTimer() {

  totalTimer.start();

#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif

  preprocessTimer.start();

  model.addMixedModel(&modelsToProcess->at(0));
  
  if (potential.getRigidHS()) model.preprocess();

  if (!model.isEmpty()) {
    modelBoundingSphere = BoundingSphereGenerator<double>::generate(model);
  }

  preprocessTimer.stop();  

  for (MixedModel<double> &mtp : *modelsToProcess) {
    Model m;
    m.addMixedModel(&mtp);
    allModels.push_back(m);
  }
}

Zeno::~Zeno() {
  totalTimer.stop();
}
  
Zeno::Status
Zeno::doWalkOnSpheres(ParametersWalkOnSpheres * parametersWalkOnSpheres,
		      ParametersResults * parametersResults) {

  if (model.isEmpty()) {
    return Status::EmptyModel;
  }
  
  initializeTimer.start();

  computeDefaultParameters(parametersWalkOnSpheres);
  computeDefaultParameters(parametersResults);

  long long numWalksInProcess = 
	computeNumInProcess(mpiSize, 
			    mpiRank, 
			    parametersWalkOnSpheres->getTotalNumWalks());
  
  std::vector<RandomNumberGenerator> threadRNGs;

  setupRNGs(parametersWalkOnSpheres->getNumThreads(),
	    parametersWalkOnSpheres->getSeed(),
	    &threadRNGs);

  delete resultsZeno;

  resultsZeno = nullptr;

  initializeTimer.stop();
  
  getWalkOnSpheresResults(numWalksInProcess,			
			  *parametersWalkOnSpheres,
			  *parametersResults,
			  parametersWalkOnSpheres->getLaunchSphere(),
			  model,
			  &threadRNGs,
			  &resultsZeno);

  return Status::Success;
}

Zeno::Status
Zeno::doInteriorSampling
(ParametersInteriorSampling * parametersInteriorSampling,
 ParametersResults * parametersResults) {

  if (model.isEmpty()) {
    return Status::EmptyModel;
  }
  
  initializeTimer.start();

  computeDefaultParameters(parametersInteriorSampling);
  computeDefaultParameters(parametersResults);

  long long numSamplesInProcess = 
	computeNumInProcess(mpiSize, 
			    mpiRank, 
			    parametersInteriorSampling->getTotalNumSamples());
  
  std::vector<RandomNumberGenerator> threadRNGs;

  setupRNGs(parametersInteriorSampling->getNumThreads(),
	    parametersInteriorSampling->getSeed(),
	    &threadRNGs);

  delete resultsInterior;

  resultsInterior = nullptr;

  initializeTimer.stop();
  
  getInteriorResults(numSamplesInProcess,			
		     *parametersInteriorSampling,
		     *parametersResults,
		     parametersInteriorSampling->getLaunchSphere(),
		     model,
		     &threadRNGs,
		     &resultsInterior);
  
  return Status::Success;
}

Zeno::Status
Zeno::doVirialSampling
(ParametersVirial * parametersVirial,
 ParametersResults * parametersResults) {

  if (model.isEmpty()) {
    return Status::EmptyModel;
  }
 
  initializeTimer.start();

  computeDefaultParameters(parametersVirial);
  computeDefaultParameters(parametersResults);

  long long numStepsInProcess =
	computeNumInProcess(mpiSize,
			    mpiRank,
			    parametersVirial->getSteps());
 
  std::vector<RandomNumberGenerator> threadRNGs;

  setupRNGs(parametersVirial->getNumThreads(),
	    parametersVirial->getSeed(),
	    &threadRNGs);

  delete resultsVirial;

  resultsVirial = nullptr;

  initializeTimer.stop();

  std::vector<BoundingSphere> boundingSpheres;
  for (Model m : allModels) {
    BoundingSphere modelBoundingSphere = BoundingSphereGenerator<double>::generate(model);
    Sphere<double> boundingSphere(modelBoundingSphere.getCenter(),
                                  modelBoundingSphere.getRadius());
    boundingSpheres.push_back(boundingSphere);
  }
 
  getVirialResults(numStepsInProcess,
		   *parametersVirial,
		   boundingSpheres,
		   potential,
		   &threadRNGs,
		   &resultsVirial);
 
  return Status::Success;
}

void
Zeno::getResults(ParametersResults * parametersResults,
		 Results * results) const {

  computeDefaultParameters(parametersResults);

  if (potential.getFlexible()) {
      parametersResults->setVirialFlexible();
  }
  
  ResultsCompiler resultsCompiler(*parametersResults);

  resultsCompiler.compile(resultsZeno,
  			  resultsInterior,
  			  resultsVirial,
  			  false,
			  results);
}

void
Zeno::getWalkOnSpheresHitPoints
(std::vector<Vector3<double> > const * * points,
 std::vector<Vector3<char> > const * * charges) const {

  if (resultsZeno != nullptr) {
    resultsZeno->gatherHitPoints();

    if (points != nullptr) {
      *points = resultsZeno->getPoints();
    }

    if (charges != nullptr) {
      *charges = resultsZeno->getCharges();
    }
  }
}
  
void
Zeno::getInteriorSamplingHitPoints
(std::vector<Vector3<double> > const * * points) const {

  if (resultsInterior != nullptr) {
    resultsInterior->gatherHitPoints();

    if (points != nullptr) {
      *points = resultsInterior->getPoints();
    }
  }
}

double
Zeno::getInitializeTime() const {
  return initializeTimer.getTime();
}

double
Zeno::getPreprocessTime() const {
  return preprocessTimer.getTime();
}

double
Zeno::getWalkOnSpheresTime() const {
  return walkOnSpheresTimer.getTime();
}

double
Zeno::getWalkOnSpheresReductionTime() const {
  return walkOnSpheresReductionTimer.getTime();
}

double
Zeno::getInteriorSamplingTime() const {
  return interiorSamplingTimer.getTime();
}

double
Zeno::getInteriorSamplingReductionTime() const {
  return interiorSamplingReductionTimer.getTime();
}

double
Zeno::getVirialTime() const {
  return virialTimer.getTime();
}

double
Zeno::getVirialReductionTime() const {
  return virialReductionTimer.getTime();
}

double
Zeno::getTotalTime() const {
  return totalTimer.getTime();
}

long long 
Zeno::computeNumInProcess(int mpiSize, int mpiRank,
			  long long totalNumSamples) const {

  long long numSamplesInProcess = totalNumSamples / mpiSize;

  if (mpiRank < totalNumSamples % mpiSize) {
    numSamplesInProcess ++;
  }

  return numSamplesInProcess;
}

void
Zeno::computeDefaultParameters(ParametersWalkOnSpheres * parameters) const {

  // The default skin thickness is set to this factor times the model bounding
  // sphere radius.  This sets the skin thickness based on the problem scale.
  const double defaultSkinThicknessFactor = 0.000001;

  // The default launch radius is set to the model bounding sphere radius plus
  // this factor times the skin thickness.  This prevents potential numerical
  // problems by preventing the launch sphere from being too close to the model
  // or the skin.
  const double defaultLaunchRadiusFactor = 2;

  if (!parameters->getSkinThicknessWasSet()) {
    parameters->setSkinThickness(modelBoundingSphere.getRadius() * 
				 defaultSkinThicknessFactor);
  }
  
  if (!parameters->getLaunchCenterWasSet()) {
    parameters->setLaunchCenter(modelBoundingSphere.getCenter());
  }
  
  if (!parameters->getLaunchRadiusWasSet()) {
    parameters->setLaunchRadius(modelBoundingSphere.getRadius() +
				defaultLaunchRadiusFactor *
				parameters->getSkinThickness());
  }

  if (!parameters->getLaunchSphere().contains(modelBoundingSphere)) {
    std::cerr << std::endl
	      << "*** Warning ***" << std::endl
	      << "User-specified launch sphere may not contain the entire model"
	      << std::endl;
  }
}

void
Zeno::computeDefaultParameters(ParametersInteriorSampling * parameters) const {

  if (!parameters->getLaunchCenterWasSet()) {
    parameters->setLaunchCenter(modelBoundingSphere.getCenter());
  }
  
  if (!parameters->getLaunchRadiusWasSet()) {
    parameters->setLaunchRadius(modelBoundingSphere.getRadius());
  }

  if (!parameters->getLaunchSphere().contains(modelBoundingSphere)) {
    std::cerr << std::endl
	      << "*** Warning ***" << std::endl
	      << "User-specified launch sphere may not contain the entire model"
	      << std::endl;
  }
}

void
Zeno::computeDefaultParameters(ParametersVirial * parameters) const {

}

void
Zeno::computeDefaultParameters(ParametersResults * parameters) const {

  if (!parameters->getLengthScaleWasSet()) {
    parameters->setLengthScale(1, Units::Length::L);
  }
}

void
Zeno::setupRNGs(int numThreads,
		int seed,
		std::vector<RandomNumberGenerator> * threadRNGs) const {

  int numStreams = numThreads * mpiSize;

  threadRNGs->reserve(numThreads);

  for (int threadNum = 0;
       threadNum < numThreads;
       threadNum++) {
    
    int streamNum = mpiRank * numThreads + threadNum;

    threadRNGs->emplace_back(streamNum, numStreams, seed);
  }
}

long long 
Zeno::estimateTotalNum(double requestedError,
		       long long numSoFar,
		       Uncertain<double> const & currentValue) {

  double requestedVariance = 
    std::pow((requestedError / 100) * currentValue.getMean(), 2);

  long long estimatedTotalNum = 
    ceil(currentValue.getVariance() * numSoFar / requestedVariance);

  return estimatedTotalNum;
}

void
Zeno::getWalkOnSpheresResults
(long long numWalksInProcess,			
 ParametersWalkOnSpheres const & parametersWalkOnSpheres,
 ParametersResults const & parametersResults,
 BoundingSphere const & boundingSphere,
 Model const & model,
 std::vector<RandomNumberGenerator> * threadRNGs,
 ResultsZeno * * resultsZeno) {

  bool saveHitPoints = parametersWalkOnSpheres.getSaveSurfacePoints();

  if (parametersWalkOnSpheres.getTotalNumWalksWasSet() &&
      !parametersWalkOnSpheres.getMaxErrorCapacitanceWasSet() &&
      !parametersWalkOnSpheres.getMaxErrorPolarizabilityWasSet()) {

    *resultsZeno =
      new ResultsZeno(boundingSphere,
		      parametersWalkOnSpheres.getNumThreads(),
		      saveHitPoints);

    doWalkOnSpheres(parametersWalkOnSpheres,
		    numWalksInProcess,
		    boundingSphere, 
		    model,
		    threadRNGs,
		    *resultsZeno);

    walkOnSpheresReductionTimer.start();
    (*resultsZeno)->reduce();
    walkOnSpheresReductionTimer.stop();
  }
  else if (parametersWalkOnSpheres.getMaxErrorCapacitanceWasSet() ||
	   parametersWalkOnSpheres.getMaxErrorPolarizabilityWasSet()) {

    *resultsZeno =
      new ResultsZeno(boundingSphere,
		      parametersWalkOnSpheres.getNumThreads(),
		      saveHitPoints);

    ResultsCompiler resultsCompiler(parametersResults);

    Results results;

    long long estimatedNumWalksRemaining =
      parametersWalkOnSpheres.getMinTotalNumWalks();

    while (estimatedNumWalksRemaining > 0) {

      long long estimatedNumWalksRemainingInProcess = 
	computeNumInProcess(mpiSize, 
			    mpiRank, 
			    estimatedNumWalksRemaining);

      doWalkOnSpheres(parametersWalkOnSpheres,
		      estimatedNumWalksRemainingInProcess,
		      boundingSphere,
		      model,
		      threadRNGs,
		      *resultsZeno);

      walkOnSpheresReductionTimer.start();
      (*resultsZeno)->reduce();
      walkOnSpheresReductionTimer.stop();

      resultsCompiler.compile(*resultsZeno,
			      NULL,
			      NULL,
			      false,
			      &results);

      long long estimatedTotalNumWalks = 0;

      long long capacitanceEstimatedTotalNumWalks = 
  	estimateTotalNum(parametersWalkOnSpheres.getMaxErrorCapacitance(),
  			 (*resultsZeno)->getNumWalks(),
  			 results.capacitance.value);

      estimatedTotalNumWalks = std::max(estimatedTotalNumWalks,
  					capacitanceEstimatedTotalNumWalks);

      long long polarizabilityEstimatedTotalNumWalks = 
  	estimateTotalNum(parametersWalkOnSpheres.getMaxErrorPolarizability(),
  			 (*resultsZeno)->getNumWalks(),
  			 results.meanPolarizability.value);

      estimatedTotalNumWalks = std::max(estimatedTotalNumWalks,
  					polarizabilityEstimatedTotalNumWalks);

      if (parametersWalkOnSpheres.getTotalNumWalksWasSet()) {
	estimatedTotalNumWalks =
	  std::min(estimatedTotalNumWalks,
		   parametersWalkOnSpheres.getTotalNumWalks());
      }

      estimatedNumWalksRemaining = 
	estimatedTotalNumWalks - (*resultsZeno)->getNumWalks();

      // To avoid repeatedly undershooting, don't take a smaller step than
      // MinTotalNumWalks unless the iterations are about to end
      if (estimatedNumWalksRemaining > 0 &&
	  !(parametersWalkOnSpheres.getTotalNumWalksWasSet() &&
	    estimatedTotalNumWalks ==
	    parametersWalkOnSpheres.getTotalNumWalks())) {
	
	estimatedNumWalksRemaining = 
	  std::max(estimatedNumWalksRemaining,
		   parametersWalkOnSpheres.getMinTotalNumWalks());
      }

      if (parametersWalkOnSpheres.getMaxRunTimeWasSet() &&
	  (totalTimer.getTime() > parametersWalkOnSpheres.getMaxRunTime())) {

        break;
      }
    }
  }
}

void
Zeno::doWalkOnSpheres(ParametersWalkOnSpheres const & parameters,
		      long long numWalksInProcess,
		      BoundingSphere const & boundingSphere, 
		      Model const & nearestSurfacePointFinder,
		      std::vector<RandomNumberGenerator> * threadRNGs,
		      ResultsZeno * resultsZeno) {

  walkOnSpheresTimer.start();

  const int numThreads = parameters.getNumThreads();

  std::thread * * threads = new std::thread *[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {

    long long numWalksInThread = numWalksInProcess / numThreads;

    if (threadNum < numWalksInProcess % numThreads) {
      numWalksInThread ++;
    }

    threads[threadNum] = 
      new std::thread(doWalkOnSpheresThread,
		      &parameters,
		      &boundingSphere, 
		      &nearestSurfacePointFinder,
		      threadNum,
		      numWalksInThread,
		      &totalTimer,
		      &(threadRNGs->at(threadNum)),
		      resultsZeno);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    threads[threadNum]->join();
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    delete threads[threadNum];
  }

  delete [] threads;

  walkOnSpheresTimer.stop();
}

void
Zeno::doWalkOnSpheresThread(ParametersWalkOnSpheres const * parameters,
			    BoundingSphere const * boundingSphere, 
			    Model const * nearestSurfacePointFinder,
			    int threadNum,
			    long long numWalks,
			    Timer const * totalTimer,
			    RandomNumberGenerator * randomNumberGenerator,
			    ResultsZeno * resultsZeno) {

  const double shellThickness = parameters->getSkinThickness();
  
  WalkerExterior<double, 
		 RandomNumberGenerator,
		 Model,
		 RandomSpherePointGenerator,
		 BiasedSpherePointGenerator>
    walker(randomNumberGenerator, 
	   *boundingSphere, 
           *nearestSurfacePointFinder,
	   shellThickness);

  for (long long walkNum = 0; walkNum < numWalks; walkNum++) {

    bool hitObject = false;
    int numSteps   = 0;

    Vector3<double> startPoint;
    Vector3<double> endPoint;

    walker.walk(&hitObject, &numSteps,
		&startPoint, &endPoint);

    if (hitObject) {
      resultsZeno->recordHit(threadNum, 
			     startPoint, endPoint,
			     randomNumberGenerator);
    }
    else {
      resultsZeno->recordMiss(threadNum);
    }

    if (parameters->getMaxRunTimeWasSet() &&
	(totalTimer->getTime() > parameters->getMaxRunTime())) {

      break;
    }
  }
}

void
Zeno::getInteriorResults
(long long numSamplesInProcess,			
 ParametersInteriorSampling const & parametersInteriorSampling,
 ParametersResults const & parametersResults,
 BoundingSphere const & boundingSphere,
 Model const & model,
 std::vector<RandomNumberGenerator> * threadRNGs,
 ResultsInterior * * resultsInterior) {

  bool saveInteriorPoints = 
    parametersInteriorSampling.getSaveInteriorPoints() ||
    false;

  if (parametersInteriorSampling.getTotalNumSamplesWasSet() &&
      !parametersInteriorSampling.getMaxErrorVolumeWasSet()) {

    *resultsInterior =
      new ResultsInterior(boundingSphere,
			  parametersInteriorSampling.getNumThreads(),
			  saveInteriorPoints);

    doInteriorSampling(parametersInteriorSampling,
		       numSamplesInProcess,
		       boundingSphere,
		       model,
		       threadRNGs,
		       *resultsInterior);

    interiorSamplingReductionTimer.start();
    (*resultsInterior)->reduce();
    interiorSamplingReductionTimer.stop();
  }
  else if (parametersInteriorSampling.getMaxErrorVolumeWasSet()) {

    *resultsInterior =
      new ResultsInterior(boundingSphere,
			  parametersInteriorSampling.getNumThreads(),
			  saveInteriorPoints);

    ResultsCompiler resultsCompiler(parametersResults);

    Results results;

    long long estimatedNumSamplesRemaining =
      parametersInteriorSampling.getMinTotalNumSamples();

    while (estimatedNumSamplesRemaining > 0) {

      long long estimatedNumSamplesRemainingInProcess = 
	computeNumInProcess(mpiSize, 
			    mpiRank, 
			    estimatedNumSamplesRemaining);

      doInteriorSampling(parametersInteriorSampling,
			 estimatedNumSamplesRemainingInProcess,
			 boundingSphere, 
			 model,
			 threadRNGs,
			 *resultsInterior);

      interiorSamplingReductionTimer.start();
      (*resultsInterior)->reduce();
      interiorSamplingReductionTimer.stop();

      resultsCompiler.compile(NULL,
			      *resultsInterior,
			      NULL,
			      false,
			      &results);

      long long estimatedTotalNumSamples = 
	estimateTotalNum(parametersInteriorSampling.getMaxErrorVolume(),
			 (*resultsInterior)->getNumSamples(),
			 results.volume.value);

      if (parametersInteriorSampling.getTotalNumSamplesWasSet()) {
	estimatedTotalNumSamples =
	  std::min(estimatedTotalNumSamples,
		   parametersInteriorSampling.getTotalNumSamples());
      }

      estimatedNumSamplesRemaining = 
	estimatedTotalNumSamples - (*resultsInterior)->getNumSamples();

      // To avoid repeatedly undershooting, don't take a smaller step than
      // MinTotalNumSamples unless the iterations are about to end 
      if (estimatedNumSamplesRemaining > 0 &&
	  !(parametersInteriorSampling.getTotalNumSamplesWasSet() &&
	    estimatedTotalNumSamples ==
	    parametersInteriorSampling.getTotalNumSamples())) {
	
	estimatedNumSamplesRemaining = 
	  std::max(estimatedNumSamplesRemaining,
		   parametersInteriorSampling.getMinTotalNumSamples());
      }

      if (parametersInteriorSampling.getMaxRunTimeWasSet() &&
	  (totalTimer.getTime() > parametersInteriorSampling.getMaxRunTime())) {

        break;
      }
    }
  }

  if (false) {
    interiorSamplingReductionTimer.start();
    (*resultsInterior)->gatherHitPoints();
    interiorSamplingReductionTimer.stop();
  }
}

void
Zeno::doInteriorSampling(ParametersInteriorSampling const & parameters,
			 long long numSamplesInProcess,
			 BoundingSphere const & boundingSphere, 
			 Model const & insideOutsideTester,
			 std::vector<RandomNumberGenerator> * threadRNGs,
			 ResultsInterior * resultsInterior) {

  interiorSamplingTimer.start();

  const int numThreads = parameters.getNumThreads();

  std::thread * * threads = new std::thread *[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {

    long long numSamplesInThread = numSamplesInProcess / numThreads;

    if (threadNum < numSamplesInProcess % numThreads) {
      numSamplesInThread ++;
    }

    threads[threadNum] = 
      new std::thread(doInteriorSamplingThread,
		      &parameters,
		      &boundingSphere, 
		      &insideOutsideTester,
		      threadNum,
		      numSamplesInThread,
		      &totalTimer,
		      &(threadRNGs->at(threadNum)),
		      resultsInterior);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    threads[threadNum]->join();
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    delete threads[threadNum];
  }

  delete [] threads;

  interiorSamplingTimer.stop();
}

void
Zeno::doInteriorSamplingThread(ParametersInteriorSampling const * parameters,
			       BoundingSphere const * boundingSphere, 
			       Model const * insideOutsideTester,
			       int threadNum,
			       long long numSamples,
			       Timer const * totalTimer,
			       RandomNumberGenerator * randomNumberGenerator,
			       ResultsInterior * resultsInterior) {

  SamplerInterior<double, 
		 RandomNumberGenerator,
		 Model,
		 RandomBallPointGenerator>
    sampler(randomNumberGenerator, 
	    *boundingSphere, 
	    *insideOutsideTester);

  for (long long sampleNum = 0; sampleNum < numSamples; sampleNum++) {

    bool hitObject = false;

    Vector3<double> hitPoint;

    sampler.sample(&hitObject,
		   &hitPoint);

    if (hitObject) {
      resultsInterior->recordHit(threadNum,
				 hitPoint);
    }
    else {
      resultsInterior->recordMiss(threadNum);
    }

    if (parameters->getMaxRunTimeWasSet() &&
	(totalTimer->getTime() > parameters->getMaxRunTime())) {

      break;
    }
  }
}

void
Zeno::getVirialResults
(long long numStepsInProcess,
 ParametersVirial const & parametersVirial,
 std::vector<BoundingSphere> & boundingSpheres,
 Potential<double> const & potential,
 std::vector<RandomNumberGenerator> * threadRNGs,
 ResultsVirial * * resultsVirial) {

  double refDiameter = parametersVirial.getReferenceDiameter();
  int nFactorial = 1;
  for (int i=2; i<=parametersVirial.getOrder(); i++) nFactorial *= i;
  double refIntegral = nFactorial*std::pow(4.0*M_PI*refDiameter*refDiameter*refDiameter/3.0,parametersVirial.getOrder()-1)/2;
  int nVirialValues = potential.getFlexible() ?
       (parametersVirial.getOrder() == 2 ? 1 : 2) :
       (parametersVirial.getNumDerivatives()+1);
  *resultsVirial = new ResultsVirial(parametersVirial.getNumThreads(),
                                     nVirialValues,
                                     refIntegral);
  (*resultsVirial)->setOrder(parametersVirial.getOrder());

  doVirialSampling(parametersVirial,
                   numStepsInProcess,
                   boundingSpheres,
                   potential,
                   threadRNGs,
                   *resultsVirial,
                   refDiameter);

  virialReductionTimer.start();
  (*resultsVirial)->reduce();
  virialReductionTimer.stop();
}

/// Launches a set of virial-coefficient samples in each of a set of parallel
/// threads.
///
void
Zeno::doVirialSampling(ParametersVirial const & parameters,
                 long long stepsInProcess,
                 std::vector<BoundingSphere> & boundingSpheres,
                 Potential<double> const & potential,
                 std::vector<RandomNumberGenerator> * threadRNGs,
                 ResultsVirial * resultsVirial,
                 double refDiameter) {

    virialTimer.start();

    const int numThreads = parameters.getNumThreads();
    //std::cout<<numThreads<<std::endl;
    std::thread * * threads = new std::thread *[numThreads];

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {

        long long stepsInThread = stepsInProcess / numThreads;

        if (threadNum < stepsInProcess % numThreads) {
            stepsInThread ++;
        }

        threads[threadNum] =
                new std::thread(doVirialSamplingThread,
                                &parameters,
                                &boundingSpheres,
                                &allModels,
                                potential,
                                threadNum,
                                stepsInThread,
                                &totalTimer,
                                &(threadRNGs->at(threadNum)),
                                resultsVirial,
                                refDiameter);
    }

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        threads[threadNum]->join();
    }

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        delete threads[threadNum];
    }

    delete [] threads;

    virialTimer.stop();
}

void
Zeno::doVirialSamplingThread(ParametersVirial const * parameters,
			                 std::vector<BoundingSphere> * boundingSpheres,
			                 std::vector<Model> * models,
			                 Potential<double> const & potential,
			                 int threadNum,
			                 long long stepsInThread,
			                 Timer const * totalTimer,
			                 RandomNumberGenerator * randomNumberGenerator,
			                 ResultsVirial * resultsVirial,
			                 double refDiameter) {

    std::vector <int> numParticles;
    bool needRotate = false;
    if ((int)models->size() == parameters->getOrder()) {
      // one of each
      for (int i=0; i<(int)models->size(); i++) {
        numParticles.push_back(1);
        needRotate = needRotate || models->at(i).getSpheres()->size() > 1;
      }
    }
    else if (models->size() == 1) {
      needRotate = models->at(0).getSpheres()->size() > 1;
      numParticles.push_back(parameters->getOrder());
    }
    else {
      std::cerr << "Number of SPECIES must match virial coefficient order" << std::endl;
      exit(1);
    }
    double temperature = parameters->getTemperature();
    int numDerivatives = parameters->getNumDerivatives();
    IntegratorMSMC<double, RandomNumberGenerator> refIntegrator(threadNum,
                                                                 totalTimer,
                                                                 randomNumberGenerator,
                                                                 boundingSpheres,
                                                                 numParticles,
                                                                 models);

    ClusterSumChain<double> clusterSumRef(refIntegrator.getParticles(), refDiameter, 0.0, 1.0);
    ClusterSumWheatleyRecursion<double> clusterSumTargetRecursion(refIntegrator.getParticles(), &potential, temperature, 0);
    ClusterSumFlexible<double> clusterSumTargetFlex(refIntegrator.getParticles(), &potential, temperature);
    MCMoveChainVirial<double, RandomNumberGenerator> mcMoveChain(refIntegrator, &clusterSumRef, refDiameter);
    MCMoveRotate<double , RandomNumberGenerator> mcMoveRotateRef(refIntegrator, &clusterSumRef);
    MCMoveBondStretch<double , RandomNumberGenerator> mcMoveStretchRef(refIntegrator, &clusterSumRef, potential, temperature);
    MCMoveBondAngle<double , RandomNumberGenerator> mcMoveAngleRef(refIntegrator, &clusterSumRef, potential, temperature);
    MCMoveBondTorsion<double , RandomNumberGenerator> mcMoveTorsionRef(refIntegrator, &clusterSumRef, potential, temperature);
    refIntegrator.addMove(&mcMoveChain, 1.0);
    bool flexModel = potential.getFlexible();
    if (needRotate) {
      refIntegrator.addMove(&mcMoveRotateRef, 1.0);
      if (potential.getBondStyle() != Fixed) refIntegrator.addMove(&mcMoveStretchRef, 1.0);
      if (potential.getAngleStyle() != AngleFixed) refIntegrator.addMove(&mcMoveAngleRef, 1.0);
      if (potential.getAnyTorsion() && potential.getAngleStyle() != AngleNone) refIntegrator.addMove(&mcMoveTorsionRef, 1.0);
    }
    refIntegrator.setCurrentValue(clusterSumRef.getValues());
    ClusterSum<double> * clusterSumTarget = flexModel ? (ClusterSum<double>*)&clusterSumTargetFlex : (ClusterSum<double>*)&clusterSumTargetRecursion;

    IntegratorMSMC<double, RandomNumberGenerator> targetIntegrator(threadNum,
                                                                totalTimer,
                                                                randomNumberGenerator,
                                                                boundingSpheres,
                                                                numParticles,
                                                                models);

    ClusterSumChain<double> clusterSumRefT(targetIntegrator.getParticles(), refDiameter, 0.0, 1.0);
    ClusterSumWheatleyRecursion<double> clusterSumTargetRecursionT(targetIntegrator.getParticles(), &potential, temperature, numDerivatives);
    ClusterSumFlexible<double> clusterSumTargetFlexT(targetIntegrator.getParticles(), &potential, temperature);
    ClusterSum<double> * clusterSumTargetT = flexModel ? (ClusterSum<double>*)&clusterSumTargetFlexT : (ClusterSum<double>*)&clusterSumTargetRecursionT;
    MCMoveTranslate<double, RandomNumberGenerator> mcMoveTranslate(targetIntegrator, clusterSumTargetT);
    MCMoveRotate<double , RandomNumberGenerator> mcMoveRotateTarget(targetIntegrator, clusterSumTargetT);
    MCMoveBondStretch<double , RandomNumberGenerator> mcMoveStretchTarget(targetIntegrator, clusterSumTargetT, potential, temperature);
    MCMoveBondAngle<double , RandomNumberGenerator> mcMoveAngleTarget(targetIntegrator, clusterSumTargetT, potential, temperature);
    MCMoveBondTorsion<double , RandomNumberGenerator> mcMoveTorsionTarget(targetIntegrator, clusterSumTargetT, potential, temperature);
    targetIntegrator.addMove(&mcMoveTranslate, 1.0);
    if (needRotate) {
      targetIntegrator.addMove(&mcMoveRotateTarget, 1.0);
      if (potential.getBondStyle() != Fixed) targetIntegrator.addMove(&mcMoveStretchTarget, 1.0);
      if (potential.getAngleStyle() != AngleFixed) targetIntegrator.addMove(&mcMoveAngleTarget, 1.0);
      if (potential.getAnyTorsion() && potential.getAngleStyle() != AngleNone) targetIntegrator.addMove(&mcMoveTorsionTarget, 1.0);
    }
    targetIntegrator.setCurrentValue(clusterSumTargetT->getValues());

    VirialAlpha<double,RandomNumberGenerator> virialAlpha(refIntegrator, targetIntegrator,
            clusterSumRef, *clusterSumTarget, clusterSumRefT, *clusterSumTargetT);
    virialAlpha.run();

    double* alphaStats = virialAlpha.getAlphaStatistics();
    /*printf("alpha: %e  %e\n", alphaStats[0], alphaStats[1]);
    printf("alpha block correlation: %f\n", alphaStats[2]);
    printf("alpha span: %f\n", alphaStats[3]);*/

    VirialProduction<double, RandomNumberGenerator> virialProduction(refIntegrator,targetIntegrator,
            clusterSumRef, *clusterSumTarget, clusterSumRefT, *clusterSumTargetT, alphaStats[0],
            resultsVirial->getRefIntegral());
    virialProduction.getRefMeter()->setBlockSize(1);
    virialProduction.getTargetMeter()->setBlockSize(std::max(stepsInThread / 1000, 1LL));
    virialProduction.runSteps(stepsInThread);
    //virialProduction.printResults(NULL);

    resultsVirial->putData(threadNum, virialProduction.getRefMeter(), virialProduction.getTargetMeter());
    resultsVirial->putOverlapRatio(threadNum, virialProduction.getAlphaStats()[0], std::pow(virialProduction.getAlphaStats()[1],2));
    for (int iVal = 0; iVal < clusterSumTargetT->numValues()-1; iVal++) {
        resultsVirial->putVirialCoefficient(threadNum, iVal, virialProduction.getFullStats()[iVal][0], std::pow(virialProduction.getFullStats()[iVal][1],2));
    }
}
