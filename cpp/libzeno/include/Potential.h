// ================================================================
//
// This software is not subject to copyright protection in the United States
// and is considered to be in the public domain. Permission to freely use,
// copy, modify, and distribute this software and its documentation without fee
// is hereby granted, provided that this notice and disclaimer of warranty
// appears in all copies.
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
// Authors: Andrew Schultz <ajs42@buffalo.edu>
// Created: Wed Feb 22 2023
//
// ================================================================

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <map>
#include <set>
#include <string>
#include <limits>
#include <iostream>
#include <cmath>

#include "Virials/Particle.h"

// ================================================================

namespace zeno {

struct BondedPair {
  int bondType;
  int sphere1;
  int sphere2;
};

struct AngleTriplet {
  int angleType;
  int sphere1;
  int sphere2; // central sphere
  int sphere3;
};

enum BondStyle {Fixed, Harmonic, FENE};

enum AngleStyle {AngleNone, AngleFixed, AngleHarmonic};

enum NonbondStyle {HardSphere, LennardJones, WCA};

struct BondedPartner {
  int sphere2;
  int bondType;
};

template <class T>
class Potential {
 public:
  Potential();
  Potential(const Potential &p);
  ~Potential();

  void setBondStyleName(std::string bondStyle);
  void addBondPair(int bondType, int sphere1, int sphere2);
  void setBondCoeff(int bondType, double* coeffs);

  void setAngleStyleName(std::string angleStyle);
  // sphere2 is the central sphere
  void addAngleTriplet(int angleType, int sphere1, int sphere2, int sphere3);
  void setAngleCoeff(int angleType, double* coeffs);

  void setNonbondStyleName(std::string nonbondStyle);
  void setNonbondCoeff(int sphereType, double* coeffs);

  const BondStyle getBondStyle() const;
  std::vector<BondedPair> getBondPairs() const;
  const std::vector<double*> * getBondCoeffs() const;

  const AngleStyle getAngleStyle() const;
  std::vector<AngleTriplet> getAngleTriplets() const;
  const std::vector<double*> * getAngleCoeffs() const;

  const NonbondStyle getNonbondStyle() const;
  const std::vector<double*> * getNonbondCoeffs1() const;

  void initialize(int numSpheres);

  double energy1(Particle<T> * particle) const;
  double energy2(Particle<T> * particle1, Particle<T> * particle2) const;

  const std::vector<std::vector<BondedPartner>> * getBondedPartners() const;

  const bool getHasTorsion() const;
  const bool getFlexible() const;

 private:
  bool empty;
  bool hasTorsion;

  BondStyle bondStyle;
  AngleStyle angleStyle;
  NonbondStyle nonbondStyle;

  std::vector<BondedPair> bondPairs;
  std::vector<double*> bondCoeffs;
  std::vector<std::vector<BondedPartner>> bondedPartners;

  std::vector<AngleTriplet> angleTriplets;
  std::vector<double*> angleCoeffs;

  std::vector<double*> nonbondCoeffs1;
  std::vector<std::vector<double*>> nonbondCoeffs;

  double** nonbondedScaling;

  double uHarmonic(Vector3<T> ri, Vector3<T> rj, double* c) const;
  double uFENE(Vector3<T> ri, Vector3<T> rj, double* c) const;

  double uAngleHarmonic(Vector3<T> ri, Vector3<T> rj, Vector3<T> rk, double* c) const;

  double uHardSphere1(Vector3<T> ri, Vector3<T> rj, double sigma) const;
  double uHardSphere(Vector3<T> ri, Vector3<T> rj, double* c) const;
  double uLennardJones(Vector3<T> ri, Vector3<T> rj, double* c) const;
  double uWCA(Vector3<T> ri, Vector3<T> rj, double* c) const;
};

template <class T>
Potential<T>::Potential()
: empty(true),
  hasTorsion(false),
  bondStyle(BondStyle::Fixed),
  angleStyle(AngleStyle::AngleFixed),
  nonbondStyle(NonbondStyle::HardSphere),
  nonbondedScaling(nullptr) {

}

template <class T>
Potential<T>::Potential(const Potential &p) : empty(true),
                                           hasTorsion(false),
                                           bondStyle(p.getBondStyle()),
                                           angleStyle(p.getAngleStyle()),
                                           nonbondStyle(p.getNonbondStyle()),
                                           bondPairs(p.getBondPairs()),
                                           angleTriplets(p.getAngleTriplets()) {
  for (int i=0; i<(int)p.getBondCoeffs()->size(); i++) {
    double* c = p.getBondCoeffs()->at(i);
    if (!c) {
      bondCoeffs.push_back(c);
    }
    else {
      if (bondStyle == BondStyle::Harmonic) {
        double* cNew = new double[2];
        cNew[0] = c[0];
        cNew[1] = c[1];
        bondCoeffs.push_back(cNew);
      }
      else if (bondStyle == BondStyle::FENE) {
        double* cNew = new double[4];
        cNew[0] = c[0];
        cNew[1] = c[1];
        cNew[2] = c[2];
        cNew[3] = c[3];
        bondCoeffs.push_back(cNew);
      }
    }
  }
  for (int i=0; i<(int)p.getAngleCoeffs()->size(); i++) {
    double* c = p.getAngleCoeffs()->at(i);
    if (!c) {
      angleCoeffs.push_back(c);
    }
    else {
      if (angleStyle == AngleStyle::AngleHarmonic) {
        double* cNew = new double[2];
        cNew[0] = c[0];
        cNew[1] = c[1];
        angleCoeffs.push_back(cNew);
      }
    }
  }

  for (int i=0; i<(int)p.getNonbondCoeffs1()->size(); i++) {
    double* c = p.getNonbondCoeffs1()->at(i);
    if (!c) {
      nonbondCoeffs1.push_back(c);
    }
    else {
      if (nonbondStyle == NonbondStyle::HardSphere) {
        double* cNew = new double[1];
        cNew[0] = c[0];
        nonbondCoeffs1.push_back(cNew);
      }
      else if (nonbondStyle == NonbondStyle::LennardJones) {
        double* cNew = new double[2];
        cNew[0] = c[0];
        cNew[1] = c[1];
        nonbondCoeffs1.push_back(cNew);
      }
      else if (nonbondStyle == NonbondStyle::WCA) {
        double* cNew = new double[2];
        cNew[0] = c[0];
        cNew[1] = c[1];
        nonbondCoeffs1.push_back(cNew);
      }
    }
  }

  initialize(p.bondedPartners.size());
}

template <class T>
Potential<T>::~Potential() {
  for (int i=0; i<(int)bondCoeffs.size(); i++) delete [] bondCoeffs[i];
  for (int i=0; i<(int)angleCoeffs.size(); i++) delete [] angleCoeffs[i];
  for (int i=0; i<(int)nonbondCoeffs.size(); i++) {
    for (int j=i; j<(int)nonbondCoeffs[i].size(); j++) delete [] nonbondCoeffs[i][j];
  }
  for (int i=0; i<(int)bondedPartners.size(); i++) delete [] nonbondedScaling[i];
  delete[] nonbondedScaling;
}

/// Set the bond style for the potential.
///
template <class T>
void
Potential<T>::setBondStyleName(std::string bondStyleName) {
  if (bondStyleName.compare("fixed") == 0) {
    bondStyle = BondStyle::Fixed;
  }
  else if (bondStyleName.compare("harmonic") == 0) {
    bondStyle = BondStyle::Harmonic;
  }
  else if (bondStyleName.compare("fene") == 0) {
    bondStyle = BondStyle::FENE;
  }
  else {
    std::cerr << "Unrecognized bond style " << bondStyleName << std::endl;

    exit(1);
  }
}

/// Add a bond to the potential.
///
template <class T>
void
Potential<T>::addBondPair(int bondType, int sphere1, int sphere2) {
  struct BondedPair bp = {bondType, sphere1-1, sphere2-1};
  bondPairs.push_back(bp);
}

/// Set the bond potential parameters for the given bond type.
///
template <class T>
void
Potential<T>::setBondCoeff(int bondType, double* coeffs) {
  if ((int)bondCoeffs.size() <= bondType) {
    bondCoeffs.resize(bondType+1, nullptr);
  }
  bondCoeffs[bondType] = coeffs;
}

/// Set the angle style for the potential.
///
template <class T>
void
Potential<T>::setAngleStyleName(std::string angleStyleName) {
  if (angleStyleName.compare("none") == 0) {
    angleStyle = AngleStyle::AngleNone;
  }
  else if (angleStyleName.compare("fixed") == 0) {
    angleStyle = AngleStyle::AngleFixed;
  }
  else if (angleStyleName.compare("harmonic") == 0) {
    angleStyle = AngleStyle::AngleHarmonic;
  }
  else {
    std::cerr << "Unrecognized angle style " << angleStyleName << std::endl;

    exit(1);
  }
}

/// Add a bond to the potential.
///
template <class T>
void
Potential<T>::addAngleTriplet(int angleType, int sphere1, int sphere2, int sphere3) {
  struct AngleTriplet bt = {angleType, sphere1-1, sphere2-1, sphere3-1};
  angleTriplets.push_back(bt);
}

/// Set the bond potential parameters for the given bond type.
///
template <class T>
void
Potential<T>::setAngleCoeff(int angleType, double* coeffs) {
  if ((int)angleCoeffs.size() <= angleType) {
    angleCoeffs.resize(angleType+1, nullptr);
  }
  angleCoeffs[angleType] = coeffs;
}

/// Set the bond potential parameters for the given bond type.
///
template <class T>
void
Potential<T>::setNonbondCoeff(int sphereType, double* coeffs) {
  if ((int)nonbondCoeffs1.size() <= sphereType) {
    nonbondCoeffs1.resize(sphereType+1, nullptr);
  }
  nonbondCoeffs1[sphereType] = coeffs;
}

/// Set the bond style for the model.
///
template <class T>
void
Potential<T>::setNonbondStyleName(std::string nonbondStyleName) {
  if (nonbondStyleName.compare("hs") == 0) {
    nonbondStyle = NonbondStyle::HardSphere;
  }
  else if (nonbondStyleName.compare("lj") == 0) {
    nonbondStyle = NonbondStyle::LennardJones;
  }
  else if (nonbondStyleName.compare("wca") == 0) {
    nonbondStyle = NonbondStyle::WCA;
  }
  else {
    std::cerr << "Unrecognized nonbond style " << nonbondStyleName << std::endl;

    exit(1);
  }
}

template <class T>
const BondStyle
Potential<T>::getBondStyle() const {
  return bondStyle;
}

template <class T>
const AngleStyle
Potential<T>::getAngleStyle() const {
  return angleStyle;
}

template <class T>
const NonbondStyle
Potential<T>::getNonbondStyle() const {
  return nonbondStyle;
}

template <class T>
std::vector<BondedPair>
Potential<T>::getBondPairs() const {
  return bondPairs;
}

template <class T>
const std::vector<double*> *
Potential<T>::getBondCoeffs() const {
  return &bondCoeffs;
}

template <class T>
std::vector<AngleTriplet>
Potential<T>::getAngleTriplets() const {
  return angleTriplets;
}

template <class T>
const std::vector<double*> *
Potential<T>::getAngleCoeffs() const {
  return &angleCoeffs;
}

template <class T>
const std::vector<double*> *
Potential<T>::getNonbondCoeffs1() const {
  return &nonbondCoeffs1;
}

template <class T>
void
Potential<T>::initialize(int numSpheres) {

  nonbondedScaling = new double*[numSpheres];
  for (int i=0; i<numSpheres; i++) {
    nonbondedScaling[i] = new double[numSpheres];
    for (int j=0; j<numSpheres; j++) {
      nonbondedScaling[i][j] = 1;
    }
  }

  bondedPartners.resize(numSpheres);

  for (BondedPair bp : bondPairs) {
    BondedPartner partner1 = {bp.sphere2, bp.bondType};
    bondedPartners[bp.sphere1].push_back(partner1);
    BondedPartner partner2 = {bp.sphere1, bp.bondType};
    bondedPartners[bp.sphere2].push_back(partner2);
    nonbondedScaling[bp.sphere1][bp.sphere2] = 0;
    nonbondedScaling[bp.sphere2][bp.sphere1] = 0;
  }

  for (AngleTriplet at : angleTriplets) {
    if (bondStyle == Fixed) {
      // need to use angle triplets to construct topology
      // these won't be used to compute energy, just to do MC moves
      BondedPartner partner2 = {at.sphere1, -1};
      bondedPartners[at.sphere2].push_back(partner2);
      partner2 = {at.sphere3, -1};
      bondedPartners[at.sphere2].push_back(partner2);
    }
    nonbondedScaling[at.sphere1][at.sphere2] = 0;
    nonbondedScaling[at.sphere1][at.sphere3] = 0;
    nonbondedScaling[at.sphere2][at.sphere3] = 0;
    nonbondedScaling[at.sphere2][at.sphere1] = 0;
    nonbondedScaling[at.sphere3][at.sphere1] = 0;
    nonbondedScaling[at.sphere3][at.sphere2] = 0;
  }

  if (bondStyle != Fixed || angleStyle != AngleFixed) {
    // check that our model is fully connected
    std::vector<int> stack;
    std::set<int> seen;
    if (bondedPartners[0].size() > 0) {
      seen.insert(0);
      for (BondedPartner bp0 : bondedPartners[0]) {
        seen.insert(bp0.sphere2);
        stack.push_back(bp0.sphere2);
      }

      while (stack.size() > 0) {
        int n = stack[stack.size()-1];
        stack.pop_back();
        for (BondedPartner bpn : bondedPartners[n]) {
          int s2 = bpn.sphere2;
          if (seen.count(s2) == 0) {
            seen.insert(s2);
            stack.push_back(s2);
          }
        }
      }
    }
    if ((int)seen.size() != numSpheres) {
      std::cerr << "Connectivity missing for some spheres.  From first sphere (1), cannot reach";
      for (int i=0; i<numSpheres; i++) {
        if (seen.count(i) == 0) std::cerr << " " << (i+1);
      }
      std::cerr << std::endl;
      exit(1);
    }
  }

  if (angleStyle != AngleFixed && numSpheres >= 4) {
    // can we find an atom with 2+ bond partners where 1 of the partners also has 2+ partners
    for (int a=0; a<numSpheres && !hasTorsion; a++) {
      if (bondedPartners[a].size() >= 2) {
        for (int bb=0; bb<(int)bondedPartners[a].size(); bb++) {
          int b = bondedPartners[a][bb].sphere2;
          if (bondedPartners[b].size() >= 2) {
            hasTorsion = true;
            break;
          }
        }
      }
    }
  }

  nonbondCoeffs.resize(nonbondCoeffs1.size());
  for (int i=0; i<(int)nonbondCoeffs1.size(); i++) {
    nonbondCoeffs[i].resize(nonbondCoeffs1.size());
  }
  for (int i=0; i<(int)nonbondCoeffs1.size(); i++) {
    double* inbc = nonbondCoeffs1[i];
    if (!inbc) {
      // there may be no spheres of this type, so ignore now and hope for the best later
      continue;
    }
    empty = false;
    for (int j=0; j<=i; j++) {
      double* jnbc = nonbondCoeffs1[j];
      if (!jnbc) continue;
      double* cNew = nullptr;
      if (nonbondStyle == NonbondStyle::HardSphere) {
        cNew = new double[1];
        // additive
        cNew[0] = (inbc[0] + jnbc[0]) / 2;
      }
      else if (nonbondStyle == NonbondStyle::LennardJones ||
               nonbondStyle == NonbondStyle::WCA) {
        cNew = new double[2];
        // Lorentz-Berthelot
        cNew[0] = (inbc[0] + jnbc[0]) / 2; // sigma
        cNew[1] = std::sqrt(inbc[1] * jnbc[1]); // epsilon
      }
      nonbondCoeffs[i][j] = nonbondCoeffs[j][i] = cNew;
    }
  }
}

template <class T>
double
Potential<T>::energy1(Particle<T> * particle) const {
  if ((bondStyle == Fixed && angleStyle == AngleFixed) || particle->numSpheres() == 1) return 0;
  double uTot = 0;
  for (int i=0; i<particle->numSpheres(); i++) {
    if (bondStyle != Fixed) {
      for (BondedPartner bp : bondedPartners[i]) {
        if (bp.sphere2 < i) continue;
        switch (bondStyle) {
          case Harmonic:
            uTot += uHarmonic(particle->getSpherePosition(i),
                              particle->getSpherePosition(bp.sphere2),
                              bondCoeffs[bp.bondType]);
            break;
          case FENE:
            uTot += uFENE(particle->getSpherePosition(i),
                          particle->getSpherePosition(bp.sphere2),
                          bondCoeffs[bp.bondType]);
            break;
          default:
            fprintf(stderr, "Unknown bond type\n");
            exit(1);
            break;
        }
      }
    }

    if (empty) {
      double ri = particle->getModel()->getSpheres()->at(i).getRadius();
      for (int j=i+1; j<particle->numSpheres(); j++) {
        // exclude bonded pairs
        int rj = particle->getModel()->getSpheres()->at(j).getRadius();

        uTot += uHardSphere1(particle->getSpherePosition(i),
                             particle->getSpherePosition(j),
                             ri+rj);
      }
    }
    else {

      int iType = (int)(particle->getModel()->getSpheres()->at(i).getRadius());
      for (int j=i+1; j<particle->numSpheres(); j++) {
        // exclude bonded pairs
        int jType = (int)(particle->getModel()->getSpheres()->at(j).getRadius());
        double ulj, nbs, uwca;
        switch (nonbondStyle) {
          case HardSphere:
            if (nonbondedScaling[i][j] == 0) break;
            uTot += uHardSphere(particle->getSpherePosition(i),
                                particle->getSpherePosition(j),
                                nonbondCoeffs[iType][jType]);
            break;
          case LennardJones:
            nbs = nonbondedScaling[i][j];
            if (nbs > 0) {
              ulj = uLennardJones(particle->getSpherePosition(i),
                                    particle->getSpherePosition(j),
                                    nonbondCoeffs[iType][jType]);
              uTot += nbs*ulj;
            }
            break;
          case WCA:
            nbs = nonbondedScaling[i][j];
            if (nbs > 0) {
              uwca = uWCA(particle->getSpherePosition(i),
                         particle->getSpherePosition(j),
                         nonbondCoeffs[iType][jType]);
              uTot += nbs*uwca;
            }
            break;
          default:
            fprintf(stderr, "Unknown nonbond type\n");
            exit(1);
            break;
        }
      }
      if (std::isinf(uTot)) return uTot;
    }
  }
  if (angleStyle != AngleFixed) {
    for (AngleTriplet at : angleTriplets) {
      switch (angleStyle) {
        case AngleHarmonic:
          uTot += uAngleHarmonic(particle->getSpherePosition(at.sphere1),
                                 particle->getSpherePosition(at.sphere2),
                                 particle->getSpherePosition(at.sphere3),
                                 angleCoeffs[at.angleType]);
          break;
        default:
          fprintf(stderr, "Unknown angle type\n");
          exit(1);
          break;
      }
    }
  }
  return uTot;
}

template <class T>
double
Potential<T>::energy2(Particle<T> * particle1, Particle <T> * particle2) const {
  if (bondStyle == Fixed && angleStyle == AngleFixed && nonbondStyle == HardSphere && empty) {
    Vector3<T> x = particle1->getBoundingSpherePosition();
    Vector3<T> y = particle2->getBoundingSpherePosition();
    Vector3<T> distCenterVec = x - y;
    T distCenterSqr = distCenterVec.getMagnitudeSqr();
    T radiusX =  particle1->getBoundingSphereRadius();
    T radiusY =  particle2->getBoundingSphereRadius();
    if(distCenterSqr > ((radiusX + radiusY) * (radiusX + radiusY))) {
        return 0;
    }
  }

  double uTot = 0;
  for (int i=0; i<particle1->numSpheres(); i++) {

    if (empty) {
      double ri = particle1->getModel()->getSpheres()->at(i).getRadius();
      for (int j=0; j<particle2->numSpheres(); j++) {
        // exclude bonded pairs
        double rj = particle2->getModel()->getSpheres()->at(j).getRadius();

        uTot += uHardSphere1(particle1->getSpherePosition(i),
                             particle2->getSpherePosition(j),
                             ri+rj);
      }
    }
    else {

      int iType = (int)(particle1->getModel()->getSpheres()->at(i).getRadius());
  
      for (int j=0; j<particle2->numSpheres(); j++) {
        // exclude bonded pairs
        int jType = (int)(particle2->getModel()->getSpheres()->at(j).getRadius());
        switch (nonbondStyle) {
          case HardSphere:
            uTot += uHardSphere(particle1->getSpherePosition(i),
                                particle2->getSpherePosition(j),
                                nonbondCoeffs[iType][jType]);
            break;
          case LennardJones:
            uTot += uLennardJones(particle1->getSpherePosition(i),
                                  particle2->getSpherePosition(j),
                                  nonbondCoeffs[iType][jType]);
            break;
          case WCA:
            uTot += uWCA(particle1->getSpherePosition(i),
                         particle2->getSpherePosition(j),
                         nonbondCoeffs[iType][jType]);
            break;
          default:
            fprintf(stderr, "Unknown nonbond type\n");
            exit(1);
            break;
        }
      }
      if (std::isinf(uTot)) return uTot;
    }
  }
  return uTot;
}

template <class T>
double
Potential<T>::uHarmonic(Vector3<T> ri, Vector3<T> rj, double* c) const {
  Vector3<T> dr = rj - ri;
  double dr1 = std::sqrt(dr.getMagnitudeSqr()) - c[0];
  return c[1]*dr1*dr1;
}

template <class T>
double
Potential<T>::uFENE(Vector3<T> ri, Vector3<T> rj, double* c) const {
  Vector3<T> dr = rj - ri;
  double r2 = dr.getMagnitudeSqr();
  double R0sq = c[1]*c[1];
  if (r2 > R0sq) return INFINITY;

  double s2 = c[2]*c[2]/r2;
  double s6 = s2*s2*s2;
  double s12 = s6*s6;

  return -0.5*c[0]*R0sq*std::log(1 - r2/R0sq) + 4*c[3]*(s12 - s6) + c[3];
}

template <class T>
double
Potential<T>::uAngleHarmonic(Vector3<T> ri, Vector3<T> rj, Vector3<T> rk, double* c) const {
  Vector3<T> drij = ri - rj;
  Vector3<T> drkj = rk - rj;
  double theta = std::acos(drij.dot(drkj) / std::sqrt(drij.getMagnitudeSqr() * drkj.getMagnitudeSqr()));
  double dtheta = theta - c[0];
  return c[1]*dtheta*dtheta;
}

template <class T>
double
Potential<T>::uHardSphere(Vector3<T> ri, Vector3<T> rj, double* c) const {
  Vector3<T> dr = rj - ri;
  double r2 = dr.getMagnitudeSqr();
  return r2 > c[0]*c[0] ? 0 : INFINITY;
}

template <class T>
double
Potential<T>::uHardSphere1(Vector3<T> ri, Vector3<T> rj, double sigma)  const {
  Vector3<T> dr = rj - ri;
  double r2 = dr.getMagnitudeSqr();
  return r2 > sigma*sigma ? 0 : INFINITY;
}

template <class T>
double
Potential<T>::uLennardJones(Vector3<T> ri, Vector3<T> rj, double* c) const {
  Vector3<T> dr = rj - ri;
  double s2 = c[0]*c[0]/dr.getMagnitudeSqr();
  if (s2 > 1e50) return INFINITY;
  double s6 = s2*s2*s2;
  double s12 = s6*s6;
  return 4*c[1]*(s12 - s6);
}

template <class T>
double
Potential<T>::uWCA(Vector3<T> ri, Vector3<T> rj, double* c) const {
  Vector3<T> dr = rj - ri;
  double s2 = c[0]*c[0]/dr.getMagnitudeSqr();
  if (s2 > 1e50) return INFINITY;
  if (s2 < 0.793700525984100) return 0;
  double s6 = s2*s2*s2;
  double s12 = s6*s6;
  return 4*c[1]*(s12 - s6) + c[1];
}

template <class T>
const std::vector<std::vector<BondedPartner>> *
Potential<T>::getBondedPartners() const {
   return &bondedPartners;
}

template <class T>
const bool
Potential<T>::getHasTorsion() const {
  return hasTorsion;
}

template <class T>
const bool
Potential<T>::getFlexible() const {
  return bondStyle != Fixed || angleStyle != AngleFixed || (hasTorsion && angleStyle != AngleFixed && angleStyle != AngleNone);
}

}
#endif
