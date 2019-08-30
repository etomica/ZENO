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

#include "OverlapTester.h"

///Creates a class to check if two particles are overalapped.
///

template <class T>
OverlapTester<T>::
OverlapTester(){
}

template <class T>
OverlapTester<T>::
  ~OverlapTester() {
}

///Checks if two particles obtained as parameters are overlapped.
///
template <class T>
bool
OverlapTester<T>::
isOverlapped(Particle<T> * a, Particle<T> * b) const {
    Vector3<T> x = a->getBoundingSphere()->getCenter() + a->getCenter();
    Vector3<T> y = b->getBoundingSphere()->getCenter() + b->getCenter();
    Vector3<T> distCenterVec = x - y;
    T distCenterSqr = distCenterVec.getMagnitudeSqr();
    T radiusX =  a->getBoundingSphere()->getRadius();
    T radiusY =  b->getBoundingSphere()->getRadius();
    if(distCenterSqr > ((radiusX + radiusY)* (radiusX + radiusY)))
    {
        return false;
    }

    for(int i = 0; i < a->numSpheres(); ++i)
    {
        Vector3<T> x = a->setFromSpherePosition(i);
        T radiusX = a->getModel()->getSpheres()->at(i).getRadius();
        for(int j = 0; j < b->numSpheres(); ++j)
        {
            Vector3<T> y = a->setFromSpherePosition(j);
            Vector3<T> distCenterVec = x - y;
            T distCenterSqr = distCenterVec.getMagnitudeSqr();
            T radiusY = b->getModel()->getSpheres()->at(i).getRadius();
            if(distCenterSqr < ((radiusX + radiusY)* (radiusX + radiusY)))
            {
                return true;
            }
        }
    }
    return false;
}

