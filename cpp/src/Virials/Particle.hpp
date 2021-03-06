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

#include "Particle.h"

///  Constructs the class to define a particle as an assembly of spheres with mutable center and orientation.
///

template <class T>
Particle<T>::
Particle(MixedModel<T> const & model, Sphere<T> const & boundingSphere) : model(model), boundingSphere(boundingSphere)
              {
}

template <class T>
Particle<T>::
  ~Particle() {
}

/// Obtain number of spheres in the particle.
///
template <class T>
int
Particle<T>::
numSpheres(){
    return model.getSpheres() -> size();
}

/// Obtain center of particle. Note: This is a constant, hence cannot be modified directly.
///
template <class T>
const Vector3<T>
Particle<T>::
getCenter() const {
    return center;
}

/// Sets the center of particle equal to the vector3 obtained as a function parameter.
///
template <class T>
void
Particle<T>::
setCenter(Vector3<T> v) {
    center = v;
}

/// Translates the particle by a step obtained as a function parameter.
///
template <class T>
void
Particle<T>::
translateBy(Vector3<T> step){
    center += step;
}

/// Rotates the particle about an axis and an by an angle, both obtained as a function parameters.
///
template <class T>
void
Particle<T>::
rotateBy(Vector3<T> axis, T angle){
    Matrix3x3<T> rotation;
    rotation.setAxisAngle(axis, angle);
    rotation.transform(orientation);
}

/// Sets the sphere at index of particle based on new position of particle.
///
template <class T>
const Vector3<T>
Particle<T>::
getSpherePosition( int index) const {
    Vector3<T> position = model.getSpheres() -> at(index).getCenter();
    orientation.transform(position);
    position += center;
    return position;
}

/// Obtain the assembly of spheres which constitute a particle.
///
template <class T>
MixedModel<T> const *
Particle<T>::
getModel(){
    return &model;
}

/// Returns a bounding sphere around the assembly of spheres.
///
template <class T>
Sphere<T> const *
Particle<T>::
getBoundingSphere(){
    return &boundingSphere;
}

template <class T>
const Vector3<T>
Particle<T>::
getBoundingSpherePosition() const {
    Vector3<T> position = boundingSphere.getCenter();
    orientation.transform(position);
    position += center;
    return position;
}
