// This file is part of libigl_stuff, sample libigl routines for constructing solids from surfaces.
// 
// Copyright (C) 2016 Esben Klarlund <eklarlund@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LIBIGL_STUFF_LAPLACE_SMOOTH
#define LIBIGL_STUFF_LAPLACE_SMOOTH
#include <Eigen/Core>


// Smooth vertex attributes using uniform Laplacian
// Inputs:
//   Ain  #V by #A eigen Matrix of mesh vertex attributes (each vertex has #A attributes)
//   F    #F by 3 eigne Matrix of face (triangle) indices
// Output:
//   Aout #V by #A eigen Matrix of mesh vertex attributes
template <typename DerivedV, typename DerivedF, typename DerivedVS, typename DerivedNS>
void laplace_smooth(
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedVS> & V_smooth,
	Eigen::PlainObjectBase<DerivedNS> & N,
	double smooth_factor);
//because this function uses Eigen::conservativeResize, along with libigl functions calling PlainObjectBase, it must use PlainObjectBase functions


#ifndef IGL_STATIC_LIBRARY
#  include "laplace_smooth.cpp"
#endif

#endif
