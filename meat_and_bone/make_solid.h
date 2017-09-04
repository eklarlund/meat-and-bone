// This file is part of libigl_stuff, sample libigl routines for constructing solids from surfaces.
// 
// Copyright (C) 2016 Esben Klarlund <eklarlund@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LIBIGL_STUFF_THICKEN
#define LIBIGL_STUFF_THICKEN
#include <Eigen/Core>
#include <vector>

// Turns a surface into a 3d solid by creating a second surface below the original by putting a vertex a a distance (offset) below each
// vertex in the original surface, then attaching the two surfaces at the border
//
// Inputs:
//   V  #V by dim matrix of vertex coordinates
//   F  #F by simplex_size  matrix of indices of simplex corners into V
//	 borderLoop  #N vector of indexes in F on the border, ordered in a loop
//	 smoothf  double that determines how much the normals will be smoothed before making secind surface
//	 offset  distance between old and new surface
//	 
// 	 offset Thickness of solid
// Output:
//   V_out #V by dim matrix of vertex coordinates of solid
//	 F  #F by simplex_size  matrix of indices of simplex corners into V od solid
//
using namespace std;

template <typename DerivedV, typename DerivedF, typename DerivedBL, typename DerivedV_out, typename DerivedF_out>

bool make_solid(
	const Eigen::PlainObjectBase<DerivedV> &V,
	const Eigen::PlainObjectBase<DerivedF> &F,
	const Eigen::PlainObjectBase<DerivedBL> &borderLoop,
	Eigen::PlainObjectBase<DerivedV_out> & V_out,
	Eigen::PlainObjectBase<DerivedF_out> & F_out,
	const double smoothf,
	const double offset);



#ifndef IGL_STATIC_LIBRARY
#  include "make_solid.cpp"
#endif

#endif
