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


// Smooth vertex attributes using uniform Laplacian
// Inputs:
//   Ain  #V by #A eigen Matrix of mesh vertex attributes (each vertex has #A attributes)
//   F    #F by 3 eigne Matrix of face (triangle) indices
// Output:
//   Aout #V by #A eigen Matrix of mesh vertex attributes
using namespace std;


bool make_solid(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::VectorXi &borderLoop,
	Eigen::MatrixXd & V_out,
	Eigen::MatrixXi & F_out,
	const double smoothf,
	const double offset);



#ifndef IGL_STATIC_LIBRARY
#  include "thicken.cpp"
#endif

#endif
