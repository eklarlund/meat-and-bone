// This file is part of libigl_stuff, sample libigl routines for constructing solids from surfaces.
// 
// Copyright (C) 2016 Esben Klarlund <eklarlund@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LIBIGL_STUFF_SKIRT
#define LIBIGL_STUFF_SKIRT
#include <Eigen/Core>
#include <vector>


// creates a skirt around surface
//
// Inputs:
//   V  #V by dim matrix of vertex coordinates
//   F  #F by simplex_size  matrix of indices of simplex corners into V
//	 borderNormals	#V by dim matrix of normals at each point in V
//	 borderLoop #N vector of indexes in F on the border, ordered in a loop
//	 skirts  number of skirts
//	 displacement  distance skirt reaches along the surface plane
//	 offset  distance skirt goes downwards
//	 drop  extra distance that the first skirt reaches
//	 
// 	 offset Thickness of solid
// Output:
//   V_plus #V by dim matrix of vertex coordinates of skirted surface
//	 F_plus  #F by simplex_size  matrix of indices of skirted surface
//

using namespace std;
template <typename DerivedV, typename DerivedF, typename DerivedBN, typename DerivedBB, typename DerivedV_plus, typename DerivedF_plus,
	typename DerivedBL>
bool make_skirt(
	const Eigen::PlainObjectBase<DerivedV> & V,
	const Eigen::PlainObjectBase<DerivedF> & F,
	const Eigen::PlainObjectBase<DerivedBN> & borderNormals,
	const Eigen::PlainObjectBase<DerivedBB> & borderBitangents,
	Eigen::PlainObjectBase<DerivedV_plus> & V_plus,
	Eigen::PlainObjectBase<DerivedF_plus> & F_plus,
	Eigen::PlainObjectBase<DerivedBL> & borderLoop,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop
);

template <typename DerivedV, typename DerivedF, typename DerivedN, typename DerivedV_plus, typename DerivedF_plus, typename DerivedBL>
bool make_skirt(
	const Eigen::PlainObjectBase<DerivedV> & V,
	const Eigen::PlainObjectBase<DerivedF> & F,
	const Eigen::PlainObjectBase<DerivedN> & N,
	Eigen::PlainObjectBase<DerivedV_plus> & V_plus,
	Eigen::PlainObjectBase<DerivedF_plus> & F_plus,
	Eigen::PlainObjectBase<DerivedBL> & borderLoop,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop
);



#endif
