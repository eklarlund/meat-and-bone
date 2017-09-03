// This file is part of libigl_stuff, sample libigl routines for constructing solids from surfaces.
// 
// Copyright (C) 2016 Esben Klarlund <eklarlund@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LIBIGL_STUFF_INTERPOLATE
#define LIBIGL_STUFF_INTERPOLATE
#include <Eigen/Core>
#include <vector>


// Uses ray shooting to find the distance between each point on surface A to surface B along the normals of surface A, if the normal does not intersect, then DOUBLE_MAX is returned
//
// Inputs:
//   V  #V by dim matrix of vertex coordinates in A
//	 N_smooth  #V by 3 matrix containing normal vectors of A
//   V_ref  #V by dim matrix of vertex coordinates in B
//   F_ref  #F by simplex_size  matrix of indices of simplex corners into V in B
//	 
// Output:
//   distances	#V dim vector of the distances from each point in a to B
//	 hits	the number of normals from A that intersect with B
//	 misses	 the number of normals from A that do not intersect with B
//
template <typename DerivedV, typename DerivedN, typename DerivedV_ref, typename DerivedF_ref, typename DerivedD>
void find_distance(
	const Eigen::MatrixBase<DerivedV> & V,
	const Eigen::MatrixBase<DerivedN> & N_smooth,
	const Eigen::MatrixBase<DerivedV_ref> & V_ref,
	const Eigen::MatrixBase<DerivedF_ref> & F_ref,
	Eigen::MatrixBase<DerivedD> & distances,
	int & hits,
	int & misses);

using namespace std;

// Find_distance will not always be able to find an intersection between the ray and the other surface, this function interpolates the missing distances
//
// Inputs:
//   adj	adjaceny matrix 
//	 misses	 number of points to be interpolated
// Output:
//   distances	#V dim vector of the distances from each point in a to B
//
void interpolate_surfaces(
	const vector<vector<int>> adj,
	Eigen::VectorXd & distances,
	int & misses);

double displace(
	double x,
	double xmin,
	double xmax,
	double ymin,
	double ymax);

// moves each point in V based on the value in distances, the set of distances from find_distance
//
// Inputs:
//   V  #V by dim matrix of vertex coordinates
//	 N_smooth  #V by 3 dim matrix of normal vectors
//	 xmin  all distances smaller than xmin will be rounded up to xmin
//	 xmax  all distances larger than xmax will be rounded up to xmax
//	 ymin  the amount by which a vertex corresponding xmin will be moved
//	 ymax  the amount by which a vertex corresponding xmax will be moved
//	 
// Output:
//   V  #V by dim matrix of vertex coordinates in A
//

void displace_vertices(
	Eigen::MatrixXd & V,
	const Eigen::MatrixXd &N_smooth,
	const Eigen::VectorXd &distances,
	double xmin,
	double xmax,
	double ymin,
	double ymax);

#ifndef IGL_STATIC_LIBRARY
#  include "interpolate.cpp"
#endif

#endif
