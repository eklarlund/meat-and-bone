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

void find_distance(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXd & N_smooth,
	const Eigen::MatrixXd & V_ref,
	const Eigen::MatrixXi & F_ref,
	Eigen::VectorXd & distances,
	int & hits,
	int & misses);

using namespace std;

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
