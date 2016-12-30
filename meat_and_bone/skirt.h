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


// Smooth vertex attributes using uniform Laplacian
// Inputs:
//   Ain  #V by #A eigen Matrix of mesh vertex attributes (each vertex has #A attributes)
//   F    #F by 3 eigne Matrix of face (triangle) indices
// Output:
//   Aout #V by #A eigen Matrix of mesh vertex attributes
using namespace std;


bool make_skirt(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & N_smooth,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop);

bool make_skirt(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & N_smooth,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	Eigen::VectorXi & borderLoop_out,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop
);


bool explore_border_loop(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXd &Normals,
	Eigen::VectorXi & borderStatus,
	const  vector< vector<int>> adj,
	Eigen::VectorXi & borderLoop,
	Eigen::MatrixXd & bitangent,
	const int v_border_1st);

int locate_borders(
	const vector<vector<int>> adj,
	const Eigen::VectorXi & adjTriangles,
	const int n,
	Eigen::VectorXi & borderStatus,
	int & v_border_1st);

void find_adj_faces(
	const Eigen::MatrixXi & F,
	Eigen::VectorXi & adjTriangles);


void make_bitangents(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & borderLoop,
	const Eigen::MatrixXd & normals,
	Eigen::MatrixXd & bitangent);

void add_skirt_layer(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & bitangent,
	const Eigen::VectorXi & borderLoop,
	const Eigen::VectorXi & borderLoop_orig,
	const Eigen::MatrixXd & N_smooth,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	const double delta_norm,
	const double delta_bitan);

bool norms_are_oriented(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & bitangent,
	const Eigen::VectorXi & borderLoop);


#ifndef IGL_STATIC_LIBRARY
#  include "skirt.cpp"
#endif

#endif
