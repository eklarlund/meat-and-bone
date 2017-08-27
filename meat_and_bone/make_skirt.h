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
	const Eigen::MatrixXd & borderNormals,
	const Eigen::MatrixXd & borderBitangents,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	Eigen::VectorXi & borderLoop,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop
);

bool make_skirt(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & N,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	Eigen::VectorXi & borderLoop,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop
);

void make_bitangents(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & borderLoop,
	const Eigen::MatrixXd & normals,
	Eigen::MatrixXd & bitangent);

void add_skirt_layer(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & borderLoop,
	const Eigen::MatrixXd & borderNormals,
	const Eigen::MatrixXd & borderBitangents,
	const double delta_norm,
	const double delta_bitan,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus
);

bool norms_are_oriented(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & bitangent,
	const Eigen::VectorXi & borderLoop);

#endif
