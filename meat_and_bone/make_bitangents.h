
#ifndef LIBIGL_STUFF_MAKE_BITANGENTS
#define LIBIGL_STUFF_MAKE_BITANGENTS
#include <Eigen/Core>

//Creates bitangents (normals pointing outwards along the edge of the surface) for a surface
//
// Inputs:
//   V    #V by 3 eigen Matrix of mesh vertices
//   F    #F by 3 eigne Matrix of face (triangle) indices
//	 N	  #V by 3 eigen Matrix of normal vectors for each element in V
//   borderLoop  Eigen vector containing the indices of the border elements in V, ordered in a Loop, 
//	             will be flipped to match the so that the normal cross the curve tangent points outwards if nessesary
//	 
// Output:
//	 borderNormals  Eigen matrix containing normals at each point in borderLoop
//   borderBitangents  #B by 3 matrix of border bitangents corresponding to borderLoop

template <typename DerivedV, typename DerivedF>

void make_normals_bitangents(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXd &N,
	Eigen::VectorXi &borderLoop,
	Eigen::MatrixXd &borderNormals,
	Eigen::MatrixXd &borderBitangents);

void make_bitangents(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & borderLoop,
	const Eigen::MatrixXd & normals,
	Eigen::MatrixXd & bitangent);

void mabtest(

);

#ifndef IGL_STATIC_LIBRARY
#  include "make_bitangents.cpp"
#endif

#endif