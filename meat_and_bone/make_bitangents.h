
#ifndef LIBIGL_STUFF_MAKE_BITANGENTS
#define LIBIGL_STUFF_MAKE_BITANGENTS
#include <Eigen/Core>
#include <vector>

//Finds bitangents (normals pointing outwards along the edge of the surface) and normal vectors on the edge of a surface
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

template <typename DerivedV, typename DerivedF, typename DerivedN, typename DerivedBL, typename DerivedBN, typename DerivedBB>
void make_normals_bitangents(
	const Eigen::MatrixBase<DerivedV> &V,
	const Eigen::MatrixBase<DerivedF> &F,
	const Eigen::MatrixBase<DerivedN> &N,
	Eigen::MatrixBase<DerivedBL> &borderLoop,
	Eigen::MatrixBase<DerivedBN> &borderNormals,
	Eigen::MatrixBase<DerivedBB> &borderBitangents);



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
template <typename DerivedV, typename DerivedF, typename DerivedBL, typename DerivedN, typename DerivedBB>
void make_bitangents(
	const Eigen::MatrixBase<DerivedV> &V,
	const Eigen::MatrixBase<DerivedF> &F,
	Eigen::MatrixBase<DerivedBL> &borderLoop,
	const Eigen::MatrixBase<DerivedN> &normals,
	Eigen::MatrixBase<DerivedBB> &borderBiTangents);


#ifndef IGL_STATIC_LIBRARY
#include "make_bitangents.cpp"
#endif

#endif