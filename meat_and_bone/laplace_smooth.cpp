#include "laplace_smooth.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/orientable_patches.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_vertex_normals.cpp>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/orient_outward.h>
#include <igl/orient_outward.cpp>
#include <vector>

template <typename DerivedV, typename DerivedF>

void laplace_smooth(
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedV> & N_smooth,
	double smooth_factor)
 {

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);


	// Diagonal per-triangle "mass matrix"
	Eigen::VectorXd dblA;
	igl::doublearea(V, F, dblA);
	// Place areas along diagonal #dim times
	const auto & T = 1.0 * (dblA.replicate(3, 1)*0.5).asDiagonal();

	// Recompute just mass matrix on each step
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M); //Replace V1 with U if looping
																// Solve (M-delta*L) U = M*U
	const auto & S = (M - smooth_factor*L);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
	assert(solver.info() == Eigen::Success);
	Eigen::MatrixXd U;
	U = solver.solve(M * V).eval();
	// Compute centroid and subtract (also important for numerics)

	igl::doublearea(U, F, dblA);
	double area = 0.5 * dblA.sum();
	U.array() /= sqrt(area);


	//Flip flipped normals
	Eigen::MatrixXi C;
	igl::orientable_patches(F, C);
	Eigen::MatrixXi FF; //F matrix, but oriented so that all normals point away from center of mass, inconsistent orientation
	Eigen::VectorXi I;
	igl::orient_outward(U, F, C, FF, I);



	igl::per_vertex_normals(U, FF, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, N_smooth);//Use vertices of smoothed surface
}


#ifdef IGL_STATIC_LIBRARY
	template void laplace_smooth
	<Eigen::Matrix<double, -1, -1, 0, -1, -1>,
	Eigen::Matrix<int, -1, -1, 0, -1, -1> >(
		Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&,
		Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&,
		Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&,
		double);
#endif