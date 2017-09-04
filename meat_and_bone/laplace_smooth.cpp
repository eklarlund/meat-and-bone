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

template <typename DerivedV, typename DerivedF, typename DerivedVS, typename DerivedNS>
void laplace_smooth(
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::PlainObjectBase<DerivedF>& F,
	Eigen::PlainObjectBase<DerivedVS>& V_smooth,
	Eigen::PlainObjectBase<DerivedNS> & N_smooth,
	double smooth_factor)
 {

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// Recompute just mass matrix on each step
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M); //Replace V1 with U if looping
																// Solve (M-delta*L) U = M*U
	const auto & S = (M - smooth_factor*L);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
	assert(solver.info() == Eigen::Success);
	V_smooth = solver.solve(M * V).eval();

	//Flip flipped normals
	Eigen::MatrixXi C;
	igl::orientable_patches(F, C);
	assert(C.maxCoeff() == 0);
	Eigen::MatrixXi FF; //F matrix, but oriented so that all normals point away from center of mass, inconsistent orientation
	Eigen::VectorXi I;
	igl::orient_outward(V_smooth, F, C, FF, I);
	igl::per_vertex_normals(V_smooth, FF, N_smooth);
}


#ifdef IGL_STATIC_LIBRARY
	template void laplace_smooth
	<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd>(
		const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		const Eigen::PlainObjectBase<Eigen::MatrixXi> &,
		Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		Eigen::PlainObjectBase<Eigen::MatrixXd>&,
		double);
#endif