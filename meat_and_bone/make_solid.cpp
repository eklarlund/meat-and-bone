//#include <Eigen/core>
#include "laplace_smooth.h"
#include "make_solid.h"
#include <iostream>

using namespace std;

template <typename DerivedV, typename DerivedF, typename DerivedBL, typename DerivedV_out, typename DerivedF_out>

bool make_solid(
	const Eigen::PlainObjectBase<DerivedV> &V,
	const Eigen::PlainObjectBase<DerivedF> &F,
	const Eigen::PlainObjectBase<DerivedBL> &borderLoop,
	Eigen::PlainObjectBase<DerivedV_out> & V_out,
	Eigen::PlainObjectBase<DerivedF_out> & F_out,
	const double smoothf,
	const double offset) {

	const int n = V.rows();
	const int num_faces = F.rows();
	const int num_border_v = borderLoop.rows();
	Eigen::MatrixXd V_smooth, N;
	laplace_smooth(V, F, V_smooth, N, smoothf);

	Eigen::MatrixXd V_new = V;
	Eigen::MatrixXi F_new = F;

	V_out.resize(2 * V.rows(), V.cols());
	F_out.resize(2 * F.rows(), F.cols());

	for (int i = 0; i < n; i++) {
		V_new.row(i) = V.row(i) - offset * N.row(i);
	}

	for (int i = 0; i < num_faces; i++) {
		for (int j = 0; j < F.cols(); j++) {
			F_new(i, j) = F(i, 2 - j);
		}
	}

	V_out << V, V_new;
	if (offset > 0)
		F_out << F, (F_new.array() + V.rows());
	else
		F_out << F_new, (F.array() + V.rows());

	int f_curr = F_out.rows();
	F_out.conservativeResize(f_curr + 2 * num_border_v, F_out.cols());
	const int V_n = V.rows();

	for (int i = 0; i < num_border_v; i++)
	{
		int next_i = (i + 1) % num_border_v;
		Eigen::Vector3i new_face_1, new_face_2;
		if (offset > 0)
		{
			new_face_1 = Eigen::Vector3i(borderLoop(next_i), borderLoop(i), borderLoop(i) + n);
			new_face_2 = Eigen::Vector3i(borderLoop(next_i) + n, borderLoop(next_i), borderLoop(i) + n);
		}
		else
		{
			new_face_1 = Eigen::Vector3i(borderLoop(i) + n, borderLoop(i), borderLoop(next_i));
			new_face_2 = Eigen::Vector3i(borderLoop(i) + n, borderLoop(next_i), borderLoop(next_i) + n);
		}
		F_out.row(f_curr++) = new_face_1;
		F_out.row(f_curr++) = new_face_2;
	}
	return 1;
}

#ifdef IGL_STATIC_LIBRARY
template bool make_solid
<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXi>(
	const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
	const Eigen::PlainObjectBase<Eigen::MatrixXi> &,
	const Eigen::PlainObjectBase<Eigen::VectorXi> &,
	Eigen::PlainObjectBase<Eigen::MatrixXd> &,
	Eigen::PlainObjectBase<Eigen::MatrixXi> &,
	const double,
	const double);
#endif
