#include "make_skirt.h"
#include "make_bitangents.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/orientable_patches.h>
#include <iostream>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/orient_outward.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.cpp>
#include <igl/orient_outward.cpp>
#include <igl/adjacency_list.h>
#include <vector>

using namespace std;

int triangle_number(int i)
{
	return i * (i + 1) / 2;
}
template <typename DerivedV, typename DerivedF, typename DerivedBL, typename DerivedBN, typename DerivedBB, typename DerivedV_plus, typename DerivedF_plus>
void add_skirt_layer(
	const Eigen::PlainObjectBase<DerivedV> &V,
	const Eigen::PlainObjectBase<DerivedF> &F,
	const Eigen::PlainObjectBase<DerivedBL> &borderLoop,
	const Eigen::PlainObjectBase<DerivedBN> &borderNormals,
	const Eigen::PlainObjectBase<DerivedBB> &borderBitangents,
	const double delta_norm,
	const double delta_bitan,
	Eigen::PlainObjectBase<DerivedV_plus> &V_plus,
	Eigen::PlainObjectBase<DerivedF_plus> &F_plus)
{
	const int num_border_v = borderLoop.rows();
	int f_curr = F.rows();
	const int n = V.rows();

	V_plus = V;
	F_plus = F;
	V_plus.conservativeResize(V_plus.rows() + num_border_v, 3);
	F_plus.conservativeResize(F_plus.rows() + 2 * num_border_v, 3);

	for (int i = 0; i < num_border_v; i++)
	{ //Adds 1st layer of edge
		int i_new = (i + 1) % num_border_v;
		int v_curr = borderLoop(i);
		V_plus.row(n + i) = V.row(v_curr) + delta_bitan * borderBitangents.row(i) + delta_norm * borderNormals.row(i);
		Eigen::Vector3i new_face(v_curr, n + i, n + i_new);
		Eigen::Vector3i new_face2(v_curr, n + i_new, borderLoop(i_new));
		F_plus.row(f_curr++) = new_face;
		F_plus.row(f_curr++) = new_face2;
	}
}

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
		const double drop)
{
	if (skirts < 1)
	{
		cerr << "\n MUST ADD POSITIVE NUMBER OF SKIRTS";
		return false;
	}

	Eigen::MatrixXd V_prev = V;
	Eigen::MatrixXi F_prev = F;
	int n = V.rows();
	int num_border_v = borderLoop.size();
	double denominator = triangle_number(skirts);
	double delta_bitan = displacement / skirts;
	double delta_norm = offset / denominator + drop;
	double total_delta_norm = 0;
	cout << "Deltabitan = " << delta_bitan << endl;

	for (int i = 0; i < skirts; i++)
	{
		cout << "\ndelta bitan, norm: " << delta_bitan << ", " << delta_norm;

		add_skirt_layer(V_prev, F_prev, borderLoop, borderNormals, borderBitangents,
			delta_norm, delta_bitan, V_plus, F_plus);
		for (int j = 0; j < num_border_v; j++)
		{
			borderLoop(j) = (i)*num_border_v + n + j;
		}
		// TODO copy operation is expensive, move before loop
		V_prev.conservativeResize(V_plus.rows(), V_plus.cols());
		V_prev.conservativeResize(F_plus.rows(), F_plus.cols());
		V_prev = V_plus;
		F_prev = F_plus;
		total_delta_norm += delta_norm;
		delta_norm += offset / denominator - ((i == 0) ? drop : 0);

		cout << "\n skirt layer " << i << " added\n";
		cout << "total delta_norm: " << total_delta_norm << endl;
	}

	return true;
}

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
) {
	const int n = V.rows();
	const int n_border = borderLoop.rows();
	Eigen::MatrixXd borderNormals(n_border, 3);
	Eigen::MatrixXd borderBitangents(n_border, 3);
	make_normals_bitangents(
		V,
		F,
		N,
		borderLoop,
		borderNormals,
		borderBitangents);

	return make_skirt(V, F, borderNormals,
		borderBitangents,
		V_plus,
		F_plus,
		borderLoop,
		skirts,
		displacement,
		offset,
		drop);
}

#ifdef IGL_STATIC_LIBRARY
template bool make_skirt<
	Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
	(
		const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		const Eigen::PlainObjectBase<Eigen::MatrixXi> &,
		const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		Eigen::PlainObjectBase<Eigen::MatrixXi> &,
		Eigen::PlainObjectBase<Eigen::VectorXi> &,
		const int,
		const double,
		const double,
		const double);

template bool make_skirt<
	Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>(
		const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		const Eigen::PlainObjectBase<Eigen::MatrixXi> &,
		const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		Eigen::PlainObjectBase<Eigen::MatrixXd> &,
		Eigen::PlainObjectBase<Eigen::MatrixXi> &,
		Eigen::PlainObjectBase<Eigen::VectorXi> &,
		const int,
		const double,
		const double,
		const double);

template void add_skirt_layer<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi>
(
	const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
	const Eigen::PlainObjectBase<Eigen::MatrixXi> &,
	const Eigen::PlainObjectBase<Eigen::VectorXi> &,
	const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
	const Eigen::PlainObjectBase<Eigen::MatrixXd> &,
	const double,
	const double,
	Eigen::PlainObjectBase<Eigen::MatrixXd> &,
	Eigen::PlainObjectBase<Eigen::MatrixXi> &);
#endif