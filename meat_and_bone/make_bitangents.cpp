#include "make_bitangents.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <iostream>

using namespace std;



bool norms_are_oriented(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXd &borderBiTangents,
	const Eigen::VectorXi &borderLoop)
{
	int corrects = 0, incorrects = 0;
	for (int times = 0; times < borderBiTangents.rows() / 2; times++)
	{
		int v_1st = borderLoop(times);
		int v_2nd = borderLoop(times + 1);
		int v_3rd = -1;

		int f_index;
		for (int i = 0; v_3rd == -1; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (v_1st == F(i, j))
				{
					int corner2 = (j + 1) % 3;
					int corner3 = (j + 2) % 3;
					if (v_2nd == F(i, corner2))
					{
						v_3rd = F(i, corner3);
						f_index = i;
					}
					else if (v_2nd == F(i, corner3))
					{
						v_3rd = F(i, corner2);
						f_index = i;
					}
				}
			}
		}

		Eigen::Vector3d normal = borderBiTangents.row(times);
		Eigen::Vector3d v_going_out = V.row(v_1st) - V.row(v_3rd);

		if (normal.dot(v_going_out) > 0)
		{
			corrects++;
		}
		else
		{
			incorrects++;
		}
	}
	cout << "Normals pointing out: " << corrects << endl;
	cout << "Normals pointing in: " << incorrects << endl;
	if (corrects >= incorrects)
	{
		cout << "Orientation: correct" << endl;
		return true;
	}
	else
	{
		cout << "Orientation: incorrect" << endl;
		return false;
	}
}

template <typename DerivedV, typename DerivedF, typename DerivedN, typename DerivedBL, typename DerivedBN, typename DerivedBB>
void make_normals_bitangents(
	const Eigen::MatrixBase<DerivedV> &V,
	const Eigen::MatrixBase<DerivedF> &F,
	const Eigen::MatrixBase<DerivedN> &N,
	Eigen::MatrixBase<DerivedBL> &borderLoop,
	Eigen::MatrixBase<DerivedBN> &borderNormals,
	Eigen::MatrixBase<DerivedBB> &borderBitangents)
{
	make_bitangents(V, F, borderLoop, N, borderBitangents);

	// Pack the rows indexed by borderLoop into borderNormals	
	if (borderNormals.rows() != borderLoop.rows()) {
		borderNormals.resize(borderLoop.rows());
	}
	for (int i = 0; i < borderLoop.rows(); i++) {
		borderNormals.row(i) = N.row(borderLoop(i));
	}
}


template <typename DerivedV, typename DerivedF, typename DerivedBL, typename DerivedN, typename DerivedBB>
void make_bitangents(
	const Eigen::MatrixBase<DerivedV> &V,
	const Eigen::MatrixBase<DerivedF> &F,
	Eigen::MatrixBase<DerivedBL> &borderLoop,
	const Eigen::MatrixBase<DerivedN> &normals,
	Eigen::MatrixBase<DerivedBB> &borderBiTangents)
{

	Eigen::Vector3d curve_tangent; //Vector from previous edge to next edge
	Eigen::Vector3d norm_curr;	 //current normal
	int num_border_v = borderLoop.rows();

	for (int i = 0; i < num_border_v; i++)
	{ //Makes Cross products
		norm_curr = normals.row(borderLoop(i));
		const int i_next = (i + 1) % num_border_v;
		const int i_prev = (i + num_border_v - 1) % num_border_v;
		curve_tangent = (V.row(borderLoop(i_prev)) - V.row(borderLoop(i_next)));
		borderBiTangents.row(i) = norm_curr.cross(curve_tangent).normalized();
	}

	if (!norms_are_oriented(V, F, borderBiTangents, borderLoop))
	{
		cout << endl << "Flippping bitans" << endl;
		borderLoop.reverseInPlace();
		borderBiTangents = borderBiTangents.colwise().reverse().eval();
		borderBiTangents = -borderBiTangents;
	}
}



#ifdef IGL_STATIC_LIBRARY
template void make_normals_bitangents<
	Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd>(
		const Eigen::MatrixBase<Eigen::MatrixXd> &,
		const Eigen::MatrixBase<Eigen::MatrixXi> &,
		const Eigen::MatrixBase<Eigen::MatrixXd> &,
		Eigen::MatrixBase<Eigen::VectorXi> &,
		Eigen::MatrixBase<Eigen::MatrixXd> &,
		Eigen::MatrixBase<Eigen::MatrixXd> &);


template void make_bitangents<
	Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd>(
		const Eigen::MatrixBase<Eigen::MatrixXd> &,
		const Eigen::MatrixBase<Eigen::MatrixXi> &,
		Eigen::MatrixBase<Eigen::VectorXi> &,
		const Eigen::MatrixBase<Eigen::MatrixXd> &,
		Eigen::MatrixBase<Eigen::MatrixXd> &);

#endif