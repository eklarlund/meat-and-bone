#include "make_skirt.h"
#include "make_bitangents.h"
#include "make_bitangents.cpp"
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

/*void make_normals_bitangents(const Eigen::MatrixXd &V,
							 const Eigen::MatrixXi &F,
							 const Eigen::MatrixXd &N,
							 Eigen::VectorXi &borderLoop,
							 Eigen::MatrixXd &borderNormals,
							 Eigen::MatrixXd &borderBitangents)
{
	make_bitangents(V, F, borderLoop, N, borderBitangents);

	// Pack the rows indexed by borderLoop into borderNormals
	igl::slice(N, borderLoop, 1, borderNormals);

	if (!norms_are_oriented(V, F, borderBitangents, borderLoop))
	{
		cout << endl
			 << "Flippping bitans" << endl;
		borderLoop.reverseInPlace();
		borderBitangents = borderBitangents.colwise().reverse().eval();
		borderBitangents = -borderBitangents;
		borderNormals = borderNormals.colwise().reverse().eval();
	};
}
*/
bool make_skirt(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXd &borderNormals,
	const Eigen::MatrixXd &borderBitangents,
	Eigen::MatrixXd &V_plus,
	Eigen::MatrixXi &F_plus,
	Eigen::VectorXi &borderLoop,
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

bool make_skirt(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXd &N,
	Eigen::MatrixXd &V_plus,
	Eigen::MatrixXi &F_plus,
	Eigen::VectorXi &borderLoop,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop)
{
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


/*void make_bitangents(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::VectorXi &borderLoop,
	const Eigen::MatrixXd &normals,
	Eigen::MatrixXd &borderBiTangents)
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
}
*/
void add_skirt_layer(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::VectorXi &borderLoop,
	const Eigen::MatrixXd &borderNormals,
	const Eigen::MatrixXd &borderBitangents,
	const double delta_norm,
	const double delta_bitan,
	Eigen::MatrixXd &V_plus,
	Eigen::MatrixXi &F_plus)
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
/*
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
}*/

