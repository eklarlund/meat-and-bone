#include "interpolate.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/orientable_patches.h>
#include <iostream>
#include <igl/per_vertex_normals.h>
#include <igl/per_vertex_normals.cpp>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <cfloat>
#include <vector>
#include <igl/ray_mesh_intersect.h>
#include <igl/ray_mesh_intersect.cpp>

using namespace std;
template <typename DerivedV, typename DerivedN, typename DerivedV_ref, typename DerivedF_ref, typename DerivedD>

void find_distance(
	const Eigen::MatrixBase<DerivedV> & V,
	const Eigen::MatrixBase<DerivedN> & N_smooth,
	const Eigen::MatrixBase<DerivedV_ref> & V_ref,
	const Eigen::MatrixBase<DerivedF_ref> & F_ref,
	Eigen::MatrixBase<DerivedD> & distances,
	int & hits,
	int & misses) {

	int num_dots = N_smooth.rows() / 20 > 2 ? N_smooth.rows() / 20 : 2;
	distances.fill(DBL_MAX);
	igl::Hit hit;

	for (int i = 0; i < N_smooth.rows(); i++)
	{
		if (igl::ray_mesh_intersect(V.row(i), N_smooth.row(i), V_ref, F_ref, hit))
		{
			hits++;
			distances(i) = hit.t;
			if (distances(i) < 0)
			{
//				std::cout << "distance<0\n";
			}
		}
		else
		{
			misses++;
		}
		if (i%num_dots == 0)
		{
			//std::cout << ".";
		}
	}
	//std::cout << "ray shooting done\n";
}

void interpolate_surfaces(
	const vector<vector<int>> adj,
	Eigen::VectorXd & distances,
	int & misses) {


	while (misses > 0)
	{
		for (int i = 0; i < distances.rows(); i++)
		{
			if (distances(i) == DBL_MAX) {
				int num_norms = 0;
				int sum_norms = 0;
				for (int j = 0; j < adj[i].size(); j++)
				{
					if (distances(adj[i][j]) != DBL_MAX) {
						sum_norms += distances(adj[i][j]);
						num_norms++;
					}
				}
				if (num_norms != 0)
				{
					distances(i) = sum_norms / num_norms;
					misses--;
				}
			}
		}
	}
	double avgD = distances.sum() / distances.rows();
}

double displace(
	double x, 
	double xmin,
	double xmax,
	double ymin,
	double ymax
)
{//calculates the magnitude of displacement, uses linear displacement with maximum and minimum
	double xRange = xmax - xmin, yRange = ymax - ymin;
	if (x > xmax)
	{
		return ymax;
	}
	else if (x < xmin)
	{
		return ymin;
	}
	else
		return yRange*(x - xmin) / xRange*(x - xmin) / xRange*(3 - 2 * (x - xmin) / xRange) + ymin;
}

void displace_vertices(
	Eigen::MatrixXd & V,
	const Eigen::MatrixXd & N_smooth,
	const Eigen::VectorXd & distances,
	double xmin,
	double xmax,
	double ymin,
	double ymax) {

	const int n = V.rows();
	for (int i = 0; i < n; i++)
	{
		V.row(i) += displace(distances(i), xmin, xmax, ymin, ymax)*N_smooth.row(i);
	}

}

#ifdef IGL_STATIC_LIBRARY
template void find_distance<
	Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd>(
		const Eigen::MatrixBase<Eigen::MatrixXd> &,
		const Eigen::MatrixBase<Eigen::MatrixXd> &,
		const Eigen::MatrixBase<Eigen::MatrixXd> &,
		const Eigen::MatrixBase<Eigen::MatrixXi> &,
		Eigen::MatrixBase<Eigen::VectorXd> &,
		int &,
		int &);
#endif