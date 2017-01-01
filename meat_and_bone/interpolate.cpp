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

void find_distance(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXd & N_smooth,
	const Eigen::MatrixXd & V_ref,
	const Eigen::MatrixXi & F_ref,
	Eigen::VectorXd & distances, 
	int & hits,
	int & misses) {

	int num_dots = N_smooth.rows() / 20;
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
				std::cout << "distance<0\n";
			}
		}
		else
		{
			misses++;
		}
		if (i%num_dots == 0)
		{
			std::cout << ".";
		}
	}
	std::cout << "ray shooting done\n";
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
	std::cout << "\nAverage: " << avgD;
	std::cout << "\nmax: " << distances.maxCoeff();
	std::cout << "\nMin: " << distances.minCoeff();
}

double displace(
	double x, 
	double xmin,
	double xmax,
	double ymin,
	double ymax
)
{//calculates the magnitude of displacement
	double xs = xmax - xmin, ys = ymax - ymin;
	if (x > xmax)
	{
		return ymax;
	}
	else if (x < xmin)
	{
		return ymin;
	}
	else
		return ys*(x - xmin) / xs*(x - xmin) / xs*(3 - 2 * (x - xmin) / xs) + ymin;
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


