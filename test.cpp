#include <cfloat>
#include <iostream>
#include <math.h>
#include <vector>

#include <igl/per_vertex_normals.h>
#include <igl/writeOBJ.h>

#include <Eigen/Dense>

#include <meat_and_bone/interpolate.h>
#include <meat_and_bone/make_solid.h>
#include <meat_and_bone/make_skirt.h>
#include <meat_and_bone/laplace_smooth.h>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

const double tolerance = .000001;
bool equals(double a, double b)
{
	return (a - b < tolerance) && (b - a < tolerance);
}

TEST_CASE("All", "[]")
{
	Eigen::MatrixXd V_work(4, 3);
	V_work << 0, 0, 0,
		0, 1, 0,
		1, 0, 0,
		2, 0, 0;
	Eigen::MatrixXi F_work(2, 3);
	F_work << 0, 2, 1,
		2, 3, 1;

	Eigen::MatrixXd V_ref(4, 3);
	V_ref << -1, -1, .6,
		0, 2, .6,
		2, 0, .6,
		4, 0, .6;
	Eigen::MatrixXi F_ref(2, 3);
	F_ref << 0, 2, 1,
		2, 3, 1;

	Eigen::MatrixXd N(4, 3);
	N << 0, 0, 1,
		0, 0, 1,
		0, 0, 1,
		0, 0, 1;

	Eigen::VectorXd distances(4);
	int hits = 0, misses = 0;
	find_distance(V_work, N, V_ref, F_ref, distances, hits, misses);

	for (int i = 0; i < distances.rows(); i++)
	{
		REQUIRE(equals(distances(i), .6)); //.6 is distance between surfaces
	}

	displace_vertices(V_work, N, distances, 0, 1, 0, 1);

	for (int i = 0; i < distances.rows(); i++)
	{
		REQUIRE(equals(V_work(i, 2), .648)); //.648 is the value of smoothstep for x = .6
	}

	Eigen::MatrixXd N_smooth, V_smooth;
	laplace_smooth(V_ref, F_ref, V_smooth, N_smooth, 1);

	for (int i = 0; i < N.rows(); i++)
	{
		REQUIRE(equals(N(i, 0), -N_smooth(i, 0)));
		REQUIRE(equals(N(i, 1), -N_smooth(i, 1)));
		REQUIRE(equals(N(i, 2), -N_smooth(i, 2)));
	}

	Eigen::MatrixXd V_plus = V_work;
	Eigen::MatrixXi F_plus;
	Eigen::VectorXi borderLoop(4);
	borderLoop << 0, 2, 3, 1;
	

	REQUIRE(make_skirt(V_work, F_work, N, V_plus, F_plus, borderLoop, 3, 1, -1, 0));
	std::cout <<V_plus(13, 2) << "      " <<  V_work(0, 2) << "     " << V_plus(0, 2) <<  "\n";
	REQUIRE(equals(V_plus(13, 2), V_work(0, 2) - 1));
	REQUIRE(equals(V_plus(14, 2), V_work(1, 2) - 1));
	REQUIRE(equals(V_plus(15, 2), V_work(2, 2) - 1));

	Eigen::MatrixXd V_out;
	Eigen::MatrixXi F_out;
	REQUIRE(make_solid(V_plus, F_plus, borderLoop, V_out, F_out, 1, .5));

	REQUIRE(equals(V_out(16, 2), V_plus(0, 2)) - .5);
}
