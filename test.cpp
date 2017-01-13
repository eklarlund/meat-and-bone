#include <cfloat>
#include <iostream>
#include <math.h> 
#include <vector>

#include <igl/per_vertex_normals.h>
#include <igl/writeOBJ.h>

#include <Eigen/Dense>

#include <meat_and_bone/interpolate.h>
#include <meat_and_bone/thicken.h>
#include <meat_and_bone/skirt.h>
#include <meat_and_bone/laplace_smooth.h>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>


TEST_CASE("Factorials are computed", "[factorial]") {
	Eigen::MatrixXd V_work(4,3);
	V_work << 0, 0, 0,
		0, 1, 0,
		1, 0, 0,
		2, 0, 0;
	Eigen::MatrixXi F_work(2,3);
	F_work << 0, 2, 1,
			2, 3, 1;

	Eigen::MatrixXd V_ref(4, 3);
	V_ref << -1, -1, 1,
		0, 2, 1,
		2, 0, 1,
		4, 0, 1;
	Eigen::MatrixXi F_ref(2, 3);
	F_ref << 0, 2, 1, 
			2, 3, 1;

	cout << 1 << endl;

	Eigen::MatrixXd N(3,3);
	N << 0, 0, 1,
		0, 0, 1,
		0, 0, 1;

	cout << 2 << endl;

	igl::writeOBJ("test_work.obj", V_work, F_work);
	igl::writeOBJ("test_ref.obj", V_ref, F_ref);

	Eigen::VectorXd distances(4);

	cout << 3 << endl;

	int hits = 0, misses = 0;
	find_distance(V_work, N, V_ref, F_ref, distances, hits, misses);
	REQUIRE(distances(1) - 1 < .000001);
}