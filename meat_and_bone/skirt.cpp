#include "skirt.h"
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

int triangle_number(int i) {
	return i*(i + 1) / 2;
}

// (TODO): needed?
bool make_skirt(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & N_smooth,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop
) {
	Eigen::VectorXi borderLoop;
	return make_skirt(V, F, N_smooth, V_plus, F_plus, borderLoop, skirts, displacement, offset, drop);
}


bool make_skirt(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & N_smooth,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	Eigen::VectorXi & borderLoop,
	const int skirts,
	const double displacement,
	const double offset,
	const double drop
)
{
	if (skirts < 1) {
		cerr << "\n MUST ADD POSITIVE NUMBER OF SKIRTS";
		return false;
	}
	std::vector<std::vector<int> > adj;
	igl::adjacency_list(F, adj);

	double denominator = triangle_number(skirts);

	const int n = V.rows();
	const int n_faces = F.rows();
	Eigen::VectorXi adjTriangles(n);//Number of triangles adjacent to Vertex in  V at same index
	Eigen::VectorXi borderStatus(n);

	find_adj_faces(F, adjTriangles);
	int v_border_1st = -1;
	int num_border_v = locate_borders(adj, adjTriangles, n, borderStatus, v_border_1st);

	if (num_border_v == -1)
		return false;

	borderLoop.conservativeResize(num_border_v);
	borderLoop.fill(-1);

	Eigen::MatrixXd bitangent(num_border_v, 3);
	if (!explore_border_loop(V,F,N_smooth,borderStatus, adj, borderLoop, bitangent, v_border_1st))
	{
		return false;
	}
	cout << "BITANGENT FINAL: " << bitangent.row(bitangent.rows() - 1) << endl;
	const Eigen::VectorXi borderLoop_orig = borderLoop;

	Eigen::MatrixXd V_prev = V;
	Eigen::MatrixXi F_prev = F;

	V_prev.conservativeResize(V.rows() + skirts * num_border_v, 3);
	F_prev.conservativeResize(F.rows() + 2 * skirts * num_border_v, 3);

	double delta_bitan = displacement / skirts;
	double delta_norm = offset / denominator + drop;
	cout << "Deltabitan = " << delta_bitan <<  endl;

	
	for (int i = 0; i < skirts; i++)
	{
		cout << "\ndelta bitan, norm: " << delta_bitan << ", " << delta_norm;
		int F_curr = n_faces + num_border_v*i*2;
		add_skirt_layer(n + i * num_border_v, F_curr, V_prev, F_prev, bitangent, borderLoop, borderLoop_orig,
			N_smooth, V_plus, F_plus, delta_norm, delta_bitan);
		for (int j = 0; j < num_border_v; j++) {
			borderLoop(j) = (i)*num_border_v + n + j;
		}
		V_prev = V_plus;
		F_prev = F_plus;

		delta_norm += (i + 1) * offset / denominator - ((i == 0) ? drop : 0);

		cout << "\n skirt layer " << i << " added\n";
	}
	return true;
}





bool explore_border_loop(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXd &Normals,
	Eigen::VectorXi & borderStatus,
	const  vector< vector<int>> adj, 
	Eigen::VectorXi & borderLoop,
	Eigen::MatrixXd & bitangent,
	const int v_border_1st) {

	borderLoop.fill(-1);
	borderLoop(0) = v_border_1st;
	borderStatus(borderLoop(0)) = 2;

	int i_border = 0;//index for borderLoop
	const int num_border_v = borderLoop.rows();
	bool foundNext = true;
	while (foundNext)
	{
		foundNext = false;
		int v_curr = borderLoop(i_border);
		int v_curr_degree = adj[v_curr].size();
		for (int i = 0; i < v_curr_degree; i++)
		{//Searches for next index
			int v_cand = adj[v_curr][i]; //Index in V1 of current vertex being investigated
			if (borderStatus(v_cand) == 1)
			{
				int shared_v = 0;
				for (int k = 0; k < adj[v_curr].size(); k++) {
					for (int j = 0; j < adj[v_cand].size(); j++) {
						if (adj[v_curr][k] == adj[v_cand][j]) {
							shared_v++;
						}
					}
				}
				if (shared_v == 1)
				{
					borderStatus(v_cand) = 2; //Marks down that borderStatus has already been added to borderLoop
					borderLoop(++i_border) = v_cand;
					foundNext = true;
					break;
				}
				else
					assert(shared_v == 2);
			}

		}	
	}
	make_bitangents(V, F, borderLoop, Normals, bitangent);
	if (!norms_are_oriented(V, F, bitangent, borderLoop)) {
		cout << endl << "Flipping bitans" << endl;
		borderLoop.reverseInPlace();
		bitangent = bitangent.colwise().reverse().eval();
		bitangent = -bitangent;
	}
	//After processing all neighbors
	if (i_border + 1 == num_border_v)
	{
		cout << "\nOrdered Edge complete";
		return true;
	}
	cout << "\nOrdered Edge Error";
	cerr << "borderLoop(i_border): " << borderLoop(i_border) << "i_border: " << i_border;
	cerr << "\nv_border=" << i_border << "borderStatus(i_border): " << borderStatus(borderLoop(i_border));
	for (int k = 0; k < adj[borderLoop(i_border)].size(); k++)
	{
		int errIndex = adj[borderLoop(i_border)][k]; //index of error in V1
		cerr << "\n adj: " << errIndex << " borderStatus: " << borderStatus(errIndex);
	}
	return false;

}



int locate_borders(
	const  vector<vector<int>> adj,
	const Eigen::VectorXi & adjTriangles, 
	const int n,
	Eigen::VectorXi & borderStatus, 
	int & v_border_1st) {//Returns number of border Vertices

	int num_border_v = 0;
	borderStatus.resize(n);
	borderStatus.fill(-1);
	for (int i = 0; i < n; i++)
	{
		if (adj[i].size() == adjTriangles(i))
		{
			borderStatus(i) = 0;
		}
		else {
			assert(adj[i].size() == adjTriangles(i) + 1);
			borderStatus(i) = 1;
			num_border_v++;
			v_border_1st = i;
		}
	}
	 cout << "\n n_border_v= " << num_border_v << endl;

	if (v_border_1st == -1)
	{
		 cerr << "\nFound closed surface\n";
		return -1;
	}
	return num_border_v;
}

void find_adj_faces(
	const Eigen::MatrixXi & F, 
	Eigen::VectorXi & adjTriangles) {


	adjTriangles.fill(0);
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			int index = F(i, j);
			++adjTriangles(index);
		}
	}
}


void make_bitangents(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & borderLoop, 
	const Eigen::MatrixXd & normals,
	Eigen::MatrixXd & bitangent) {


	Eigen::Vector3d curve_tangent; //Vector from previous edge to next edge
	Eigen::Vector3d norm_curr; //current normal
	int num_border_v = borderLoop.rows();


	for (int i = 0; i < num_border_v; i++)
	{//Makes Cross products
		norm_curr = normals.row(borderLoop(i));
		const int i_next = (i + 1) % num_border_v;
		const int i_prev = (i + num_border_v - 1) % num_border_v;
		curve_tangent = (V.row(borderLoop(i_prev)) - V.row(borderLoop(i_next)));
		bitangent.row(i) = norm_curr.cross(curve_tangent).normalized();
	}
}

void add_skirt_layer(
	const int n,
	int f_curr,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & bitangent,
	const Eigen::VectorXi & borderLoop,
	const Eigen::VectorXi & borderLoop_orig,
	const Eigen::MatrixXd & N_smooth,
	Eigen::MatrixXd & V_plus,
	Eigen::MatrixXi & F_plus,
	const double delta_norm,
	const double delta_bitan
	) {

	
	const int num_border_v = borderLoop.rows();
	//int f_curr = F.rows();
	//const int n = V.rows();

	V_plus = V;
	F_plus = F;
	cout << "F_curr: " << f_curr << "F.rows: " << F.rows() << "num_border_v: " << num_border_v;
	//V_plus.conservativeResize(V_plus.rows() + num_border_v, 3);
	//F_plus.conservativeResize(F_plus.rows() + 2 * num_border_v, 3);

	for (int i = 0; i < num_border_v; i++)
	{//Adds 1st layer of edge
		int i_new = (i + 1) % num_border_v;
		int v_curr = borderLoop(i);
		int v_norm_curr = borderLoop_orig(i);
		V_plus.row(n + i) = V.row(v_curr) + delta_bitan*bitangent.row(i) + delta_norm * N_smooth.row(v_norm_curr);
		Eigen::Vector3i new_face(v_curr, n + i, n + i_new);
		Eigen::Vector3i new_face2(v_curr, n + i_new, borderLoop(i_new));
		F_plus.row(f_curr++) = new_face;
		F_plus.row(f_curr++) = new_face2;
	}
}


bool norms_are_oriented(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & bitangent,
	const Eigen::VectorXi & borderLoop
) {
	int corrects = 0, incorrects = 0;
	for (int times = 0; times < bitangent.rows() / 2; times++) {
		int v_1st = borderLoop(times);
		int v_2nd = borderLoop(times + 1);
		int v_3rd = -1;

		int f_index;
		for (int i = 0; v_3rd == -1; i++)
		{
			for (int j = 0; j < 3; j++) {
				if (v_1st == F(i, j)) {
					int corner2 = (j + 1) % 3;
					int corner3 = (j + 2) % 3;
					if (v_2nd == F(i, corner2)) {
						v_3rd = F(i, corner3);
						f_index = i;
					}
					else if (v_2nd == F(i, corner3)) {
						v_3rd = F(i, corner2);
						f_index = i;
					}
				}
			}
		}

		Eigen::Vector3d normal = bitangent.row(0);
		Eigen::Vector3d v_going_out = V.row(v_1st) - V.row(v_3rd);

		if (normal.dot(v_going_out) > 0) {
			corrects++;
		}
		else {
			incorrects++;
		}
	}
	cout << "Corrects: " << corrects << endl;
	cout << "Incorrects: " << incorrects << endl;
	if (corrects >= incorrects) {
		cout << "Orientation: correct" << endl;
		return true;
	}
	else {
		cout << "Orientation: incorrect" << endl;
		return false;
	}
}