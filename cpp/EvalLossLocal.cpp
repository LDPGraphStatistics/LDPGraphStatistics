#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "MemoryOperation.h"
#include "include/stats.hpp"

using namespace std;

string EdgeFile;
int NodeNum;
double Eps;
string Eps_s;
double EpsNsMaxDeg;
int NSType;
double PosBias;
string PosBias_s;
int ItrNum;
int Alg;
double Balloc[3];
char *Balloc_s[3];

// Initialization of statslib
stats::rand_engine_t engine(1776);

FILE *FileOpen(string filename, const char *mode) {
	FILE *fp;

	if ((fp = fopen(filename.c_str(), mode)) == NULL) {
		cout << "cannot open " << filename << endl;
		exit(-1);
	}
	return fp;
}

bool checkFileExistence(const std::string& str) {
    std::ifstream ifs(str);
    return ifs.is_open();
}

// Randomly generate 0, 1, 2, ..., size-1, and store the first num values into rndperm
void MakeRndPerm(int *rndperm, int size, int num) {
	int rnd;
	int *ordperm;
	int i, j;

	// 0, 1, 2, ..., size-1 --> ordperm
	ordperm = (int *)malloc(size * sizeof(int));
	for (i = 0; i < size; i++) {
		ordperm[i] = i;
	}

	for (i = 0; i < num; i++) {
		rnd = genrand_int32() % (size - i);
		rndperm[i] = ordperm[rnd];
		for (j = rnd + 1; j < size - i; j++) {
			ordperm[j - 1] = ordperm[j];
		}
	}

	free(ordperm);
}

// Read edges from the edge file
void ReadEdges(map<int, int> *a_mat, int *node_order){
	int node1, node2;
	int i;
	char s[1025];
	char *tok;
	FILE *fp;

	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<3;i++) fgets(s, 1024, fp);
	while(fgets(s, 1024, fp) != NULL){
		// 1st node --> node1
		tok = strtok(s, ",");
		node1 = atoi(tok);
		// 2nd node --> node2
		tok = strtok(NULL, ",");
		node2 = atoi(tok);
		if(node1 == node2) continue;
		// If both nodes exist, add the edge
		if(node_order[node1] < NodeNum && node_order[node2] < NodeNum){
			a_mat[node_order[node1]][node_order[node2]] = 1;
			a_mat[node_order[node2]][node_order[node1]] = 1;
		}
	}
	fclose(fp);
}

// Calculate the noisy max degree
double CalcNSMaxDeg(int *deg, int &max_deg, double eps, string outfile){
	double *deg_lap;
	double max_deg_ns;
	int i;
	FILE *fp;

	// Initialization
	malloc1D(&deg_lap, NodeNum);

	// Calculate the max of Deg and Deg+Lap(eps) --> max_deg, max_deg_ns
	max_deg = 0;
	max_deg_ns = 0.0;
	for(i=0;i<NodeNum;i++){
		deg_lap[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/eps, engine);
		if(max_deg < deg[i]) max_deg = deg[i];
		if(max_deg_ns < deg_lap[i]) max_deg_ns = deg_lap[i];
	}
	// Add positive bias --> max_deg_ns
	max_deg_ns += PosBias;

	// If max_deg_ns exceeds NodeNum - 1, then use NodeNum - 1
	if(max_deg_ns > (double)NodeNum - 1.0) max_deg_ns = (double)NodeNum - 1.0;

	// free
	free1D(deg_lap);

	return max_deg_ns;
}

// Calculate #triangles in the centralized model
void CalcCentTri(long long tri_num, map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, double &sen_tri){
	int glb_max_deg;
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<tuple<int, int>, int> a2_mat;
	map<tuple<int, int>, int>::iterator a2_itr;
	double EpsLap;
	int max_deg;
	double max_deg_ns;
	int i, j, k;
	FILE *fp;

	// Initialization
	tri_num_ns = tri_num;

    // Lap (#nodes)
    if(NSType == 0){
        // Global sensitivity using #nodes --> sen_tri
        glb_max_deg = NodeNum - 1;
        sen_tri = (double)glb_max_deg - 1.0;

        // Add Lap(sen_tri/Eps) --> tri_num_ns
		tri_num_ns += stats::rlaplace(0.0, sen_tri/Eps, engine);
	}
	// Lap (max degree)
	else if(NSType == 1){
        // Global sensitivity using max degree --> sen_tri
		// max(deg) --> max_deg
		max_deg = 0;
		for(i=0;i<NodeNum;i++){
			if(max_deg < deg[i]) max_deg = deg[i];
		}
		sen_tri = (double)max_deg;

		// Add Lap(sen_tri/Eps) --> tri_num_ns
		tri_num_ns += stats::rlaplace(0.0, sen_tri/Eps, engine);
	}
}

// Calculate #2-stars and #3-stars in the centralized model
void CalcCentSt(long long st2_num, long long st3_num, map<int, int> *a_mat, int *deg, string outfile, double &st2_num_ns, double &st3_num_ns, double &sen_st2, double &sen_st3){
	int glb_max_deg;
	int max_deg;
	int sen_st2_ij, sen_st3_ij;
	double EpsLap;
	double max_deg_ns;
	int i, j;
	FILE *fp;

    // Initialization
    st2_num_ns = st2_num;
    st3_num_ns = st3_num;

    // Lap (#nodes)
    if(NSType == 0){
		// Global sensitivity using #nodes --> sen_st2, sen_st3
        glb_max_deg = NodeNum - 1;
        sen_st2 = 2.0 * (double)glb_max_deg;
        sen_st3 = (double)glb_max_deg * ((double)glb_max_deg - 1.0);

		// Add Lap(sen_tri/Eps) --> tri_num_ns
		st2_num_ns += stats::rlaplace(0.0, sen_st2/Eps, engine);
		st3_num_ns += stats::rlaplace(0.0, sen_st3/Eps, engine);
	}
	// Lap (max degree)
	else if(NSType == 1){
        // Global sensitivity using max degree --> sen_st2, sen_st3
        sen_st2 = sen_st3 = 0;
		// max(deg) --> max_deg
		max_deg = 0;
		for(i=0;i<NodeNum;i++){
			if(max_deg < deg[i]) max_deg = deg[i];
		}
		sen_st2 = 2.0 * (double)max_deg;
		sen_st3 = (double)max_deg * ((double)max_deg - 1.0);

		// Add Lap(sen_tri/Eps) --> st2_num_ns, st3_num_ns
		st2_num_ns += stats::rlaplace(0.0, sen_st2/Eps, engine);
		st3_num_ns += stats::rlaplace(0.0, sen_st3/Eps, engine);
	}
}

// Calculate #triangles in the non-interactive local model (RR)
void CalcNLocTri(map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, int emp){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg_ns;									// noisy degree
	long long tot_edge_num_ns;
	long long tri_num, st2_num, ed2_num, ed1_num, non_num;
	double q;
	double alp, alp_1_3, q_inv_11, q_inv_21, q_inv_31, q_inv_41;
	int i, j, k;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	malloc1D(&deg_ns, NodeNum);

	// Flip probability --> q
    q = 1.0 / (exp(Eps) + 1.0);

	// Flip 0/1 in a_mat with probability q --> a_mat_ns
	for(i=0;i<NodeNum;i++){
		for(j=i+1;j<NodeNum;j++){
			// 0 --> 1 (flip)
			if(genrand_real2() < q && a_mat[i].count(j) == 0){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(genrand_real2() >= q && a_mat[i].count(j) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
		}
	}

	// Degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) deg_ns[i] += 1;
	}

	// Total number of edges --> tot_edge_num_ns
	tot_edge_num_ns = 0;
	for(i=0;i<NodeNum;i++) tot_edge_num_ns += (long long)deg_ns[i];
	tot_edge_num_ns /= 2;

	// #triangles --> tri_num
	tri_num = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) {
			j = aitr->first;
			if (i >= j) continue;
			for (aitr2 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++) {
				k = aitr2->first;
				if (j >= k) continue;
				if(a_mat_ns[j].count(k) > 0) tri_num++;
			}
		}
	}

	// With empirical estimation
	if(emp == 1){
		// #2-stars --> st2_num
		st2_num = 0;
		for(i=0;i<NodeNum;i++){
			st2_num += ((long long)deg_ns[i] * ((long long)deg_ns[i]-1)) / 2;
		}

		// #2-edges --> ed2_num
		ed2_num = st2_num - 3*tri_num;
		// #1-edge --> ed1_num
		ed1_num = (long long)tot_edge_num_ns*(NodeNum-2) - 2*ed2_num - 3*tri_num;
		// #none --> non_num
		non_num = (long long)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num - ed2_num - ed1_num;

		alp = exp(Eps);
		alp_1_3 = (alp-1.0)*(alp-1.0)*(alp-1.0);
		q_inv_11 = (alp*alp*alp) / alp_1_3;
		q_inv_21 = - alp*alp / alp_1_3;
		q_inv_31 = alp / alp_1_3;
		q_inv_41 = - 1.0 / alp_1_3;

		tri_num_ns = (double)tri_num * q_inv_11 + (double)ed2_num * q_inv_21 + (double)ed1_num * q_inv_31 + (double)non_num * q_inv_41;
	}
	// Without empirical estimation
	else{
		tri_num_ns = (double)tri_num;
	}

	delete[] a_mat_ns;
	free1D(deg_ns);
}

// Calculate #2-stars and #3-stars in the non-interactive local model
void CalcNLocSt(long long st2_num, long long st3_num, int *deg, string outfile, double &st2_num_ns, double &st3_num_ns, double &sen_st2, double &sen_st3){
	int glb_max_deg;
	int max_deg;
	int sen_st2_ij, sen_st3_ij;
	double EpsLap;
	double max_deg_ns;
	int max_deg_ns_floor;
	int del_num;
	int *rndperm;
	int i, j, x;
	FILE *fp;

    // Initialization
    st2_num_ns = st2_num;
    st3_num_ns = st3_num;

    // Lap (#nodes)
    if(NSType == 0){
		// Global sensitivity using #nodes --> sen_st2, sen_st3
        glb_max_deg = NodeNum - 1;
        sen_st2 = (double)glb_max_deg;
        sen_st3 = (double)glb_max_deg * ((double)glb_max_deg - 1.0) / 2.0;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_st2/Eps) --> st2_num_ns
			st2_num_ns += stats::rlaplace(0.0, sen_st2/Eps, engine);
			// Add Lap(sen_st3/Eps) --> st3_num_ns
			st3_num_ns += stats::rlaplace(0.0, sen_st3/Eps, engine);
		}
	}
	// Lap (max degree)
	else if(NSType == 1){
        // Global sensitivity using max degree --> sen_st2, sen_st3
        sen_st2 = sen_st3 = 0;
		// max(deg) --> max_deg
		max_deg = 0;
		for(i=0;i<NodeNum;i++){
			if(max_deg < deg[i]) max_deg = deg[i];
		}
		sen_st2 = (double)max_deg;
		sen_st3 = (double)max_deg * ((double)max_deg - 1.0) / 2.0;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_st2/Eps) --> st2_num_ns
			st2_num_ns += stats::rlaplace(0.0, sen_st2/Eps, engine);
			// Add Lap(sen_st3/Eps) --> st3_num_ns
			st3_num_ns += stats::rlaplace(0.0, sen_st3/Eps, engine);
		}
	}
	// Lap (noisy max degree)
	else{
		// Eps = EpsNsMaxDeg + EpsLap
		EpsLap = Eps - EpsNsMaxDeg;
		// Calculate the noisy max degree --> max_deg_ns
		max_deg_ns = CalcNSMaxDeg(deg, max_deg, EpsNsMaxDeg, outfile);

		// If max(deg) exceeds max_deg_ns, then perform graph projection
		if((double)max_deg > max_deg_ns){
			max_deg_ns_floor = (int)floor(max_deg_ns);
			del_num = 0;
			// Calculate #2-stars, #3-stars again --> st2_num_ns, st3_num_ns
			st2_num_ns = st3_num_ns = 0;
			for(i=0;i<NodeNum;i++){
				if(deg[i] <= max_deg_ns_floor){
					st2_num_ns += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
					st3_num_ns += ((long long)deg[i] * ((long long)deg[i]-1) * ((long long)deg[i]-2)) / 6;
				}
				else{
					st2_num_ns += ((long long)max_deg_ns_floor * ((long long)max_deg_ns_floor-1)) / 2;
					st3_num_ns += ((long long)max_deg_ns_floor * ((long long)max_deg_ns_floor-1) * ((long long)max_deg_ns_floor-2)) / 6;
					del_num += (deg[i] - max_deg_ns_floor);
				}
			}
		}

        // Global sensitivity using noisy max degree --> sen_st2, sen_st3
        sen_st2 = max_deg_ns;
        sen_st3 = max_deg_ns * (max_deg_ns - 1.0) / 2.0;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_st2/Eps) --> st2_num_ns
			st2_num_ns += stats::rlaplace(0.0, sen_st2/EpsLap, engine);
			// Add Lap(sen_st3/Eps) --> st3_num_ns
			st3_num_ns += stats::rlaplace(0.0, sen_st3/EpsLap, engine);
		}
	}
}

// Calculate #triangles and #2-stars in the interactive local model
void CalcILocTri(map<int, int> *a_mat, int *deg, string outfile, double &tri_num_ns, double &sen_tri){
	double EpsLap, Eps1st, Eps2ndTr, Eps2ndSt, Eps2ndTrSt;
	double Balloc_sum = 0.0;
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int> *a_mat_del;			// deleted adjacency matrix after projection
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<tuple<int, int>, int> a2_mat;
	map<tuple<int, int>, int>::iterator a2_itr;
	double p, q;
	double *tri_num_u, *st2_num_u, *trist2_num_u;
	int glb_max_deg;
	int max_deg;
	int sen_st2_ij;
	double max_deg_ns;
	int max_deg_ns_floor;
	int del_num;
	int *rndperm;
	int i, j, k, x;
	FILE *fp;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	a_mat_del = new map<int, int>[NodeNum];
	malloc1D(&tri_num_u, NodeNum);
	malloc1D(&st2_num_u, NodeNum);
	malloc1D(&trist2_num_u, NodeNum);

	// Privacy budget allocation
	for(i=0;i<2;i++) Balloc_sum += Balloc[i];
	if(NSType == 0 || NSType == 1){
		Eps1st = Eps * Balloc[0] / Balloc_sum;
		Eps2ndTrSt = Eps * Balloc[1] / Balloc_sum;
	}
	else{
		EpsLap = Eps - EpsNsMaxDeg;
		Eps1st = EpsLap * Balloc[0] / Balloc_sum;
		Eps2ndTrSt = EpsLap * Balloc[1]  / Balloc_sum;
	}

	// Flip probability --> q
    q = 1.0 / (exp(Eps1st) + 1.0);
    p = 1.0 - q;

	// RR --> a_mat_ns
    // Count #noisy triangles and #noisy 2-stars for each user --> tri_num_u, st2_num_u
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
			j = aitr->first;
			for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
				k = aitr2->first;
				// it must be j < k < i.
				if (j >= k || k >= i) continue;

				st2_num_u[i] += 1.0;

				// If a_mat_ns[j][k] does not exist
				if(a_mat_ns[j].count(k) == 0){
					// Flip 0/1 in a_mat[j][k] with probability q --> a_mat_ns[j][k]
					// 0 --> 1 (flip)
					if(genrand_real2() < q && a_mat[j].count(k) == 0){
						a_mat_ns[j][k] = 1;
					}
					// 1 --> 1 (not flip)
					else if(genrand_real2() >= q && a_mat[j].count(k) == 1){
						a_mat_ns[j][k] = 1;
					}
					// 1 --> 0 (flip) or 0 --> 0 (not flip)
					else{
						a_mat_ns[j][k] = 0;
					}
				}
				if(a_mat_ns[j][k] == 1) tri_num_u[i] += 1.0;
			}
		}
	}

	// #triangles - (1 - p) * #2-stars --> trist2_num_u
	for(i=0;i<NodeNum;i++) trist2_num_u[i] = tri_num_u[i] - (1 - p) * st2_num_u[i];

    // Lap (#nodes)
	if(NSType == 0){
		// Global sensitivity using #nodes --> sen_tri
        glb_max_deg = NodeNum - 1;
		sen_tri = (double)glb_max_deg;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
			trist2_num_u[i] += stats::rlaplace(0.0, sen_tri/Eps2ndTrSt, engine);
		}
	}
	// Lap (max degree)
	else if(NSType == 1){
		// Global sensitivity using max degree --> sen_tri
		sen_tri = 0.0;
		// max(deg) --> max_deg
		max_deg = 0;
		for(i=0;i<NodeNum;i++){
			if(max_deg < deg[i]) max_deg = deg[i];
		}
		sen_tri = (double)max_deg;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
			trist2_num_u[i] += stats::rlaplace(0.0, sen_tri/Eps2ndTrSt, engine);
		}
	}
	// Lap (noisy max degree)
	else{
		// Calculate the noisy max degree --> max_deg_ns
		max_deg_ns = CalcNSMaxDeg(deg, max_deg, EpsNsMaxDeg, outfile);

		// If max(deg) exceeds max_deg_ns, then perform graph projection
		if((double)max_deg > max_deg_ns){
			max_deg_ns_floor = (int)floor(max_deg_ns);
			del_num = 0;
			// For each user
			for(i=0;i<NodeNum;i++){
				//If deg[i] exceeds max_deg_ns, then perform graph projection
				if(deg[i] > max_deg_ns){
					// Randomly generate 0, 1, ..., deg[i]-1 --> rndperm
					malloc1D(&rndperm, deg[i]);
					MakeRndPerm(rndperm, deg[i], deg[i]);

					// Randomly delete (deg[i] - max_deg_ns) edges from a_mat[i]
					x = 0;
					for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
						if(rndperm[x] >= max_deg_ns_floor){
							j = aitr->first;
							// Deleted edge --> a_mat_del[i][j]
							a_mat_del[i][j] = 1;
							del_num++;
						}
						x++;
					}
					free1D(rndperm);

					// Count #noisy triangles and #noisy 2-stars again --> tri_num_u, st2_num_u
					tri_num_u[i] = 0.0;
					st2_num_u[i] = 0.0;
					for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
						j = aitr->first;
						// Continue if the edge is deleted
						if(a_mat_del[i].count(j) == 1) continue;
						for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
							k = aitr2->first;
							// Continue if the edge is deleted
							if(a_mat_del[i].count(k) == 1) continue;
							if (j >= k || k >= i) continue;
							st2_num_u[i] += 1.0;
							if(a_mat_ns[j][k] == 1) tri_num_u[i] += 1.0;
						}
					}
				}
			}
			// #triangles - (1 - p) * #2-stars --> trist2_num_u
			for(i=0;i<NodeNum;i++) trist2_num_u[i] = tri_num_u[i] - (1 - p) * st2_num_u[i];
		}

		// Global sensitivity using noisy max degree --> sen_tri
		sen_tri = max_deg_ns;

		// Add Lap for each user
		for(i=0;i<NodeNum;i++){
			// Add Lap(sen_tri/Eps2ndTrSt) --> trist2_num_u
			trist2_num_u[i] += stats::rlaplace(0.0, sen_tri/Eps2ndTrSt, engine);
		}
	}

    // Empirical estimate --> tri_num_ns
    tri_num_ns = 0;
	for(i=0;i<NodeNum;i++) tri_num_ns += trist2_num_u[i];
	// Divide #triangles by 2p - 1 (= 1 - 2q)
	tri_num_ns /= (2.0 * p - 1.0);

	delete[] a_mat_ns;
	delete[] a_mat_del;
	free1D(tri_num_u);
	free1D(st2_num_u);
	free1D(trist2_num_u);
}

// Calculate the clustering-coefficient
double CalcClstCoef(double tri_num_ns, double st2_num_ns){
	double clst_ns;

    if(tri_num_ns < 0) tri_num_ns = 0;
    if(st2_num_ns < 0) st2_num_ns = 0;
    if(st2_num_ns == 0) clst_ns = 1.0;
    else clst_ns = 3.0 * tri_num_ns / st2_num_ns;
    if(clst_ns > 1.0) clst_ns = 1.0;
    else if(clst_ns < 0.0) clst_ns = 0.0;

	return clst_ns;
}

int main(int argc, char *argv[])
{
	int all_node_num;
	int triplet_num;
	int **node_order;
	map<int, int> *a_mat;			// adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg;									// degree
	long long tot_edge_num;
	long long tri_num, st2_num, st3_num, ed2_num, ed1_num, non_num;
	double clst;
	double tri_num_ns, sen_tri;
	double st2_num_ns, st3_num_ns, sen_st2, sen_st3;
	double clst_ns;
	double tri_re_ns, tri_l2_ns;
	double tri_re_ns_avg, tri_l2_ns_avg;
	double st2_re_ns, st2_l2_ns;
	double st2_re_ns_avg, st2_l2_ns_avg;
	double st3_re_ns, st3_l2_ns;
	double st3_re_ns_avg, st3_l2_ns_avg;
	double clst_re_ns, clst_l2_ns;
	double clst_re_ns_avg, clst_l2_ns_avg;
	int itr;
	int i, j, k, x;
	string outdir;
	string outfile;
	char s[1025], *str;
	char str_1[] = "1";
	FILE *fp;

	// Initialization of Mersennne Twister
	unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
	init_by_array(init, length);

	if (argc < 2) {
	    printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon (default: 1)] [NSType (default: 0)] [PosBias (default: 0)] [#itr (default: 1)] [alg (default: 0)] [Balloc (default: 1-1-1)])\n\n", argv[0]);
		printf("[EdgeFile]: Edge file\n");
		printf("[#nodes]: Number of nodes (-1: all)\n");
		printf("[epsilon]: Epsilon\n");
		printf("[NSType]: Noise type (0: Lap (#nodes), 1: Lap (max degree), 2: Lap (noisy max degree))\n");
		printf("[PosBias]: Positive bias (noisy max degree)\n");
		printf("[#itr]: Number of iterations\n");
		printf("[alg]: Algorithm (0: centralized, 1: non-interactive local (RR w/ emp), 2: non-interactive local (RR w/o emp), 3: interactive local)\n");
		printf("[Balloc]: Privacy budget allocation (alg=3): Eps1st-Eps2ndTr-Eps2ndSt\n");
		return -1;
	}

	EdgeFile = argv[1];

	NodeNum = -1;
	if (argc >= 3) NodeNum = atoi(argv[2]);
	Eps = 1.0;
	Eps_s = "1.0";
	if (argc >= 4){
		Eps = atof(argv[3]);
		Eps_s = argv[3];
	}
	// Epsilon for calculating the noisy max degree
	EpsNsMaxDeg = Eps / 10;
	NSType = 0;
	if (argc >= 5) NSType = atoi(argv[4]);
	PosBias = 0.0;
	PosBias_s = "0.0";
	if (argc >= 6){
		PosBias = atof(argv[5]) / EpsNsMaxDeg;
		PosBias_s = argv[5];
	}
	ItrNum = 1;
	if (argc >= 7) ItrNum = atoi(argv[6]);
	Alg = 0;
	if (argc >= 8) Alg = atoi(argv[7]);
	if (Alg < 0 || Alg > 4){
		printf("Error: incorect [Alg]\n");
		exit(-1);
	}
	if (Alg == 0 && NSType == 2){
		printf("Error: N/A (Alg = 0, NSType = 2)\n");
		exit(-1);
	}

	for(i=0;i<3;i++){
		Balloc[i] = 1.0;
		Balloc_s[i] = str_1;
	}
	if (argc >= 9){
		if((Balloc_s[0] = strtok(argv[8], "-")) == NULL){
			printf("Error: incorect [Balloc]\n");
			exit(-1);
		}
		Balloc[0] = atof(Balloc_s[0]);
		if((Balloc_s[1] = strtok(NULL, "-")) == NULL){
			printf("Error: incorect [Balloc]\n");
			exit(-1);
		}
		Balloc[1] = atof(Balloc_s[1]);
		if((Balloc_s[2] = strtok(NULL, "-")) == NULL){
			printf("Error: incorect [Balloc]\n");
			exit(-1);
		}
		Balloc[2] = atof(Balloc_s[2]);
	}

	// Total number of nodes --> all_node_num
	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<2;i++) fgets(s, 1024, fp);
	all_node_num = atoi(s);
	fclose(fp);

	// malloc
	malloc2D(&node_order, ItrNum, all_node_num);

	// Use all nodes
	if (NodeNum == -1){
		NodeNum = all_node_num;
		ItrNum = 1;
		for(j=0;j<NodeNum;j++) node_order[0][j] = j;
	}
	// Randomly generate the order of nodes --> node_order
	else{
		i = EdgeFile.find_last_of("/");
		outdir = EdgeFile.substr(0, i+1);
		outfile = outdir + "node-order_itr" + to_string(ItrNum) + ".csv";
		if(checkFileExistence(outfile)){
			fp = FileOpen(outfile, "r");
			for(j=0;j<all_node_num;j++){
				fgets(s, 1024, fp);
				strtok(s, ",");
				for(i=0;i<ItrNum;i++){
					node_order[i][j] = atoi(strtok(NULL, ","));
				}
			}
			fclose(fp);
		}
		else{
			for(i=0;i<ItrNum;i++){
				MakeRndPerm(node_order[i], all_node_num, all_node_num);
			}
			fp = FileOpen(outfile, "w");
			for(j=0;j<all_node_num;j++){
				fprintf(fp, "%d,", j);
				for(i=0;i<ItrNum;i++) fprintf(fp, "%d,", node_order[i][j]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
	}

	// #triplet --> triplet_num
	triplet_num = NodeNum*(NodeNum-1)*(NodeNum-2)/6;

	// Initialization
	malloc1D(&deg, NodeNum);
	tri_re_ns_avg = tri_l2_ns_avg = 0.0;
	st2_re_ns_avg = st2_l2_ns_avg = 0.0;
	st3_re_ns_avg = st3_l2_ns_avg = 0.0;
	clst_re_ns_avg = clst_l2_ns_avg = 0.0;

	// Output the header
	i = EdgeFile.find_last_of("/");
	outdir = EdgeFile.substr(0, i+1);
	for(i=0;i<3;i++){
		outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "_ns" + to_string(NSType) + "_pb" + PosBias_s + "_ba" + Balloc_s[0] + "-" + Balloc_s[1] + "-" + Balloc_s[2] + "_itr" + to_string(ItrNum) + ".csv";
		fp = FileOpen(outfile, "w");
		fprintf(fp, "#tri(true),#tri(est),#tri(rel-err),#tri(l2-loss),#2st(true),#2st(est),#2st(rel-err),#2st(l2-loss),#3st(true),#3st(est),#3st(rel-err),#3st(l2-loss),clst(true),clst(est),clst(rel-err),clst(l2-loss),sen_tri,sen_2st,sen_3st\n");
		fclose(fp);
	}

	// For each iteration
	for(itr=0;itr<ItrNum;itr++){
		// Initialization
		a_mat = new map<int, int>[NodeNum];

	    // Read edges from the edge file --> a_mat
    	ReadEdges(a_mat, node_order[itr]);

		// Degree --> deg
		for(i=0;i<NodeNum;i++) deg[i] = 0;
		for(i=0;i<NodeNum;i++){
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg[i] += 1;
		}

		// Total number of edges --> tot_edge_num
		tot_edge_num = 0;
		for(i=0;i<NodeNum;i++) tot_edge_num += (long long)deg[i];
		tot_edge_num /= 2;

		// #triangles --> tri_num
		tri_num = 0;
		for(i=0;i<NodeNum;i++){
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				j = aitr->first;
				if (i >= j) continue;
				for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
					k = aitr2->first;
					if (j >= k) continue;
					if(a_mat[j].count(k) > 0) tri_num++;
				}
			}
		}

	    // #2-stars, #3-stars --> st2_num, st3_num
		st2_num = st3_num = 0;
		for(i=0;i<NodeNum;i++){
			st2_num += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
			st3_num += ((long long)deg[i] * ((long long)deg[i]-1) * ((long long)deg[i]-2)) / 6;
		}

	    // clustering coefficient --> clst
		if(st2_num != 0) clst = 3.0 * (double)tri_num / (double)st2_num;
		else clst = 1.0;

 	    // #2-edges --> ed2_num
		ed2_num = st2_num - 3*tri_num;
		// #1-edge --> ed1_num
		ed1_num = (long long)tot_edge_num*(NodeNum-2) - 2*ed2_num - 3*tri_num;
	    // #none --> non_num
	    non_num = (long long)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num - ed2_num - ed1_num;

		/************************ Calculate sub-graph counts ************************/
		// Centralized
		if (Alg == 0){
			// Calculate #2-stars and #3-stars
        	CalcCentSt(st2_num, st3_num, a_mat, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
			CalcCentTri(tri_num, a_mat, deg, outfile, tri_num_ns, sen_tri);
		}
		// Non-interactive (1-round) local (RR w/ emp)
		else if (Alg == 1){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
//			CalcNLocTri(a_mat, deg, outfile, tri_num_ns, 1);
			if(NodeNum <= 10000) CalcNLocTri(a_mat, deg, outfile, tri_num_ns, 1);
			else tri_num_ns = 0.0;
			sen_tri = 0.0;
		}
		// Non-interactive (1-round) local (RR w/o emp)
		else if (Alg == 2){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
//			CalcNLocTri(a_mat, deg, outfile, tri_num_ns, 0);
			if(NodeNum <= 10000) CalcNLocTri(a_mat, deg, outfile, tri_num_ns, 0);
			else tri_num_ns = 0.0;
			sen_tri = 0.0;
		}
		// Interactive (2-rounds) local
		else if (Alg == 3){
			// Calculate #2-stars and #3-stars
        	CalcNLocSt(st2_num, st3_num, deg, outfile, st2_num_ns, st3_num_ns, sen_st2, sen_st3);
        	// Calculate #triangles
			CalcILocTri(a_mat, deg, outfile, tri_num_ns, sen_tri);
		}

		/******************** Calculate the cluster coefficient *********************/
		clst_ns = CalcClstCoef(tri_num_ns, st2_num_ns);

		/**************************** Evaluate the loss *****************************/
		// relative error --> tri_re_ns
		tri_re_ns = fabs(tri_num_ns - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);
		tri_re_ns_avg += tri_re_ns;
		// l2_loss --> tri_l2_ns
		tri_l2_ns = (tri_num_ns - (double)tri_num)*(tri_num_ns - (double)tri_num);
		tri_l2_ns_avg += tri_l2_ns;

		// relative error --> st2_re_ns
		st2_re_ns = fabs(st2_num_ns - (double)st2_num) / max((double)st2_num, 0.001 * NodeNum);
		st2_re_ns_avg += st2_re_ns;
		// l2_loss --> st2_l2_ns
		st2_l2_ns = (st2_num_ns - (double)st2_num)*(st2_num_ns - (double)st2_num);
		st2_l2_ns_avg += st2_l2_ns;

		// relative error --> st3_re_ns
		st3_re_ns = fabs(st3_num_ns - (double)st3_num) / max((double)st3_num, 0.001 * NodeNum);
		st3_re_ns_avg+= st3_re_ns;
		// l2_loss --> st3_l2_ns
		st3_l2_ns = (st3_num_ns - (double)st3_num)*(st3_num_ns - (double)st3_num);
		st3_l2_ns_avg += st3_l2_ns;

		// relative error --> clst_re_ns
		clst_re_ns = fabs(clst_ns - clst) / (double)clst;
		clst_re_ns_avg += clst_re_ns;
		// l2_loss --> clst_l2_ns
		clst_l2_ns = (clst_ns - clst)*(clst_ns - clst);
		clst_l2_ns_avg += clst_l2_ns;

		/**************************** Output the results ****************************/
		fp = FileOpen(outfile, "a");
		fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", 
		(double)tri_num, tri_num_ns, tri_re_ns, tri_l2_ns, 
		(double)st2_num, st2_num_ns, st2_re_ns, st2_l2_ns, 
		(double)st3_num, st3_num_ns, st3_re_ns, st3_l2_ns, 
		clst, clst_ns, clst_re_ns, clst_l2_ns,
		sen_tri, sen_st2, sen_st3);
		fclose(fp);

		delete[] a_mat;
	}

	/************************* Output the results (AVG) *************************/
	tri_re_ns_avg /= (double)ItrNum;
	tri_l2_ns_avg /= (double)ItrNum;
	st2_re_ns_avg /= (double)ItrNum;
	st2_l2_ns_avg /= (double)ItrNum;
	st3_re_ns_avg /= (double)ItrNum;
	st3_l2_ns_avg /= (double)ItrNum;
	clst_re_ns_avg /= (double)ItrNum;
	clst_l2_ns_avg /= (double)ItrNum;

	fp = FileOpen(outfile, "a");
	fprintf(fp, "function,AVG(rel-err),AVG(l2-loss)\n");
	fprintf(fp, "Triangles,%e,%e\n", tri_re_ns_avg, tri_l2_ns_avg);
	fprintf(fp, "2-stars,%e,%e\n", st2_re_ns_avg, st2_l2_ns_avg);
	fprintf(fp, "3-stars,%e,%e\n", st3_re_ns_avg, st3_l2_ns_avg);
	fprintf(fp, "Clst,%e,%e\n", clst_re_ns_avg, clst_l2_ns_avg);
	fclose(fp);

	// free
	free2D(node_order, ItrNum);
	free1D(deg);

	return 0;
}
