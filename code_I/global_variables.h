#ifndef VARIABLE_local_HEADER
#define VARIABLE_local_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ilcplex/cplex.h>
#include <set>
#include <queue>
#include <map>
#include <set>
#include <regex>

using namespace std;

#define INF_DIST 99999
#define MAX_SIZE 99999
#define epsilon 1e-9

#define to_file
#define relax_w
//#define relax_z
//#define relax_y
#define use_lp_solver

#define write_prob
//#define write_prob_call_back
#define print_solution
#define write_annotation
//#define print_star
//#define drawing

typedef struct inst_struct{
	char* input_file,*problem_type,*instance_type,*instance_name,*benders,*options,*posfile;
	int dimension, p, n_edges;
	int algorithm;
	double fixed_cost;
	int** c;
	vector<vector<int>> index_sorted;
	vector<vector<bool>> b;
	vector<vector<double>> pos;
	vector<vector<int>> fstar_cover;
	vector<vector<int>> fstar_facility;
	int** fs,**fs_facility;
	int** fs_pointer, **fs_pointer_facility;
	int* Js;
	double alpha;
	int* z_tilde;
	int* I_tilde;
	bool only_root = false;
	string instance_number;
	double timelimit;

	CPXENVptr env_MPIF, env_CPIF;
	CPXLPptr lp_MPIF, lp_CPIF;

	int lazy_cuts = 0;
	int lazy_cuts_connectivity = 0;
	int lazy_cuts_cost = 0;
	int lazy_cuts_coverage = 0;
	int user_cuts_coverage = 0;
	int user_cuts = 0;
	int user_cuts_cost = 0;

	double* rmatval, * obj, * lb, * ub, * rhs;
	double* x;
	double objval, best_lb, root_bound;
	char* c_type, * sense;
	char** colname;
	int* rmatbeg, * rmatind;

	int r, R;
	int status,ccnt,rcnt, nzcnt;

	int* y,*z,*v;
	int vv;
	int** w,**f,***fm;
	int teta;


	int* ySol, ** wSol;
	double** fSol;

	double opened_facilities = 0.0;

	double custom_solution_value = CPX_INFBOUND;
	bool use_custom_solution = false;

} instance;

#endif
