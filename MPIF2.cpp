#include "MPIF.h"
#include "MPIF2.h"

double xSol_v_2[MAX_SIZE];
double rmatval_v_2[MAX_SIZE];
int    rmatind_v_2[MAX_SIZE];
double rhs_v_2;

int CPXPUBLIC mylazycallback_NV_min(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p) {
(*useraction_p) = CPX_CALLBACK_DEFAULT;
instance* inst = (instance*)cbhandle;
//cout << "Lazy callback" << endl;

int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_v_2, 0, inst->ccnt - 1);
if (status != 0) {
	printf("cannot get the x\n");
	exit(-1);
}



//SEC
	//finding root subtour
vector <bool> neighbours_to_open(inst->dimension_f, false);
vector <bool> visited(inst->dimension_f, false);
vector <bool> root_component(inst->dimension_f, false);
vector <bool> not_reached(inst->dimension_f, false); //for plotting purposes
queue <int> remote_nodes;
queue <int> remote_nodes_2;
queue <int> qu;



vector<set<int>> C;
vector<set<int>> A;
C.push_back(set<int>());
A.push_back(set<int>());
//vector<bool>is_in_root_component(inst->dimension, false);

/*
for (int i = 0; i < inst->dimension; i++) {
	C.push_back(set<int>());
	A.push_back(set<int>());
}
*/
//C_0 A_0
//find C_i
fill(visited.begin(), visited.end(), false);
visited[0] = true;
C[0].insert(0);
qu.push(0);

while (!qu.empty()) {

	int jj = qu.front();
	qu.pop();

	for (auto neighbour : inst->fstar_facility[jj]) {
		if (!visited[neighbour]) {
			if (xSol_v_2[inst->y[neighbour]] > 0.999) {
				qu.push(neighbour);
				C[0].insert(neighbour);
			}
			else {
				A[0].insert(neighbour);
			}
			visited[neighbour] = true;
		}
	}
}
/*
for (int j = 1; j < inst->dimension; j++) {
	if (xSol[inst->y[j]] > 0.999) {
		//find C_i
		fill(visited.begin(), visited.end(), false);
		visited[j] = true;
		C[j].insert(j);
		qu.push(j);

		while (!qu.empty()) {

			int jj = qu.front();
			qu.pop();

			for (auto neighbour : inst->fstar_facility[jj]) {
				if (!visited[neighbour]) {
					if (neighbour == 0) {
						is_in_root_component[j] = true;
					}

					if (xSol[inst->y[neighbour]] > 0.999) {
						qu.push(neighbour);
						C[j].insert(neighbour);
					}
					else {
						A[j].insert(neighbour);
					}
					visited[neighbour] = true;
				}

			}
		}
	}
}
*/

//building node separator for each opened facility

bool network_is_feasible = true;
for (int j = 1; j < inst->dimension_f; j++) {
	if (xSol_v_2[inst->y[j]] > 0.999 && C[0].count(j) <= 0) {

		network_is_feasible = false;
		//reset
		fill(visited.begin(), visited.end(), false);
		int counter = 0;
		qu.push(j);
		visited[j] = false;

		while (!qu.empty()) {
			int i = qu.front();
			qu.pop();

			for (auto neighbour : inst->fstar_facility[i])
			{
				if (!visited[neighbour]) {

					if (A[0].count(neighbour) <= 0) {
						qu.push(neighbour);
					}
					else {
						//R intersection A_0
						rmatind_v_2[counter] = inst->y[neighbour];
						rmatval_v_2[counter] = 1.0;
						//cout << "+y" << neighbour;
						counter++;
					}

					visited[neighbour] = true;
				}
			}
		}

		rmatind_v_2[counter] = inst->y[j];
		rmatval_v_2[counter] = -1.0;
		//cout << "-y" << j;

		rhs_v_2 = 0.0;
		int nzcnt = counter + 1;



		//check if cut violated
		double lhs_check = 0, rhs_check = 0;

		for (int i = 0; i <= counter; i++) {
			lhs_check += xSol_v_2[rmatind_v_2[i]] * rmatval_v_2[i];
		}

		//cout << "lhs - rhs = " << lhs_check - rhs_check;
		if (lhs_check - rhs_check != 0.0) {
			//if(true){
			status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_v_2, 'G', rmatind_v_2, rmatval_v_2, 0);
			inst->lazy_cuts_connectivity++;
			if (status != 0) {
				printf("CPXcutcallbackadd\n");
				exit(-1);
			}

		}
	}
}

if(network_is_feasible) {

	for (int j = 0; j < inst->dimension_c; j++) {
		for (int k = 0; k < inst->dimension_f; k++) {
			if (!inst->b[j][inst->index_sorted[j][k]]) {
				double lhs = 0.0;
				for (int i = k; i >= 0; i--) {
					lhs += (inst->c_fc[inst->index_sorted[j][k]][j] - inst->c_fc[inst->index_sorted[j][i]][j]) * xSol_v_2[inst->y[inst->index_sorted[j][i]]];
				}
				lhs += xSol_v_2[inst->v[j]];
				double rhs = inst->c_fc[inst->index_sorted[j][k]][j];
				if (lhs < rhs - 0.00001) {
					int counter = 0;
					for (int i = k; i >= 0; i--) {
						rmatind_v_2[counter] = inst->y[inst->index_sorted[j][i]];
						rmatval_v_2[counter] = (inst->c_fc[inst->index_sorted[j][k]][j] - inst->c_fc[inst->index_sorted[j][i]][j]);
						counter++;
					}
					rmatind_v_2[counter] = inst->v[j];
					rmatval_v_2[counter] = 1.0;
					rhs_v_2 = inst->c_fc[inst->index_sorted[j][k]][j];

					int nzcnt = counter + 1;

					status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_v_2, 'G', rmatind_v_2, rmatval_v_2, 0);
					inst->lazy_cuts_cost++;
					if (status != 0) {
						printf("CPXcutcallbackadd\n");
						exit(-1);
					}
					inst->b[j][inst->index_sorted[j][k]] = true;

					//skipping to the next j
					break;
				}
			}
		}
	}

}


#ifdef write_prob_call_back
double per_test;
CPXLPptr lp_test;
CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
CPXwriteprob(env, lp_test, "prob.lp", NULL);
cout << "NODE LP WRITTEN Lazy callback\n\n\n";
cin.get();
#endif

(*useraction_p) = CPX_CALLBACK_SET;
return 0;
}

int CPXPUBLIC myusercallback_NV(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p) {

	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle;

	int mythread = 0;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
	int nodecnt; 
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, (void*)(&nodecnt));

	if (nodecnt == 0 && mythread == 0)		// root node best bound (thread safe) 
	{
		CPXLPptr nodelp; CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp);
		CPXgetobjval(env, nodelp, &(inst->root_bound));
		cout <<"Root bound: "<< inst->root_bound<<endl;
	}

	
	
	vector<int> rmatind;
	vector<double> rmatval;

	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_v_2, 0, inst->ccnt - 1);
	if (status != 0) {
		printf("cannot get the x\n");
		exit(-1);
	}

	if (wherefrom != CPX_CALLBACK_MIP_CUT_LAST) return 0;

	for (int j = 0; j < inst->dimension_c; j++) {
	
		for (int k = 0; k < inst->dimension_f; k++) {
			if (!inst->b[j][inst->index_sorted[j][k]]) {
				double lhs = 0.0;
				for (int i = k; i >= 0; i--) {
					lhs += (inst->c_fc[inst->index_sorted[j][k]][j] - inst->c_fc[inst->index_sorted[j][i]][j]) * xSol_v_2[inst->y[inst->index_sorted[j][i]]];
				}
				lhs += xSol_v_2[inst->v[j]];
				double rhs = inst->c_fc[inst->index_sorted[j][k]][j];
				if (lhs < rhs - 0.00001) {
					int counter = 0;
					for (int i = k; i >= 0; i--) {
						rmatind_v_2[counter] = inst->y[inst->index_sorted[j][i]];
						rmatval_v_2[counter] = (inst->c_fc[inst->index_sorted[j][k]][j] - inst->c_fc[inst->index_sorted[j][i]][j]);
						counter++;
					}
					rmatind_v_2[counter] = inst->v[j];
					rmatval_v_2[counter] = 1.0;
					rhs_v_2 = inst->c_fc[inst->index_sorted[j][k]][j];

					int nzcnt = counter + 1;

					//cout << "User callback" << endl;
					status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_v_2, 'G', rmatind_v_2, rmatval_v_2, 0);
					inst->user_cuts_cost++;
					if (status != 0) {
						printf("CPXcutcallbackadd\n");
						exit(-1);
					}
					inst->b[j][inst->index_sorted[j][k]] = true;


					//skipping to the next j
					break;
				}
			}
		}


		/*
		int k = 0;
		double sum_k = xSol_v_2[inst->y[inst->index_sorted[j][k]]];

		while (sum_k < 1) {
			k++;
			sum_k += xSol_v_2[inst->y[inst->index_sorted[j][k]]];
		}

		rmatind.push_back(inst->v[j]);
		rmatval.push_back(1.0);

		for (int i = 0; i < k; i++) {

			rmatind.push_back(inst->y[inst->index_sorted[j][i]]);
			rmatval.push_back(inst->c[j][k] - inst->c[j][inst->index_sorted[j][i]]);

		}		
	
		int nzcnt = rmatind.size();
		int rhs_user = inst->c[j][k];

		status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_user, 'G', rmatind.data(), rmatval.data(), 0);
		if (status != 0) {
			printf("CPXcutcallbackadd\n");
			exit(-1);
		}
		else {
			inst->user_cuts_cost++;
		}

		rmatind.clear();
		rmatval.clear();
	*/
	}	
	
#ifdef write_prob_call_back
		double per_test;
		CPXLPptr lp_test;
		CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
		CPXwriteprob(env, lp_test, "prob.lp", NULL);
		cout << "NODE LP WRITTEN User callback\n\n\n";
		cin.get();
#endif
	
		(*useraction_p) = CPX_CALLBACK_SET;
	

	return 0;
}

void add_cardinal_constraint(instance* inst) {

	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	//inst->rhs[0] = 2.0;
	inst->rhs[0] = 3.0;

	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'G';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension_f;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));


	for (int j = 0; j < inst->dimension_f; j++) {
		inst->rmatind[j] = inst->y[j];
		inst->rmatval[j] = 1.0;
	}

	inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
	if (inst->status != 0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);

}

double calculate_only_root_case(instance* inst)
{

	double cost = 0.0;
	for (int i = 0; i < inst->dimension_c; i++) {
		cost += inst->demand[i] * inst->c_fc[0][i];
	}
	
	return cost;
}

double calculate_only_2(instance* inst)
{
	
	double best_cost = CPX_INFBOUND;
	
	if (inst->fstar_facility[0].size() > 1)
		for (auto second_facility : inst->fstar_facility[0])
			if (second_facility != 0)
			{
				double cost = inst->fixed_cost_ind[second_facility];
				for (int i = 0; i < inst->dimension_c; i++) {
					cost += min(inst->demand[i] * inst->c_fc[0][i], inst->demand[i] * inst->c_fc[second_facility][i]);
				}

				if (best_cost > cost) best_cost = cost;
			}


	return best_cost;
}

void build_model_MPIF_NV(instance* inst) {
	cout << "Only y and v variables" << endl;

	//CPLEX environment
	inst->env_MPIF = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//CPLEX problem
	inst->lp_MPIF = CPXcreateprob(inst->env_MPIF, &(inst->status), "MPIF");
	if (inst->status != 0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}

	//adding variables
	inst->ccnt = inst->dimension_f + inst->dimension_c;
	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(10, sizeof(char)); }
	inst->y = (int*)calloc(inst->dimension_f, sizeof(int));
	inst->v = (int*)calloc(inst->dimension_c, sizeof(int));

	int counter = 0;
	for (int j = 0; j < inst->dimension_f; j++) {

		if (inst->fixed_cost_ind[j] == 1)
			inst->obj[counter] = 0.0;
		else
			inst->obj[counter] = inst->fixed_cost_ind[j];
		inst->lb[counter] = 0.0;
		inst->ub[counter] = 1.0;
#ifdef relax_y
		inst->c_type[counter] = 'C';
#else
		inst->c_type[counter] = 'B';
#endif // relax_y
		sprintf(inst->colname[counter], "y.%d", j);
		inst->y[j] = counter;
		counter++;
	}

	// y_0 = 1
	inst->lb[inst->y[0]] = 1.0;
	inst->obj[inst->y[0]] = 0.0;

	for (int j = 0; j < inst->dimension_c; j++) {
		inst->obj[counter] = inst->demand[j];
		inst->lb[counter] = 0.0;
		inst->ub[counter] = CPX_INFBOUND;
		inst->c_type[counter] = 'C';
		sprintf(inst->colname[counter], "v.%d", j);
		inst->v[j] = counter;
		counter++;
	}

	//finding root subtour
	vector <bool> reachable(inst->dimension_f, false);
	vector <bool> visited(inst->dimension_f, false);
	queue <int> qu;

	//root node subtour
	visited[0] = true;
	reachable[0] = true;
	qu.push(0);

	while (!qu.empty()) {
		int i = qu.front();
		qu.pop();

		//for (int j = inst->fs_pointer_facility[i][1]; j <= inst->fs_pointer_facility[i][2] && j > -1; j++) 
		for (auto neighbour : inst->fstar_facility[i]) {
			//int neighbour = inst->fs_facility[j][1];
			if (!visited[neighbour]) {
				qu.push(neighbour);
				reachable[neighbour] = true;
				visited[neighbour] = true;
			}
		}
	}

	for (int k = 0; k < inst->dimension_f; k++) {
		if (!reachable[k]) {
			inst->ub[inst->y[k]] = 0.0;
		}
	}

	inst->status = CPXnewcols(inst->env_MPIF, inst->lp_MPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
	if (inst->status != 0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

	//optimization sense - minimization
	CPXchgobjsen(inst->env_MPIF, inst->lp_MPIF, CPX_MIN);

	//adding constraints
	if (inst->fixed_cost_ind[0] == 1) {
		inst->rcnt = 1;
		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->rhs[0] = B;
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
		inst->sense[0] = 'L';
		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatbeg[0] = 0;
		inst->nzcnt = inst->dimension_f;

		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->dimension_f; i++)
		{

			inst->rmatval[i] = 1.0;
			inst->rmatind[i] = inst->y[i];

		}

		inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			std::printf("error in CPXaddrows\n");
			exit(-1);
		}

		free(inst->rmatval);
		free(inst->rmatind);
	}
	
	add_cardinal_constraint(inst);

	//sort distances
	calculate_sorted_c(inst);

	//initial v constraints
	initialize_v(inst);


#ifdef write_prob
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	inst->status = CPXwriteprob(inst->env_MPIF, inst->lp_MPIF, "MPIF_project_out_w.lp", NULL);
	if (inst->status != 0) {
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//exit(-1);
#endif
}

void solve_model_MPIF_NV(instance* inst) {
	CPXsetintparam(inst->env_MPIF, CPX_PARAM_SCRIND, CPX_ON);
	CPXsetintparam(inst->env_MPIF, CPX_PARAM_THREADS, 1);
	inst->status = CPXsetdblparam(inst->env_MPIF, CPX_PARAM_TILIM, inst->timelimit);
	if (inst->status)
	{
		printf("error for CPX_PARAM_EPRHS\n");
	}

	CPXsetintparam(inst->env_MPIF, CPX_PARAM_MIPCBREDLP, CPX_OFF);        // let MIP callbacks work on the original model
	CPXsetintparam(inst->env_MPIF, CPX_PARAM_PRELINEAR, CPX_OFF);         // assure linear mappings between the presolved and original models
	CPXsetintparam(inst->env_MPIF, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);

	//CPXsetintparam(inst->env_MPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);	
	   
	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_MPIF, mylazycallback_NV_min, inst);
	if (inst->status)
	{
		printf("error for CPXsetlazyconstraintcallbackfunc\n");
	}

	
	inst->status = CPXsetusercutcallbackfunc(inst->env_MPIF, myusercallback_NV, inst);
	if (inst->status)
	{
		printf("error for CPXsetuserconstraintcallbackfunc\n");
	}
	

	inst->lazy_cuts_connectivity = 0;
	inst->lazy_cuts_cost = 0;
	inst->user_cuts_cost = 0;

	double solution_time;


	clock_t time_start = clock();

	cout << "\nCPXmipopt:\n";

	double root_solution = calculate_only_root_case(inst);
	double two_solution = calculate_only_2(inst);

	double cut_up = min(root_solution, two_solution);

	CPXsetdblparam(inst->env_MPIF, CPX_PARAM_CUTUP, cut_up);

	inst->status = CPXmipopt(inst->env_MPIF, inst->lp_MPIF);
	if (inst->status != 0)
	{
		printf("error in CPXmipopt\n");
		//exit(-1);
	}

	clock_t time_end = clock();

	solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	inst->status = CPXgetmipobjval(inst->env_MPIF, inst->lp_MPIF, &(inst->objval));
	if (inst->status != 0)
	{
		printf("error in CPXgetmipobjval\n");
	}

	int lpstat = CPXgetstat(inst->env_MPIF, inst->lp_MPIF);

	if (((lpstat == CPXMIP_OPTIMAL) ||
		(lpstat == CPXMIP_OPTIMAL_INFEAS) ||
		//( lpstat ==  CPXMIP_OPTIMAL_RELAXED ) ||
		(lpstat == CPXMIP_OPTIMAL_TOL) ||
		(lpstat == CPXMIP_TIME_LIM_FEAS)) && (inst->objval < min(root_solution, two_solution))) {

		inst->x = (double*)calloc(inst->ccnt, sizeof(double));
		inst->status = CPXgetmipx(inst->env_MPIF, inst->lp_MPIF, inst->x, 0, inst->ccnt - 1);

		inst->opened_facilities = 0;
		for (int i = 0; i < inst->dimension_f; i++) {
			inst->opened_facilities += inst->x[inst->y[i]];
		}
		inst->custom_solution_value = inst->objval;
	}
	else {
		inst->use_custom_solution = true;

		if (root_solution <= two_solution) {
			inst->custom_solution_value = root_solution;
			inst->opened_facilities = 1;
		}
		else {
			inst->custom_solution_value = two_solution;
			inst->opened_facilities = 2;
		}
	}
	

	print_solution_MPIF(inst);
	write_log_MPIF(inst, solution_time);
	
	printf("\n\nMIP solution value ->\t\%f", inst->objval);
	cout << "Time elapsed: " << solution_time;
	
	cout << "Lazy cuts connectivity: " << inst->lazy_cuts_connectivity << "  Lazy cuts cost: " << inst->lazy_cuts_cost<< " User cuts cost"  <<inst->user_cuts_cost<<endl;
}
