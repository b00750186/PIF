#include "CPIF.h"
#include "CPIF2.h"

double xSol_CPIF[MAX_SIZE];
double obj_value_CPIF;
double rmatval_CPIF[MAX_SIZE];
int    rmatind_CPIF[MAX_SIZE];
double rhs_CPIF;

double xSol_CPIF_user[MAX_SIZE];
double obj_value_CPIF_user;
double rmatval_CPIF_user[MAX_SIZE];
int    rmatind_CPIF_user[MAX_SIZE];
double rhs_CPIF_user;

double xSol_CPIF_SEC[MAX_SIZE];
double obj_value_CPIF_SEC;
double rmatval_CPIF_SEC[MAX_SIZE];
int    rmatind_CPIF_SEC[MAX_SIZE];
double rhs_CPIF_SEC;

int CPXPUBLIC mydummycallback_CPIF(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	return 0;
}

int CPXPUBLIC mylazycallback_NZ(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p) {
	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle;
	//cout << "Lazy callback" << endl;

	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_CPIF, 0, inst->ccnt - 1);
	if (status != 0) {
		cout << "error getting x at node" << endl;
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


	//root node subtour
	visited[0] = true;
	root_component[0] = true;
	qu.push(0);

	while (!qu.empty()) {
		int i = qu.front();
		qu.pop();

		for (auto neighbour : inst->fstar_facility[i])
		{
			if (!visited[neighbour]) {
				if (xSol_CPIF[inst->y[neighbour]] > 0.999) {
					qu.push(neighbour);
					root_component[neighbour] = true;
				}
				else {
					neighbours_to_open[neighbour] = true;
				}
				visited[neighbour] = true;
			}
		}
	}

	//opened facilities not connected to the root
	for (int j = 0; j < inst->dimension_f; j++) {
		if (xSol_CPIF[inst->y[j]] > 0.999 && !root_component[j]) {
			remote_nodes.push(j);
			remote_nodes_2.push(j);
		}
	}

	//cut for each opened facility not connected to root component

	int counter = 0;
	for (int j = 0; j < inst->dimension_f; j++) {
		if (neighbours_to_open[j]) {
			rmatind_CPIF[counter] = inst->y[j];
			rmatval_CPIF[counter] = 1.0;
			counter++;
		}
	}
	int nzcnt = counter + 1;
	rhs_CPIF = 0;


	//connection from the root to remote facilities
	while (!remote_nodes.empty()) {
		rmatind_CPIF[counter] = inst->y[remote_nodes.front()];
		rmatval_CPIF[counter] = -1.0;
		not_reached[remote_nodes.front()] = true;
		remote_nodes.pop();

		status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF, 'G', rmatind_CPIF, rmatval_CPIF, 0);
		inst->lazy_cuts_connectivity++;
		if (status != 0) {
			cout << "CPXcutcallbackadd" << endl;
			exit(-1);
		}
	}

	//connection from remote facilities to the root
	//reset flags
	vector<bool> visited_opened(inst->dimension_f, false);
	queue<int> opened_facility_in_connected_component;
	vector<bool> neighbours_of_conn_comp(inst->dimension_f, false);
	vector<bool> visited_in_comp(inst->dimension_f, false);

	while (!remote_nodes_2.empty()) {
		int k = remote_nodes_2.front();
		remote_nodes_2.pop();


		if (!visited_opened[k]) {

			visited_opened[k] = true;

			fill(neighbours_to_open.begin(), neighbours_to_open.end(), false);
			fill(visited_in_comp.begin(), visited_in_comp.end(), false);

			qu.push(k);
			opened_facility_in_connected_component.push(k);

			while (!qu.empty()) {
				int i = qu.front();
				qu.pop();

				for (auto neighbour : inst->fstar_facility[i]) {
					if (!visited_in_comp[neighbour]) {
						visited_in_comp[neighbour] = true;

						if (xSol_CPIF[inst->y[neighbour]] > 0.999) {
							qu.push(neighbour);

							//queue for cut generation
							opened_facility_in_connected_component.push(neighbour);

							//flag for global cycle for all opened facilities not connected to the root
							visited_opened[neighbour] = true;
						}
						else {
							neighbours_to_open[neighbour] = true;
						}

					}
				}
			}

			//cut for each opened facility in the connected component
			counter = 0;
			for (int j = 0; j < inst->dimension_f; j++) {
				if (neighbours_to_open[j]) {
					rmatind_CPIF[counter] = inst->y[j];
					rmatval_CPIF[counter] = 1.0;
					counter++;
				}
			}
			int nzcnt = counter + 1;
			rhs_CPIF = 0;

			while (!opened_facility_in_connected_component.empty()) {
				rmatind_CPIF[counter] = inst->y[opened_facility_in_connected_component.front()];
				rmatval_CPIF[counter] = -1.0;

				opened_facility_in_connected_component.pop();


				status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF, 'G', rmatind_CPIF, rmatval_CPIF, 0);
				inst->lazy_cuts_connectivity++;
				if (status != 0) {
					cout << "CPXcutcallbackadd" << endl;
					exit(-1);
				}
			}

		}


	}


#ifdef write_prob_call_back
	double per_test;
	CPXLPptr lp_test;
	CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
	CPXwriteprob(env, lp_test, "prob_NZ.lp", NULL);
	cout << "Node LP Lazy callback is written" << endl;
	cin.get();
#endif

	(*useraction_p) = CPX_CALLBACK_SET;
	return 0;
}

//using only minimal separator
int CPXPUBLIC mylazycallback_NZ_min(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p) {
	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle;
	//cout << "Lazy callback" << endl;

	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_CPIF, 0, inst->ccnt - 1);
	if (status != 0) {
		cout << "error getting x at node" << endl;
		exit(-1);
	}

	int mythread = 0;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
	int nodecnt;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, (void*)(&nodecnt));

	if (nodecnt == 0 && mythread == 0)		// root node best bound (thread safe)
	{
		CPXLPptr nodelp; CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp);
		CPXgetobjval(env, nodelp, &(inst->root_bound));
		cout << "Root bound: " << inst->root_bound << endl;
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
				if (xSol_CPIF[inst->y[neighbour]] > 0.999) {
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

	for (int j = 1; j < inst->dimension_f; j++) {
		if (xSol_CPIF[inst->y[j]] > 0.999 && C[0].count(j) <= 0) {

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
							rmatind_CPIF[counter] = inst->y[neighbour];
							rmatval_CPIF[counter] = 1.0;
							//cout << "+y" << neighbour;
							counter++;
						}

						visited[neighbour] = true;
					}
				}
			}

			rmatind_CPIF[counter] = inst->y[j];
			rmatval_CPIF[counter] = -1.0;
			//cout << "-y" << j;

			rhs_CPIF = 0.0;
			int nzcnt = counter + 1;



			//check if cut violated
			double lhs_check = 0, rhs_check = 0;

			for (int i = 0; i <= counter; i++) {
				lhs_check += xSol_CPIF[rmatind_CPIF[i]] * rmatval_CPIF[i];
			}

			//cout << "lhs - rhs = " << lhs_check - rhs_check;
			if (lhs_check - rhs_check != 0.0) {
				//if(true){
				status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF, 'G', rmatind_CPIF, rmatval_CPIF, 0);
				inst->lazy_cuts_connectivity++;
				if (status != 0) {
					cout << "CPXcutcallbackadd" << endl;
					exit(-1);
				}

			}


		}
	}



#ifdef write_prob_call_back
	double per_test;
	CPXLPptr lp_test;
	CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
	CPXwriteprob(env, lp_test, "prob_NZ_min.lp", NULL);
	cout << "Node LP Lazy callback is written" << endl;
	//cin.get();
#endif

	(*useraction_p) = CPX_CALLBACK_SET;
	return 0;
}

void add_singleflow_constraints_CPIF(instance* inst) {

	cout <<"Adding single commodity constraints" << endl;
	//\sum_{i \in I \setminus \{0\}} f_{0 i}=\sum_{i \in I \setminus 0} y_{i}-1
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	vector<int> rmatind = vector<int>();
	vector<double> rmatval = vector<double>();

	int counter = 0;
	for (int i : inst->fstar_facility[0]) if (i != 0) {
		rmatind.push_back(inst->f[0][i]);
		rmatval.push_back(1.0);
		counter++;
		rmatind.push_back(inst->f[i][0]);
		rmatval.push_back(-1.0);
		counter++;
	}

	for (int i = 1; i < inst->dimension_f; i++) {
		rmatind.push_back(inst->y[i]);
		rmatval.push_back(-1.0);
		counter++;
	}

	inst->nzcnt = counter;
	inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL);
	if (inst->status != 0)
	{
		cout << "error in CPXaddrows" << endl;
		exit(-1);
	}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);
	rmatind.clear();
	rmatval.clear();
	// \sum_{i \in I \setminus \{j\}} f_{i j}-\sum_{i \in I \setminus \{j\}} f_{j i} = y_{j}
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;



	for (int j = 1; j < inst->dimension_f; j++) {
		counter = 0;
		if (inst->fstar_facility[j].size() > 1) {
			for (int i : inst->fstar_facility[j]) {
				if (i != j) {
					rmatind.push_back(inst->f[i][j]);
					rmatval.push_back(1.0);
					counter++;

					rmatind.push_back(inst->f[j][i]);
					rmatval.push_back(-1.0);
					counter++;
				}
			}
		}

		rmatind.push_back(inst->y[j]);
		rmatval.push_back(-1.0);
		inst->nzcnt = counter + 1;
		inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL);
		if (inst->status != 0)
		{
			cout << "error in CPXaddrows" << endl;
			exit(-1);
		}
		rmatind.clear();
		rmatval.clear();
	}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	//f_{i j} + f_{j i} \leq M y_{i}
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = 3;
	rmatind.push_back(0);
	rmatind.push_back(0);
	rmatind.push_back(0);
	rmatval.push_back(0);
	rmatval.push_back(0);
	rmatval.push_back(0);

	double M = inst->dimension_f - 1;

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j : inst->fstar_facility[i]) if (i != j) {

			counter = 0;

			rmatind[0] = inst->f[i][j];
			rmatval[0] = 1.0;
			rmatind[1] = inst->f[j][i];
			rmatval[1] = 1.0;
			rmatind[2] = inst->y[i];
			rmatval[2] = -M;
			inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL);
			if (inst->status != 0)
			{
				cout << "error in CPXaddrows" << endl;
				exit(-1);
			}
		}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);
}

void add_multiflow_constraints_CPIF(instance* inst) {

	cout << "Adding multiflow constraints 1" << endl;

	vector<int> rmatind;
	vector<double> rmatval;
	for (int i = 1; i < inst->dimension_f; i++)
		for (int k = 0; k < inst->dimension_f; k++)
			if (inst->fstar_facility[k].size() > 0) {

				inst->rcnt = 1;
				inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
				inst->rhs[0] = 0.0;
				inst->sense = (char*)calloc(inst->rcnt, sizeof(char));
				inst->sense[0] = 'E';
				inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
				inst->rmatbeg[0] = 0;

				for (int j : inst->fstar_facility[k]) {
					if (k != j && inst->c_f[k][j] <= inst->r) {
						rmatind.push_back(inst->fm[i][j][k]);
						rmatval.push_back(1.0);

						rmatind.push_back(inst->fm[i][k][j]);
						rmatval.push_back(-1.0);
					}
				}

				//y
				if (i == k) {
					rmatind.push_back(inst->y[i]);
					rmatval.push_back(-1.0);
				}
				else if (k == 0) {
					rmatind.push_back(inst->y[i]);
					rmatval.push_back(1.0);
				}

				inst->nzcnt = rmatind.size();

				inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL);
				if (inst->status != 0)
				{
					cout << "error in CPXaddrows" << endl;
					exit(-1);
				}

				rmatind.clear();
				rmatval.clear();

				free(inst->rmatbeg);
				free(inst->rhs);
				free(inst->sense);
			}

	//flow between facilities j and k
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = 3;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	cout << "Adding multiflow constraints 2" << endl;
	for (int i = 1; i < inst->dimension_f; i++)
		for (int k = 0; k < inst->dimension_f; k++) {
			for (int j : inst->fstar_facility[k]) {
				{
					if (j != k && inst->c_f[j][k] <= inst->c_f[j][k]) {
						inst->rmatind[0] = inst->fm[i][j][k];
						inst->rmatval[0] = 1.0;
						inst->rmatind[1] = inst->fm[i][k][j];
						inst->rmatval[1] = 1.0;
						inst->rmatind[2] = inst->y[k];
						inst->rmatval[2] = -1.0;
						inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
						if (inst->status != 0)
						{
							cout << "error in CPXaddrows" << endl;
							exit(-1);
						}
					}
				}
			}

		}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);

}

void add_pmed_constraint_CPIF(instance* inst) {
	// p median

	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = inst->p;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension_f;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));


	for (int j = 0; j < inst->dimension_f; j++) {
		inst->rmatind[j] = inst->y[j];
		inst->rmatval[j] = 1.0;
	}

	inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
	if (inst->status != 0)
	{
		cout << "error in CPXaddrows" << endl;
		exit(-1);
	}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);

}

//single commmodity flow model
void build_model_CPIF_SZ(instance* inst)
{
	//CPLEX environment
	inst->env_CPIF = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		cout << "error in CPXopenCPLEX" << endl;
		exit(-1);
	}

	//CPLEX problem
	inst->lp_CPIF = CPXcreateprob(inst->env_CPIF, &(inst->status), "CPIF");
	if (inst->status != 0)
	{
		cout << "error in CPXcreateprob" << endl;
		exit(-1);
	}

	//adding variables
	inst->ccnt = inst->dimension_f * inst->dimension_f + inst->dimension_c + inst->dimension_f;
	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(15, sizeof(char)); }

	inst->y = (int*)calloc(inst->dimension_f, sizeof(int));
	inst->z = (int*)calloc(inst->dimension_c, sizeof(int));
	inst->f = (int**)calloc(inst->dimension_f, sizeof(int*));
	for (int i = 0; i < inst->dimension_f; i++) { inst->f[i] = (int*)calloc(inst->dimension_f, sizeof(int)); }

	int counter = 0;
	for (int j = 0; j < inst->dimension_f; j++) {
		inst->obj[counter] = inst->fixed_cost_ind[j] * inst->alpha;
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
		inst->obj[counter] = -inst->demand[j];
		inst->lb[counter] = 0.0;
		inst->ub[counter] = 1.0;

#ifdef  relax_z
		inst->c_type[counter] = 'C';
#else
		inst->c_type[counter] = 'B';
#endif //  relax_z
		sprintf(inst->colname[counter], "z.%d", j);
		inst->z[j] = counter;
		counter++;
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++)
			if (j != i && inst->c_f[j][i] <= inst->r)
		{
			inst->obj[counter] = 0.0;
			inst->lb[counter] = 0.0;

			inst->ub[counter] = CPX_INFBOUND;

			inst->c_type[counter] = 'C';
			sprintf(inst->colname[counter], "f.%d.%d", i, j);
			inst->f[i][j] = counter;
			counter++;
		}

	inst->ccnt = counter;
#ifdef use_lp_solver
	inst->status = CPXnewcols(inst->env_CPIF, inst->lp_CPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, NULL, inst->colname);
#else
	inst->status = CPXnewcols(inst->env_CPIF, inst->lp_CPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
#endif

	if (inst->status != 0)
	{
		cout << "error in CPXnewcols" << endl;
		exit(-1);
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

	//optimization sense - minimization
	CPXchgobjsen(inst->env_CPIF, inst->lp_CPIF, CPX_MIN);

	//adding offset
	CPXchgobjoffset(inst->env_CPIF, inst->lp_CPIF, inst->total_demand);

	//adding constraints
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

		inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout << "error in CPXaddrows" << endl;
			exit(-1);
		}

		free(inst->rmatval);
		free(inst->rmatind);
	}

	// z_k <= sum a_ki y_i
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;

	vector<int> rmatind;
	vector<double> rmatval;

	cout << "Adding coverage constraints"<<endl;
	for (int k = 0; k < inst->dimension_c; k++) {
		counter = 0;
		for (int i:inst->fstar_cover[k]) {
			rmatind.push_back(inst->y[i]);
			rmatval.push_back(-1.0);
			counter++;
		}

		rmatind.push_back(inst->z[k]);
		rmatval.push_back(1.0);
		inst->nzcnt = counter + 1;
		inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL);
		if (inst->status != 0)
		{
			cout << "error in CPXaddrows" << endl;
			exit(-1);
		}
		rmatind.clear();
		rmatval.clear();
	}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	add_singleflow_constraints_CPIF(inst);

#ifdef write_prob
	inst->status = CPXwriteprob(inst->env_CPIF, inst->lp_CPIF, "CPIFSZ.lp", NULL);
	if (inst->status != 0) {
		cout << "error in CPXwriteprob" << endl;
		exit(-1);
	}
	//cin.get();
#endif
}

//node separator
void build_model_CPIF_NZ(instance* inst)
{
	//CPLEX environment
	inst->env_CPIF = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		cout << "error in CPXopenCPLEX" << endl;
		exit(-1);
	}

	//CPLEX problem
	inst->lp_CPIF = CPXcreateprob(inst->env_CPIF, &(inst->status), "CPIF_NZ");
	if (inst->status != 0)
	{
		cout << "error in CPXcreateprob" << endl;
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
	inst->z = (int*)calloc(inst->dimension_c, sizeof(int));

	int counter = 0;
	for (int j = 0; j < inst->dimension_f; j++) {
		inst->obj[counter] = inst->fixed_cost_ind[j] * inst->alpha;
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
		inst->obj[counter] = -inst->demand[j];
		inst->lb[counter] = 0.0;
		inst->ub[counter] = 1.0;

#ifdef  relax_z
		inst->c_type[counter] = 'C';
#else
		inst->c_type[counter] = 'B';
#endif //  relax_z
		sprintf(inst->colname[counter], "z.%d", j);
		inst->z[j] = counter;
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

		for (auto neighbour: inst->fstar_facility[i]) {
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


#ifdef use_lp_solver
	inst->status = CPXnewcols(inst->env_CPIF, inst->lp_CPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, NULL, inst->colname);
#else
	inst->status = CPXnewcols(inst->env_CPIF, inst->lp_CPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
#endif

	if (inst->status != 0)
	{
		cout << "error in CPXnewcols" << endl;
		exit(-1);
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

	//optimization sense - minimization
	CPXchgobjsen(inst->env_CPIF, inst->lp_CPIF, CPX_MIN);


	//adding offset
	CPXchgobjoffset(inst->env_CPIF, inst->lp_CPIF, inst->total_demand);
	//adding constraints
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

		inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout << "error in CPXaddrows" << endl;
			exit(-1);
		}

		free(inst->rmatval);
		free(inst->rmatind);
	}
	// z_k <= sum a_ki y_i
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension_f + 1;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	for (int k = 0; k < inst->dimension_c; k++) {
		counter = 0;
		for (int i = 0; i < inst->dimension_f; i++) {
			inst->rmatind[counter] = inst->y[i];
			if (inst->c_fc[i][k] > inst->R)
				inst->rmatval[counter] = 0.0;
			else
				inst->rmatval[counter] = -1.0;
			counter++;
		}

		inst->rmatind[counter] = inst->z[k];
		inst->rmatval[counter] = 1.0;

		inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout << "error in CPXaddrows" << endl;
			exit(-1);
		}
	}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);

	add_cardinal_constraint_CPIF(inst);

#ifdef write_prob

	inst->status = CPXwriteprob(inst->env_CPIF, inst->lp_CPIF, "CPIFNZ.lp", NULL);
	if (inst->status != 0) {
		cout << "error in CPXwriteprob" << endl;
		exit(-1);
	}

#endif
}

//multi commmodity flow model
void build_model_CPIF_MZ(instance* inst)
{
	//CPLEX environment
	inst->env_CPIF = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		cout << "error" << endl;
		exit(-1);
	}

	//CPLEX problem
	inst->lp_CPIF = CPXcreateprob(inst->env_CPIF, &(inst->status), "CPIF");
	if (inst->status != 0)
	{
		cout << "error" << endl;
		exit(-1);
	}

	//adding variables
	inst->ccnt = inst->dimension_f * inst->dimension_f * inst->dimension_f + inst->dimension_f + inst->dimension_c;
	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(25, sizeof(char)); }

	inst->y = (int*)calloc(inst->dimension_f, sizeof(int));
	inst->z = (int*)calloc(inst->dimension_c, sizeof(int));

	inst->fm = (int***)calloc(inst->dimension_f, sizeof(int**));
	for (int i = 0; i < inst->dimension_f; i++) { inst->fm[i] = (int**)calloc(inst->dimension_f, sizeof(int*)); }
	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++)
		{
			inst->fm[i][j] = (int*)calloc(inst->dimension_f, sizeof(int));
		}

	int counter = 0;
	for (int j = 0; j < inst->dimension_f; j++) {

		inst->obj[counter] = inst->fixed_cost_ind[j] * inst->alpha;
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

		inst->obj[counter] = -inst->demand[j];
		inst->lb[counter] = 0.0;
		inst->ub[counter] = 1.0;

#ifdef  relax_z
		inst->c_type[counter] = 'C';
#else
		inst->c_type[counter] = 'B';
#endif //  relax_z
		sprintf(inst->colname[counter], "z.%d", j);
		inst->z[j] = counter;
		counter++;
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++)
			for (int k = 0; k < inst->dimension_f; k++)
				if (j != k && inst->c_f[j][k] <= inst->r)
				{
					inst->obj[counter] = 0.0;
					inst->lb[counter] = 0.0;

					inst->ub[counter] = 1.0;

					inst->c_type[counter] = 'C';
					sprintf(inst->colname[counter], "f_m.%d.%d.%d", i, j, k);
					inst->fm[i][j][k] = counter;
					counter++;
				}

	inst->ccnt = counter;
#ifdef use_lp_solver
	inst->status = CPXnewcols(inst->env_CPIF, inst->lp_CPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, NULL, inst->colname);
#else
	inst->status = CPXnewcols(inst->env_CPIF, inst->lp_CPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
#endif

	if (inst->status != 0)
	{
		cout << "error in CPXnewcols" << endl;
		exit(-1);
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

	for (int i = 0; i < inst->ccnt; i++) { free(inst->colname[i]); }
	free(inst->colname);

	//optimization sense - minimization
	CPXchgobjsen(inst->env_CPIF, inst->lp_CPIF, CPX_MIN);

	//adding offset
	CPXchgobjoffset(inst->env_CPIF, inst->lp_CPIF, inst->total_demand);

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

		inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout << "error in CPXaddrows" << endl;
			exit(-1);
		}

		free(inst->rmatval);
		free(inst->rmatind);
	}

	// z_k <= sum a_ki y_i
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;

	vector<int> rmatind;
	vector<double> rmatval;

	for (int k = 0; k < inst->dimension_c; k++) {
		counter = 0;
		for (int i : inst->fstar_cover[k]) {
			rmatind.push_back(inst->y[i]);
			rmatval.push_back(-1.0);
			counter++;
		}

		rmatind.push_back(inst->z[k]);
		rmatval.push_back(1.0);
		inst->nzcnt = counter + 1;
		inst->status = CPXaddrows(inst->env_CPIF, inst->lp_CPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL);
		if (inst->status != 0)
		{
			cout << "error in CPXaddrows" << endl;
			exit(-1);
		}
		rmatind.clear();
		rmatval.clear();
	}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	add_multiflow_constraints_CPIF(inst);

#ifdef write_prob

	inst->status = CPXwriteprob(inst->env_CPIF, inst->lp_CPIF, "CPIF.lp", NULL);
	if (inst->status != 0) {
		cout << "error in CPXwriteprob" << endl;
		exit(-1);
	}

#endif
}

//solve model by standard solver or automatic Benders
void solve_model_CPIF(instance* inst) {

	CPXsetintparam(inst->env_CPIF, CPX_PARAM_SCRIND, CPX_ON);

	CPXsetintparam(inst->env_CPIF, CPX_PARAM_THREADS, 1);
	//CPXsetintparam(inst->env_CPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);

	inst->status = CPXsetdblparam(inst->env_CPIF, CPX_PARAM_TILIM, inst->timelimit);
	if (inst->status)
	{
		cout << "error in CPX_PARAM_EPRHS" << endl;
	}

	clock_t time_start = clock();

	/*if (strncmp(inst->options, "MZ", 2) == 0) {
		CPXsetintparam(inst->env_CPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
	}*/
	if (strncmp(inst->benders, "autobenders", 4) == 0) {
		CPXsetintparam(inst->env_CPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
		inst->status = CPXsetintparam(inst->env_CPIF, CPXPARAM_Benders_Strategy, CPX_BENDERSSTRATEGY_FULL);


		inst->status = CPXwritebendersannotation(inst->env_CPIF, inst->lp_CPIF, "benders.ann");
		if (inst->status) {
			cout << "error" << endl;
		}

		inst->status = CPXbendersopt(inst->env_CPIF, inst->lp_CPIF);
		if (inst->status) {
			cout << "error" << endl;
			exit(-1);
		}
	}
	else
		if (strncmp(inst->benders, "nobenders", 4) == 0) {

			if (inst->status)
			{
				cout << "error in CPXsetlazyconstraintcallbackfunc" << endl;
			}

#ifdef use_lp_solver
			inst->status = CPXlpopt(inst->env_CPIF, inst->lp_CPIF);
			if (inst->status != 0)
			{
				cout << "error in CPXlpopt" << endl;
				exit(-1);
			}
#else
			//dummy callback
			inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_CPIF, mydummycallback_CPIF, inst);
			inst->status = CPXmipopt(inst->env_CPIF, inst->lp_CPIF);
			if (inst->status != 0)
			{
				cout << "error in CPXmipopt" << endl;
				exit(-1);
			}
#endif
		}


	clock_t time_end = clock();



	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

#ifdef use_lp_solver
	inst->status = CPXgetobjval(inst->env_CPIF, inst->lp_CPIF, &(inst->objval));
#else
	inst->status = CPXgetmipobjval(inst->env_CPIF, inst->lp_CPIF, &(inst->objval));
#endif
	if (inst->status != 0)
	{
		cout << "error in CPXgetmipobjval" << endl;
	}

	cout << "CPLEX objval: " << inst->objval << endl;

	cout << "Time elapsed: " << solution_time;

	print_solution_CPIF(inst);

	write_log_CPIF(inst, solution_time);
}

void solve_model_CPIF_NZ(instance* inst)
{
	CPXsetintparam(inst->env_CPIF, CPX_PARAM_SCRIND, CPX_ON);
	CPXsetintparam(inst->env_CPIF, CPX_PARAM_MIPCBREDLP, CPX_OFF);      
	CPXsetintparam(inst->env_CPIF, CPX_PARAM_PRELINEAR, CPX_OFF);            
	CPXsetintparam(inst->env_CPIF, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
	//CPXsetintparam(inst->env_CPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);

	CPXsetintparam(inst->env_CPIF, CPX_PARAM_THREADS, 1);

	inst->status = CPXsetdblparam(inst->env_CPIF, CPX_PARAM_TILIM, inst->timelimit);
	//inst->status = CPXsetdblparam(inst->env_CPIF, CPX_PARAM_EPAGAP, 1);


	inst->lazy_cuts_connectivity = 0;
	inst->lazy_cuts_coverage = 0;
	inst->user_cuts_coverage = 0;

	double solution_time;

	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_CPIF, mylazycallback_NZ_min, inst);
	if (inst->status)
	{
		cout << "error in CPXsetlazyconstraintcallbackfunc" << endl;
	}



	clock_t time_start = clock();

	double root_solution = calculate_only_root_case_CPIF(inst);
	double two_facility_solution = calculate_2_facility_sol_CPIF(inst);

	inst->status = CPXmipopt(inst->env_CPIF, inst->lp_CPIF);
	if (inst->status != 0)
	{
		cout << "error in CPXmipopt" << endl;
		//exit(-1);
	}

	clock_t time_end = clock();

	solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	inst->status = CPXgetmipobjval(inst->env_CPIF, inst->lp_CPIF, &(inst->objval));
	if (inst->status != 0)
	{
		cout << "error in CPXgetmipobjval" << endl;
	}

	int lpstat = CPXgetstat(inst->env_CPIF, inst->lp_CPIF);

	if (((lpstat == CPXMIP_OPTIMAL) ||
		(lpstat == CPXMIP_OPTIMAL_INFEAS) ||
		//( lpstat ==  CPXMIP_OPTIMAL_RELAXED ) ||
		(lpstat == CPXMIP_OPTIMAL_TOL) ||
		(lpstat == CPXMIP_TIME_LIM_FEAS)) && (inst->objval < min(root_solution, two_facility_solution))) {

		inst->x = (double*)calloc(inst->ccnt, sizeof(double));
		inst->status = CPXgetmipx(inst->env_CPIF, inst->lp_CPIF, inst->x, 0, inst->ccnt - 1);

		inst->opened_facilities = 0;
		for (int i = 0; i < inst->dimension_f; i++) {
			inst->opened_facilities += inst->x[inst->y[i]];
		}
		inst->custom_solution_value = inst->objval;
	}
	else {
		inst->use_custom_solution = true;
		if (root_solution <= two_facility_solution) {
			inst->custom_solution_value = root_solution;
			inst->opened_facilities = 1;
		}
		else {
			inst->custom_solution_value = two_facility_solution;
			inst->opened_facilities = 2;
		}
	}




	write_log_CPIF(inst, solution_time);
	print_solution_CPIF(inst);



	cout << "CPLEX objval: " << inst->objval << endl;

}

void clean_model_CPIF(instance* inst) {
	inst->status = CPXfreeprob(inst->env_CPIF, &(inst->lp_CPIF));
	if (inst->status != 0) { cout << "error in CPXfreeprob" << endl;}

	inst->status = CPXcloseCPLEX(&(inst->env_CPIF));
	if (inst->status != 0) { cout << "error in CPXcloseCPLEX" << endl;}
}
