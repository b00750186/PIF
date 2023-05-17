#include "CPIF.h"
#include "CPIF2.h"

double xSol_CPIF_user_2[MAX_SIZE];
double obj_value_CPIF_user_2;
double rmatval_CPIF_user_2[MAX_SIZE];
int    rmatind_CPIF_user_2[MAX_SIZE];
double rhs_CPIF_user_2;

double rmatval_CPIF_SEC_2[MAX_SIZE];
int    rmatind_CPIF_SEC_2[MAX_SIZE];
double rhs_CPIF_SEC_2;

double xSol_CPIF_2[MAX_SIZE];
double obj_value_CPIF_2;
double rmatval_CPIF_2[MAX_SIZE];
int    rmatind_CPIF_2[MAX_SIZE];
double rhs_CPIF_2;

int CPXPUBLIC myusercallback_NT(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{

	instance* inst = (instance*)cbhandle;
	//cout << "User callback started";
	if (wherefrom != CPX_CALLBACK_MIP_CUT_LAST) return 0;

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

	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_CPIF_user_2, 0, inst->ccnt - 1);
	if (status != 0) {
		cout << "cannot get the x" << endl;
		exit(-1);
	}

	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj_value_CPIF_user_2);
	if (status != 0) {
		cout << "cannot get the x" << endl;
		exit(-1);
	}

	//B1f
	//calculate tilde I_j
	inst->I_tilde = (int*)calloc(inst->dimension, sizeof(int));
	for (int i = 0; i < inst->dimension; i++) {
		if (xSol_CPIF_user_2[i] > 0)
			for (int j = inst->fs_pointer[i][1]; j <= inst->fs_pointer[i][2]; j++) {
				inst->I_tilde[inst->fs[j][1]] += xSol_CPIF_user_2[i];
			}
	}

	for (int i = 0; i < inst->dimension; i++) {
		rmatind_CPIF_user_2[i] = i;
		rmatval_CPIF_user_2[i] = 0;
		for (int j = inst->fs_pointer[i][1]; j <= inst->fs_pointer[i][2]; j++) {
			if (inst->I_tilde[inst->fs[j][1]] < 1) {
				rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
			}

			if (inst->I_tilde[inst->fs[j][1]] == 1 && inst->Js[inst->fs[j][1]] == 1) {
				rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
			}
		}
	}

	rmatind_CPIF_user_2[inst->dimension] = inst->teta; //D == theta
	rmatval_CPIF_user_2[inst->dimension] = -1.0;

	rhs_CPIF_user_2 = 0.0;

	for (int j = 0; j < inst->dimension; j++) {
		if (inst->I_tilde[j] >= 1)
			rhs_CPIF_user_2--; //demand == 1
		if (inst->I_tilde[j] == 1 && inst->Js[j] == 1)
			rhs_CPIF_user_2++; //demand == 1
	}

	//check if cut is violated
	double first_term_and_rhs = 0.0;
	double theta_tilde = 0.0;

	//first term
	for (int i = 0; i < inst->dimension; i++) {
		first_term_and_rhs += rmatval_CPIF_user_2[i] * xSol_CPIF_user_2[rmatind_CPIF_user_2[i]];
	}

	//rhs
	first_term_and_rhs -= rhs_CPIF_user_2;

	//theta tilde
	theta_tilde = xSol_CPIF_user_2[inst->teta];

	//cout << "Theta: " << theta_tilde << endl;
	//cout << "First term and rhs: " << first_term_and_rhs << endl;
	/*int dummy;
	cin>>dummy;*/

	double viol = (theta_tilde - first_term_and_rhs) / theta_tilde * 100;
	viol = max(viol, 0.0001);

	//if ((theta_tilde - first_term_and_rhs) / theta_tilde > 1) {
	//if (theta_tilde - first_term_and_rhs > 0.001) {
	if (viol > 5.0) {
		/*for (int i = 0; i < inst->dimension + 1; i++) {
			cout << rmatval_CPIF_user[i] << " " << inst->colname[rmatind_CPIF_user[i]] <<" ";
		}*/
		//cout << "User callback added";
		inst->user_cuts_coverage++;
		int nzcnt = inst->dimension + 1;
		status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_user_2, 'G', rmatind_CPIF_user_2, rmatval_CPIF_user_2, 0);
		if (status != 0) {
			cout << "CPXcutcallbackadd" << endl;
			exit(-1);
		}


#ifdef write_prob_call_back
		double per_test;
		CPXLPptr lp_test;
		CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
		CPXwriteprob(env, lp_test, "prob.lp", NULL);
		cin.get();
#endif
		(*useraction_p) = CPX_CALLBACK_SET;

	}
	//cout << " finished" << endl;

	free(inst->I_tilde);

	return 0;
}

int CPXPUBLIC myusercallback_NT_2(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p) {

	instance* inst = (instance*)cbhandle;
	//cout << "User callback started";
	if (wherefrom != CPX_CALLBACK_MIP_CUT_LAST) return 0;

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

	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_CPIF_user_2, 0, inst->ccnt - 1);
	if (status != 0) {
		cout << "cannot get the x" << endl;
	}

	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj_value_CPIF_user_2);
	if (status != 0) {
		cout << "cannot get the x" << endl;
	}

	//B2f
	//calculate tilde I_j
	inst->I_tilde = (int*)calloc(inst->dimension, sizeof(int));
	for (int i = 0; i < inst->dimension; i++) {
		if (xSol_CPIF_user_2[i] > 0.0)
			for (auto j:inst->fstar_cover[i]) {
				inst->I_tilde[j] += xSol_CPIF_user_2[i];
			}
	}

	for (int i = 0; i < inst->dimension; i++) {
		rmatind_CPIF_user_2[i] = i;
		rmatval_CPIF_user_2[i] = 0;
		for (int j:inst->fstar_cover[i]) {
			if (inst->I_tilde[j] < 1.0) {
				rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
			}
			if (inst->I_tilde[j] == 1.0) {
				rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
			}
		}
	}

	rmatind_CPIF_user_2[inst->dimension] = inst->teta; //D == theta
	rmatval_CPIF_user_2[inst->dimension] = -1.0;

	rhs_CPIF_user_2 = 0.0;

	for (int j = 0; j < inst->dimension; j++) {
		if (inst->I_tilde[j] >= 1.0)
			rhs_CPIF_user_2--; //demand == 1
		if (inst->I_tilde[j] == 1.0)
			rhs_CPIF_user_2++; //demand == 1
	}

	//check if cut is violated
	double first_term_and_rhs = 0.0;
	double theta_tilde = 0.0;

	//first term
	for (int i = 0; i < inst->dimension; i++) {
		first_term_and_rhs += rmatval_CPIF_user_2[i] * xSol_CPIF_user_2[rmatind_CPIF_user_2[i]];
	}

	//rhs
	first_term_and_rhs -= rhs_CPIF_user_2;

	//theta tilde
	theta_tilde = xSol_CPIF_user_2[inst->teta];

	//cout << "Theta: " << theta_tilde << endl;
	//cout << "First term and rhs: " << first_term_and_rhs << endl;
	/*int dummy;
	cin>>dummy;*/

	double viol = (theta_tilde - (first_term_and_rhs)) / theta_tilde * 100;
	viol = max(viol, 0.0001);

	//if ((theta_tilde - first_term_and_rhs) / theta_tilde > 1) {
	//if (theta_tilde - first_term_and_rhs > 0.001) {
	if (viol > 5.0) {
		/*for (int i = 0; i < inst->dimension + 1; i++) {
			cout << rmatval_CPIF_user[i] << " " << inst->colname[rmatind_CPIF_user[i]] <<" ";
		}*/
		//cout << "User callback added";
		inst->user_cuts_coverage++;
		int nzcnt = inst->dimension + 1;
		status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_user_2, 'G', rmatind_CPIF_user_2, rmatval_CPIF_user_2, 0);
		if (status != 0) {
			cout << "CPXcutcallbackadd" << endl;
			exit(-1);
		}


#ifdef write_prob_call_back
		double per_test;
		CPXLPptr lp_test;
		CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
		CPXwriteprob(env, lp_test, "prob.lp", NULL);
		cin.get();
#endif
		(*useraction_p) = CPX_CALLBACK_SET;

	}
	//cout << " finished" << endl;

	free(inst->I_tilde);

	return 0;
}

int CPXPUBLIC mylazycallback_CPIF_NN(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{

	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle;


	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_CPIF_2, 0, inst->ccnt - 1);
	if (status != 0) {
		cout << "cannot get the x" << endl;
		exit(-1);
	}

	//SEC
	//finding root subtour
	vector <bool> neighbours_to_open(inst->dimension, false);
	vector <bool> visited(inst->dimension, false);
	vector <bool> close_nodes(inst->dimension, false);
	vector <bool> not_reached(inst->dimension, false);
	queue <int> remote_nodes;
	queue <int> remote_nodes_2;
	queue <int> qu;

	//root node subtour
	visited[0] = true;
	close_nodes[0] = true;
	qu.push(0);
	int root_covered = 1;

	while (!qu.empty()) {
		int i = qu.front();
		qu.pop();

		for(auto neighbour:inst->fstar_facility[i]){
			if (!visited[neighbour]) {
				if (xSol_CPIF_2[inst->y[neighbour]] > 0.999) {
					qu.push(neighbour);
					close_nodes[neighbour] = true;
					root_covered++;
				}
				else {
					neighbours_to_open[neighbour] = true;
				}
				visited[neighbour] = true;
			}
		}
	}

	//number of opened facilities
	for (int j = 0; j < inst->dimension; j++) {
		if (xSol_CPIF_2[inst->y[j]] > 0.999 && !close_nodes[j]) {
			remote_nodes.push(j);
			remote_nodes_2.push(j);
		}
	}

	if (!remote_nodes.empty()) {
		int counter = 0;
		for (int j = 0; j < inst->dimension; j++) {
			if (neighbours_to_open[j]) {
				rmatind_CPIF_SEC_2[counter] = inst->y[j];
				rmatval_CPIF_SEC_2[counter] = 1.0;
				counter++;
			}
		}
		int nzcnt = counter + 1;
		rhs_CPIF_SEC_2 = 0;

		while (!remote_nodes.empty()) {
			rmatind_CPIF_SEC_2[counter] = inst->y[remote_nodes.front()];
			not_reached[remote_nodes.front()] = true;
			remote_nodes.pop();
			rmatval_CPIF_SEC_2[counter] = -1.0;

			status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_SEC_2, 'G', rmatind_CPIF_SEC_2, rmatval_CPIF_SEC_2, 0);
			if (status != 0) {
				cout << "CPXcutcallbackadd" << endl;
				exit(-1);
			}
		}

		//cut from remote connected components
		//reset flags
		vector<bool> visited_opened(inst->dimension, false);
		queue<int> opened_facility_in_connected_component;
		vector<bool> neighbours_of_conn_comp(inst->dimension, false);
		vector<bool> visited_in_comp(inst->dimension, false);

		while (!remote_nodes_2.empty()) {
			int k = remote_nodes_2.front();
			remote_nodes_2.pop();

			if (!visited_opened[k]) {

				fill(neighbours_of_conn_comp.begin(), neighbours_of_conn_comp.end(), false);
				fill(visited_in_comp.begin(), visited_in_comp.end(), false);

				qu.push(k);
				visited_opened[k] = true;
				opened_facility_in_connected_component.push(k);

				while (!qu.empty()) {
					int i = qu.front();
					qu.pop();

					for(auto neighbour: inst->fstar_facility[i]){
						if (!visited_in_comp[neighbour]) {
							if (xSol_CPIF_2[inst->y[neighbour]] > 0.999) {
								qu.push(neighbour);

								//queue for cut generation
								opened_facility_in_connected_component.push(neighbour);

								//flag for global cycle for all opened facilities not connected to the root
								visited_opened[neighbour] = true;
							}
							else {
								neighbours_of_conn_comp[neighbour] = true;
							}
							visited_in_comp[neighbour] = true;
						}
					}
				}

				//cut for each opened facility in the connected component
				counter = 0;
				for (int j = 0; j < inst->dimension; j++) {
					if (neighbours_of_conn_comp[j]) {
						rmatind_CPIF_SEC_2[counter] = inst->y[j];
						rmatval_CPIF_SEC_2[counter] = 1.0;
						counter++;
					}
				}
				int nzcnt = counter + 1;
				rhs_CPIF_SEC_2 = 0;

				while (!opened_facility_in_connected_component.empty()) {
					rmatind_CPIF_SEC_2[counter] = inst->y[opened_facility_in_connected_component.front()];
					opened_facility_in_connected_component.pop();
					rmatval_CPIF_SEC_2[counter] = -1.0;

					status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_SEC_2, 'G', rmatind_CPIF_SEC_2, rmatval_CPIF_SEC_2, 0);
					if (status != 0) {
						cout << "CPXcutcallbackadd" << endl;
						exit(-1);
					}
				}

			}
		}
		inst->lazy_cuts_connectivity++;
		(*useraction_p) = CPX_CALLBACK_SET;
	}
	else {
	//calculate z_tilde, here z_tilde[i] == 1 iff customer i is covered given the current solution
	inst->z_tilde = (int*)calloc(inst->dimension, sizeof(int));
	for (int i = 0; i < inst->dimension; i++) {
		if (xSol_CPIF_2[inst->y[i]] > 0.999) {
			for (auto j:inst->fstar_cover[i]){
				inst->z_tilde[j] = 1;
			}
		}
	}

	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj_value_CPIF_2);
	if (status != 0) {
		cout << "cannot get the x" << endl;
		exit(-1);
	}

	int nzcnt;
	double Dminus;

	//B1
	nzcnt = inst->dimension + 1;
	for (int i = 0; i < inst->dimension; i++)
	{
		rmatind_CPIF_2[i] = i;
		rmatval_CPIF_2[i] = 0.0;
		if (xSol_CPIF_2[i] >= 0.9999) {
			for (int j = 0; j < inst->dimension; j++) {
				if (inst->c[i][j] <= inst->R && inst->Js[j] == 1)
					rmatval_CPIF_2[i] += 1; //demand[j] == 1
			}
		}
		else {
			for (int j = 0; j < inst->dimension; j++) {
				if (inst->c[i][j] <= inst->R && inst->z_tilde[j] == 0)
					rmatval_CPIF_2[i] += 1; //demand[j] == 1
			}
		}

	}

	rmatind_CPIF_2[inst->dimension] = inst->teta; //theta variable
	rmatval_CPIF_2[inst->dimension] = -1.0;

	//D(J(K_tilda\Js)
	Dminus = 0;
	for (int i = 0; i < inst->dimension; i++)
		if (inst->Js[i] != 1 && inst->z_tilde[i] > 0.999)
			Dminus--; //demand == 1

	rhs_CPIF_2 = Dminus;

	//check violation
	double lhs = 0.0;
	double rhs = 0.0;

	for (int i = 0; i < nzcnt; i++) {
		lhs += rmatval_CPIF_2[i] * xSol_CPIF_2[rmatind_CPIF_2[i]];
	}

	//cout << "theta: " << xSol_CPIF_2[inst->teta] <<endl;

	rhs = rhs_CPIF_2;

	if (lhs < rhs - 0.00001) {
		status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_2, 'G', rmatind_CPIF_2, rmatval_CPIF_2, 0);
		if (status != 0) {
			cout << "CPXcutcallbackadd" << endl;
			exit(-1);
		}

		inst->lazy_cuts_coverage++;
		(*useraction_p) = CPX_CALLBACK_SET;
	}
	free(inst->z_tilde);
}

#ifdef write_prob_call_back
	double per_test;
	CPXLPptr lp_test;
	CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
	CPXwriteprob(env, lp_test, "prob_NN.lp", NULL);
	cin.get();
#endif
	//cout << " finished ";
	//cout << "Objective: " << obj_value_CPIF << endl;


	return 0;

}

int CPXPUBLIC mylazycallback_CPIF_NN_min(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{

	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle;


	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_CPIF_2, 0, inst->ccnt - 1);
	if (status != 0) {
		cout << "cannot get the x" << endl;
		exit(-1);
	}

	bool network_is_feasible = true;

	//SEC
	//finding root subtour
	vector <bool> neighbours_to_open(inst->dimension, false);
	vector <bool> visited(inst->dimension, false);
	vector <bool> root_component(inst->dimension, false);
	vector <bool> not_reached(inst->dimension, false); //for plotting purposes
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
				if (xSol_CPIF_2[inst->y[neighbour]] > 0.999) {
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

	for (int j = 1; j < inst->dimension; j++) {
		if (xSol_CPIF_2[inst->y[j]] > 0.999 && C[0].count(j) <= 0) {

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
							rmatind_CPIF_2[counter] = inst->y[neighbour];
							rmatval_CPIF_2[counter] = 1.0;
							//cout << "+y" << neighbour;
							counter++;
						}

						visited[neighbour] = true;
					}
				}
			}

			rmatind_CPIF_2[counter] = inst->y[j];
			rmatval_CPIF_2[counter] = -1.0;
			//cout << "-y" << j;

			rhs_CPIF_2 = 0.0;
			int nzcnt = counter + 1;



			//check if cut violated
			double lhs_check = 0, rhs_check = 0;

			for (int i = 0; i <= counter; i++) {
				lhs_check += xSol_CPIF_2[rmatind_CPIF_2[i]] * rmatval_CPIF_2[i];
			}

			//cout << "lhs - rhs = " << lhs_check - rhs_check;
			if (lhs_check - rhs_check != 0.0) {
				//if(true){
				status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_2, 'G', rmatind_CPIF_2, rmatval_CPIF_2, 0);
				inst->lazy_cuts_connectivity++;
				if (status != 0) {
					cout << "CPXcutcallbackadd" << endl;
					exit(-1);
				}

			}


		}
	}

	if(network_is_feasible) {
		//B1f
	//calculate tilde I_j
		inst->I_tilde = (int*)calloc(inst->dimension, sizeof(int));
		for (int i = 0; i < inst->dimension; i++) {
			if (xSol_CPIF_2[i] > 0.9999)
				for (int j = inst->fs_pointer[i][1]; j <= inst->fs_pointer[i][2]; j++) {
					inst->I_tilde[inst->fs[j][1]] += 1.0;
				}
		}

		for (int i = 0; i < inst->dimension; i++) {
			rmatind_CPIF_user_2[i] = i;
			rmatval_CPIF_user_2[i] = 0;
			for (int j = inst->fs_pointer[i][1]; j <= inst->fs_pointer[i][2]; j++) {
				if (inst->I_tilde[inst->fs[j][1]] < 1) {
					rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
				}

				if (inst->I_tilde[inst->fs[j][1]] == 1 && inst->Js[inst->fs[j][1]] == 1) {
					rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
				}
			}
		}

		rmatind_CPIF_user_2[inst->dimension] = inst->teta; //D == theta
		rmatval_CPIF_user_2[inst->dimension] = -1.0;

		rhs_CPIF_user_2 = 0.0;

		for (int j = 0; j < inst->dimension; j++) {
			if (inst->I_tilde[j] >= 1)
				rhs_CPIF_user_2--; //demand == 1
			if (inst->I_tilde[j] == 1 && inst->Js[j] == 1)
				rhs_CPIF_user_2++; //demand == 1
		}

		//check if cut is violated
		double first_term_and_rhs = 0.0;
		double theta_tilde = 0.0;

		//first term
		for (int i = 0; i < inst->dimension; i++) {
			first_term_and_rhs += rmatval_CPIF_user_2[i] * xSol_CPIF_user_2[rmatind_CPIF_user_2[i]];
		}

		//rhs
		first_term_and_rhs -= rhs_CPIF_user_2;

		//theta tilde
		theta_tilde = xSol_CPIF_2[inst->teta];

		//cout << "Theta: " << theta_tilde << endl;
		//cout << "First term and rhs: " << first_term_and_rhs << endl;
		/*int dummy;
		cin>>dummy;*/

		double viol = (theta_tilde - (first_term_and_rhs)) / theta_tilde * 100;
		viol = max(viol, 0.0001);

		//if ((theta_tilde - first_term_and_rhs) / theta_tilde > 1) {
		//if (theta_tilde - first_term_and_rhs > 0.001) {
		if (true) {
			/*for (int i = 0; i < inst->dimension + 1; i++) {
				cout << rmatval_CPIF_user[i] << " " << inst->colname[rmatind_CPIF_user[i]] <<" ";
			}*/
			//cout << "User callback added";
			inst->user_cuts_coverage++;
			int nzcnt = inst->dimension + 1;
			status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_user_2, 'G', rmatind_CPIF_user_2, rmatval_CPIF_user_2, 0);
			if (status != 0) {
				cout << "CPXcutcallbackadd" << endl;
				exit(-1);
			}

		}

		free(inst->I_tilde);
	}

#ifdef write_prob_call_back
	double per_test;
	CPXLPptr lp_test;
	CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
	CPXwriteprob(env, lp_test, "prob_NN.lp", NULL);
	cout << "NODE LP WRITTEN\n\n\n";
	cin.get();
#endif
	//cout << " finished ";
	//cout << "Objective: " << obj_value_CPIF << endl;


	return 0;

}

int CPXPUBLIC mylazycallback_CPIF_NN_min_2(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{

	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle;


	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol_CPIF_2, 0, inst->ccnt - 1);
	if (status != 0) {
		cout << "cannot get the x" << endl;
		exit(-1);
	}

	bool network_is_feasible = true;

	//SEC
	//finding root subtour
	vector <bool> neighbours_to_open(inst->dimension, false);
	vector <bool> visited(inst->dimension, false);
	vector <bool> root_component(inst->dimension, false);
	vector <bool> not_reached(inst->dimension, false); //for plotting purposes
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
				if (xSol_CPIF_2[inst->y[neighbour]] > 0.999) {
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

	for (int j = 1; j < inst->dimension; j++) {
		if (xSol_CPIF_2[inst->y[j]] > 0.999 && C[0].count(j) <= 0) {

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
							rmatind_CPIF_2[counter] = inst->y[neighbour];
							rmatval_CPIF_2[counter] = 1.0;
							//cout << "+y" << neighbour;
							counter++;
						}

						visited[neighbour] = true;
					}
				}
			}

			rmatind_CPIF_2[counter] = inst->y[j];
			rmatval_CPIF_2[counter] = -1.0;
			//cout << "-y" << j;

			rhs_CPIF_2 = 0.0;
			int nzcnt = counter + 1;



			//check if cut violated
			double lhs_check = 0, rhs_check = 0;

			for (int i = 0; i <= counter; i++) {
				lhs_check += xSol_CPIF_2[rmatind_CPIF_2[i]] * rmatval_CPIF_2[i];
			}

			//cout << "lhs - rhs = " << lhs_check - rhs_check;
			if (lhs_check - rhs_check != 0.0) {
				//if(true){
				status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_2, 'G', rmatind_CPIF_2, rmatval_CPIF_2, 0);
				inst->lazy_cuts_connectivity++;
				if (status != 0) {
					cout << "CPXcutcallbackadd" << endl;
					exit(-1);
				}

			}


		}
	}

	if (network_is_feasible) {
		
	//B2f
	//calculate tilde I_j
//calculate tilde I_j
		inst->I_tilde = (int*)calloc(inst->dimension, sizeof(int));
		for (int i = 0; i < inst->dimension; i++) {
			if (xSol_CPIF_2[i] > 0.99999)
				for (auto j : inst->fstar_cover[i]) {
					inst->I_tilde[j] += 1.0;
				}
		}

		for (int i = 0; i < inst->dimension; i++) {
			rmatind_CPIF_user_2[i] = i;
			rmatval_CPIF_user_2[i] = 0.0;
			for (int j : inst->fstar_cover[i]) {
				if (inst->I_tilde[j] < 0.999999) {
					rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
				}
				if (inst->I_tilde[j] >= 0.9999 && inst->Js[j] == 1) {
					rmatval_CPIF_user_2[i] += 1; //demand[j] == 1
				}
			}
		}

		rmatind_CPIF_user_2[inst->dimension] = inst->teta; //D == theta
		rmatval_CPIF_user_2[inst->dimension] = -1.0;

		rhs_CPIF_user_2 = 0.0;

		for (int j = 0; j < inst->dimension; j++) {
			if (inst->I_tilde[j] >= 0.99999)
				rhs_CPIF_user_2--; //demand == 1
			if (inst->Js[j] == 1 && inst->I_tilde[j] >= 0.99999)
				rhs_CPIF_user_2++; //demand == 1
		}

		//check if cut is violated
		double first_term_and_rhs = 0.0;
		double theta_tilde = 0.0;

		//first term
		for (int i = 0; i < inst->dimension; i++) {
			first_term_and_rhs += rmatval_CPIF_user_2[i] * xSol_CPIF_2[rmatind_CPIF_user_2[i]];
		}

		//rhs
		first_term_and_rhs -= rhs_CPIF_user_2;

		//theta tilde
		theta_tilde = xSol_CPIF_2[inst->teta];

		//cout << "Theta: " << theta_tilde << endl;
		//cout << "First term and rhs: " << first_term_and_rhs << endl;
		/*int dummy;
		cin>>dummy;*/

		double viol = (theta_tilde - (first_term_and_rhs)) / theta_tilde * 100;
		viol = max(viol, 0.0001);

		//if ((theta_tilde - first_term_and_rhs) / theta_tilde > 1) {
		//if (theta_tilde - first_term_and_rhs > 0.001) {
		if (true) {
			/*for (int i = 0; i < inst->dimension + 1; i++) {
				cout << rmatval_CPIF_user[i] << " " << inst->colname[rmatind_CPIF_user[i]] <<" ";
			}*/
			//cout << "User callback added";
			inst->lazy_cuts_coverage++;
			int nzcnt = inst->dimension + 1;
			status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs_CPIF_user_2, 'G', rmatind_CPIF_user_2, rmatval_CPIF_user_2, 0);
			if (status != 0) {
				cout << "CPXcutcallbackadd" << endl;
			}
		
			(*useraction_p) = CPX_CALLBACK_SET;
		}
		
		free(inst->I_tilde);
		
	}

#ifdef write_prob_call_back
	double per_test;
	CPXLPptr lp_test;
	CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
	CPXwriteprob(env, lp_test, "prob_NN.lp", NULL);
	cout << "Node LP is written";
	cin.get();
#endif
	//cout << " finished ";
	//cout << "Objective: " << obj_value_CPIF << endl;

	return 0;

}

void add_cardinal_constraint_CPIF(instance* inst) {

	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 3.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'G';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));


	for (int j = 0; j < inst->dimension; j++) {
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

double calculate_only_root_case_CPIF(instance* inst)
{

	double cost = 0.0;
	double teta= inst->fstar_cover[0].size();

	cost = inst->dimension - teta;

	return double(cost);
}

double calculate_2_facility_sol_CPIF(instance* inst) {

	set<int> covered_customers_from_root;
	set<int> covered_from_second_facility_and_root;

	int best_coverage = 0;

	for (auto i : inst->fstar_cover[0])
		covered_customers_from_root.insert(i);

	for (auto i : inst->fstar_facility[0]) {

		covered_from_second_facility_and_root = covered_customers_from_root;
		for (auto j : inst->fstar_cover[i])
			covered_from_second_facility_and_root.insert(j);

		if (best_coverage < covered_from_second_facility_and_root.size())
			best_coverage = covered_from_second_facility_and_root.size();
	}

	return double(inst->dimension - best_coverage + inst->fixed_cost * inst->alpha);
}

void build_model_CPIF_NT(instance* inst) {
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
	inst->ccnt = inst->dimension + 1;
	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(20, sizeof(char)); }

	inst->y = (int*)calloc(inst->dimension, sizeof(int));

	int counter = 0;
	for (int j = 0; j < inst->dimension; j++) {
		inst->obj[counter] = inst->fixed_cost * inst->alpha;
		//cout << "Objectiv coeff y " << j << " : " << inst->obj[counter] << endl;
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

	//master problem teta variable
	inst->obj[counter] = -1.0;
	inst->lb[counter] = 0.0;
	inst->ub[counter] = inst->dimension;
	inst->c_type[counter] = 'C';
	//inst->c_type[counter] = 'I';
	sprintf(inst->colname[counter], "teta");
	inst->teta = counter;
	counter++;

	//finding root subtour
	vector <bool> reachable(inst->dimension, false);
	vector <bool> visited(inst->dimension, false);
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

	for (int k = 0; k < inst->dimension; k++) {
		if (!reachable[k]) {
			inst->ub[inst->y[k]] = 0.0;
		}
	}

	inst->status = CPXnewcols(inst->env_CPIF, inst->lp_CPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
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
	CPXchgobjoffset(inst->env_CPIF, inst->lp_CPIF, inst->dimension);

	// p median
#ifndef ignore_p
	if (strncmp(inst->instance_type, "pmed", 4) == 0) {
		add_pmed_constraint_CPIF(inst);
	}
	else {
		add_cardinal_constraint_CPIF(inst);
	}
#endif // !ignore_p



#ifdef write_prob
	inst->status = CPXwriteprob(inst->env_CPIF, inst->lp_CPIF, "CPIF.lp", NULL);
	if (inst->status != 0) {
		cout << "error in CPXwriteprob" << endl;
		exit(-1);
	}
#endif
}

void solve_model_CPIF_NT(instance* inst) {

	/*
	double root_solution = calculate_only_root_case_CPIF(inst);
	double two_facility_solution = calculate_2_facility_sol(inst);

	FILE* pFile;
	int n;
	char name[100];
	pFile = fopen("log.txt", "a");


	if (root_solution <= two_facility_solution) {
		fprintf(pFile, "obj: \t %f \t fac: \t 1 \n", root_solution);
		cout << 1 <<"  obj = " << root_solution << endl;
	}
	else {
		fprintf(pFile, "obj: \t %f \t fac: \t 2 \n", two_facility_solution);
		cout << 2 << "  obj = " << two_facility_solution << endl;
	}
	*/



	CPXsetintparam(inst->env_CPIF, CPX_PARAM_SCRIND, CPX_ON);
	CPXsetintparam(inst->env_CPIF, CPX_PARAM_MIPCBREDLP, CPX_OFF);        
	CPXsetintparam(inst->env_CPIF, CPX_PARAM_PRELINEAR, CPX_OFF);             
	CPXsetintparam(inst->env_CPIF, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
	//CPXsetintparam(inst->env_CPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);

	CPXsetintparam(inst->env_CPIF, CPX_PARAM_THREADS, 1);

	inst->status = CPXsetdblparam(inst->env_CPIF, CPX_PARAM_TILIM, inst->timelimit);
	if (inst->status)
	{
		cout << "error in CPX_PARAM_EPRHS" << endl;
	}

	inst->lazy_cuts_connectivity = 0;
	inst->lazy_cuts_coverage = 0;
	inst->user_cuts_coverage = 0;

	double solution_time;

	/*
	inst->status = CPXsetusercutcallbackfunc(inst->env_CPIF, myusercallback_NT, inst);
	if (inst->status)
	{
		cout << "error in CPXsetlazyconstraintcallbackfunc" << endl;
	}
	*/
	/*
	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_CPIF, mylazycallback_CPIF_NN_min, inst);
	if (inst->status)
	{
		cout << "error in CPXsetlazyconstraintcallbackfunc" << endl;
	}
	*/
	
	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_CPIF, mylazycallback_CPIF_NN_min_2, inst);
	if (inst->status)
	{
		cout << "error in CPXsetlazyconstraintcallbackfunc" << endl;
	}	
	
	
	inst->status = CPXsetusercutcallbackfunc(inst->env_CPIF, myusercallback_NT_2, inst);
	if (inst->status)
	{
		cout << "error in CPXsetlazyconstraintcallbackfunc" << endl;
	}
	
	
	if (strncmp(inst->instance_type, "pmed", 4) != 0) {

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
			for (int i = 0; i < inst->dimension; i++) {
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

	}
	else {

		clock_t time_start = clock();

		cout << "\nCPXmipopt:\n";

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
			(lpstat == CPXMIP_TIME_LIM_FEAS))) {

			inst->x = (double*)calloc(inst->ccnt, sizeof(double));
			inst->status = CPXgetmipx(inst->env_CPIF, inst->lp_CPIF, inst->x, 0, inst->ccnt - 1);

			inst->opened_facilities = 0.0;
			for (int i = 0; i < inst->dimension; i++) {
				inst->opened_facilities += inst->x[inst->y[i]];
			}
			inst->custom_solution_value = inst->objval;
		}

	}

	write_log_CPIF(inst, solution_time);
	print_solution_CPIF(inst);

	cout << "CPLEX objval: " << inst->objval << endl;

}
