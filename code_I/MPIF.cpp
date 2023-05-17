#include "MPIF.h"
#include "MPIF2.h"

double xSol[MAX_SIZE];
double obj_value;
double rmatval[MAX_SIZE];
int    rmatind[MAX_SIZE];
double rhs;

int CPXPUBLIC mydummycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	return 0;
}

//using only minimal separator
int CPXPUBLIC mylazycallback_NW_min(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p) {
	(*useraction_p) = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle;
	//cout << "Lazy callback" << endl;

	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, xSol, 0, inst->ccnt - 1);
	if (status != 0) {
		cout<<"error in CPXgetcallbacknodex"<<endl;
		
	}
	/*
	for (int i = 0; i < inst->dimension; i++) {

		if (xSol[inst->y[i]] > 0.0) {
			cout << "y" << i << " = " << xSol[inst->y[i]] << "   ";
			cout << endl;
		}

	}
	*/

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
				if (xSol[inst->y[neighbour]] > 0.999) {
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
		if (xSol[inst->y[j]] > 0.999 && C[0].count(j) <= 0) {

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
					if (!visited[neighbour])  {

						if (A[0].count(neighbour) <= 0 ) {
							qu.push(neighbour);
						}
						else {
							//R intersection A_0
							rmatind[counter] = inst->y[neighbour];
							rmatval[counter] = 1.0;
							//cout << "+y" << neighbour;
							counter++;
						}

						visited[neighbour] = true;
					}
				}
			}

			rmatind[counter] = inst->y[j];
			rmatval[counter] = -1.0;
			//cout << "-y" << j;

			rhs = 0.0;
			int nzcnt = counter + 1;



			//check if cut violated
			double lhs_check = 0, rhs_check = 0;

			for (int i = 0; i <= counter; i++) {
				lhs_check += xSol[rmatind[i]] * rmatval[i];
			}

			//cout << "lhs - rhs = " << lhs_check - rhs_check;
			if (lhs_check - rhs_check != 0.0) {
			//if(true){
				status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs, 'G', rmatind, rmatval, 0);
				inst->lazy_cuts_connectivity++;
				if (status != 0) {
					cout<<"CPXcutcallbackadd"<<endl;
					
				}

			}


		}
	}


#ifdef write_prob_call_back
	double per_test;
	CPXLPptr lp_test;
	CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp_test);
	CPXwriteprob(env, lp_test, "prob_NW_min.lp", NULL);
	cout << "Node LP Lazy callback written\n\n\n";
	cin.get();
#endif

	(*useraction_p) = CPX_CALLBACK_SET;
	return 0;
}

void add_singleflow_constraints_MPIF(instance* inst) {
	//\sum_{i \in I \setminus \{0\}} f_{0 i}=\sum_{i \in I} y_{i}-1
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = -1.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension * 2 - 1;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	int counter = 0;
	for (int i = 1; i < inst->dimension; i++) {
		inst->rmatind[counter] = inst->f[0][i];
		inst->rmatval[counter] = 1.0;
		counter++;
	}

	for (int i = 0; i < inst->dimension; i++) {
		inst->rmatind[counter] = inst->y[i];
		inst->rmatval[counter] = -1.0;
		counter++;
	}
	inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
	if (inst->status != 0)
	{
		cout<<"error in CPXaddrows"<<endl;
		
	}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);

	// \sum_{i \in I \setminus \{j\}} f_{i j}-\sum_{i \in I \setminus \{j\}} f_{j i}=y_{j}
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension * 2 - 1;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));


	for (int j = 1; j < inst->dimension; j++) {
		counter = 0;
		for (int i = 0; i < inst->dimension; i++) {
			if (i != j) {
				inst->rmatind[counter] = inst->f[i][j];
				inst->rmatval[counter] = 1.0;
				counter++;

				inst->rmatind[counter] = inst->f[j][i];
				inst->rmatval[counter] = -1.0;
				counter++;
			}
		}


		inst->rmatind[inst->nzcnt - 1] = inst->y[j];
		inst->rmatval[inst->nzcnt - 1] = -1.0;

		inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout<<"error in CPXaddrows"<<endl;
			
		}
	}
	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
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
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	double M;
	if (strncmp(inst->instance_type, "pmed", 4) == 0)
		M = inst->p - 1;
	else
		M = inst->dimension - 1;

	for (int i = 0; i < inst->dimension; i++)
		for (int j = 0; j < inst->dimension; j++) {

			counter = 0;

			inst->rmatind[0] = inst->f[i][j];
			inst->rmatval[0] = 1.0;
			inst->rmatind[1] = inst->f[j][i];
			inst->rmatval[1] = 1.0;
			inst->rmatind[2] = inst->y[i];
			inst->rmatval[2] = -M;
			inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
			if (inst->status != 0)
			{
				cout<<"error in CPXaddrows"<<endl;
				
			}
		}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);
}

void add_multiflow_constraints_MPIF(instance* inst) {
	//at facility k for commmodity i
	cout << "Adding MCF constraints 1"<< endl;

	vector<int> rmatind;
	vector<double> rmatval;
	for (int i = 1; i < inst->dimension; i++)
		for (int k = 0; k < inst->dimension; k++) 
		if(inst->fstar_facility[k].size() > 0) {

			inst->rcnt = 1;
			inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
			inst->rhs[0] = 0.0;
			inst->sense = (char*)calloc(inst->rcnt, sizeof(char));
			inst->sense[0] = 'E';
			inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
			inst->rmatbeg[0] = 0;

			for (int j: inst->fstar_facility[k]) {
				if (k != j && inst->c[k][j] <= inst->r) {
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
			else if (k==0) {
				rmatind.push_back(inst->y[i]);
				rmatval.push_back(1.0);
			}

			inst->nzcnt = rmatind.size();

			inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, rmatind.data(), rmatval.data(), NULL, NULL);
			if (inst->status != 0)
			{
				cout<<"error in CPXaddrows"<<endl;
				
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

	double M;
	if (strncmp(inst->instance_type, "pmed", 4) == 0)
		M = inst->p - 1;
	else
		M = inst->dimension - 1;
	cout << "Adding MCF constraints 2" << endl;
	for (int i = 1; i < inst->dimension; i++)
	for (int k = 0; k < inst->dimension; k++) {
		for (int j:inst->fstar_facility[k]) {
			{
				if (j != k && inst->c[j][k] <= inst->c[j][k]) {
					inst->rmatind[0] = inst->fm[i][j][k];
					inst->rmatval[0] = 1.0;
					inst->rmatind[1] = inst->fm[i][k][j];
					inst->rmatval[1] = 1.0;
					inst->rmatind[2] = inst->y[k];
					inst->rmatval[2] = -1.0;
					inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
					if (inst->status != 0)
					{
						cout<<"error in CPXaddrows"<<endl;
						
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

void add_pmed_constraint(instance* inst) {
	// p median

		inst->rcnt = 1;
		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->rhs[0] = inst->p;
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
		inst->sense[0] = 'E';
		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatbeg[0] = 0;
		inst->nzcnt = inst->dimension;
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));


		for (int j = 0; j < inst->dimension; j++) {
			inst->rmatind[j] = inst->y[j];
			inst->rmatval[j] = 1.0;
		}

		inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout<<"error in CPXaddrows"<<endl;
			
		}

		free(inst->rmatbeg);
		free(inst->rmatval);
		free(inst->rmatind);
		free(inst->rhs);
		free(inst->sense);

}

void initialize_v(instance* inst) {

	inst->b = vector<vector<bool>>(inst->dimension, vector<bool>(inst->dimension, false));

	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'G';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = 2;
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	for (int j = 0; j < inst->dimension; j++) {

		inst->rhs[0] = inst->c[j][inst->index_sorted[j][1]];

		inst->rmatind[0] = inst->y[j];
		inst->rmatval[0] = inst->c[j][inst->index_sorted[j][1]];
		inst->rmatind[1] = inst->v[j];
		inst->rmatval[1] = 1.0;

		inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout<<"error in CPXaddrows"<<endl;
			
		}
	}
	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);
}

//single commmodity with w variables
void build_model_MPIF_SW(instance* inst)
{
	cout << "Single commodity flow formulation" << endl;
	//CPLEX environment
	inst->env_MPIF = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		cout<<"error in CPXopenCPLEX"<<endl;
		
	}

	//CPLEX problem
	inst->lp_MPIF = CPXcreateprob(inst->env_MPIF, &(inst->status), "MPIF");
	if (inst->status != 0)
	{
		cout<<"error in CPXcreateprob"<<endl;
		
	}

	//adding variables
	inst->ccnt = inst->dimension * inst->dimension + inst->dimension + inst->dimension * inst->dimension - inst->dimension;
	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(20, sizeof(char)); }

	inst->y = (int*)calloc(inst->dimension, sizeof(int));
	inst->w = (int**)calloc(inst->dimension, sizeof(int*));
	for (int i = 0; i < inst->dimension; i++) { inst->w[i] = (int*)calloc(inst->dimension, sizeof(int)); }
	inst->f = (int**)calloc(inst->dimension, sizeof(int*));
	for (int i = 0; i < inst->dimension; i++) { inst->f[i] = (int*)calloc(inst->dimension, sizeof(int)); }

	int counter = 0;
	for (int j = 0; j < inst->dimension; j++) {

		inst->obj[counter] = inst->fixed_cost;
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


	for (int i = 0; i < inst->dimension; i++)
		for (int j = 0; j < inst->dimension; j++)
			if(i!=j)
			{
				//inst->obj[counter] = inst->demand[j];
				inst->obj[counter] = 1.0 * inst->c[i][j];
				inst->lb[counter] = 0.0;

				if (inst->c[i][j] > inst->R)
					inst->ub[counter] = 0.0;
				else
					inst->ub[counter] = 1.0;


#ifdef relax_w
				inst->c_type[counter] = 'C';
#else
				inst->c_type[counter] = 'B';
#endif

				sprintf(inst->colname[counter], "w.%d.%d", i, j);
				inst->w[i][j] = counter;
				counter++;
			}

	for (int i = 0; i < inst->dimension; i++)
		for (int j = 0; j < inst->dimension; j++)
		{
			inst->obj[counter] = 0.0;
			inst->lb[counter] = 0.0;

			if (inst->c[i][j] > inst->r || i==j)
				inst->ub[counter] = 0.0;
			else
				inst->ub[counter] = CPX_INFBOUND;

			inst->c_type[counter] = 'C';
			sprintf(inst->colname[counter], "f.%d.%d", i, j);
			inst->f[i][j] = counter;
			counter++;
		}

#ifdef use_lp_solver
	inst->status = CPXnewcols(inst->env_MPIF, inst->lp_MPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, NULL, inst->colname);
#else
	inst->status = CPXnewcols(inst->env_MPIF, inst->lp_MPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
#endif // use_lp_solver

	if (inst->status != 0)
	{
		cout<<"error in CPXnewcols"<<endl;
		
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

   	for (int i = 0; i < inst->ccnt; i++) { free(inst->colname[i]); }
	free(inst->colname);

	//optimization sense - minimization
	CPXchgobjsen(inst->env_MPIF, inst->lp_MPIF, CPX_MIN);

	//adding constraints

	//y_{k}+\sum_{i \in I \setminus {k}} w_{k i} = 1
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 1.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension;

	for (int k = 0; k < inst->dimension; k++) {

		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->dimension; i++)
		{
			if (i != k) {
				inst->rmatval[i] = 1.0;
				inst->rmatind[i] = inst->w[k][i];
			}
		}

		inst->rmatval[inst->dimension-1] = 1.0;
		inst->rmatind[inst->dimension-1] = inst->y[k];

		inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout<<"error in CPXaddrows"<<endl;
			
		}

		free(inst->rmatval);
		free(inst->rmatind);
	}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	//w_{k i} \leq y_{i}
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = 2;

	for (int k = 0; k < inst->dimension; k++)
		for (int i = 0; i < inst->dimension; i++)
			if (k != i) {
				inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
				inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));
				inst->rmatval[0] = 1.0;
				inst->rmatval[1] = -1.0;
				inst->rmatind[0] = inst->w[k][i];
				inst->rmatind[1] = inst->y[i];
				inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
				if (inst->status != 0)
				{
					cout<<"error in CPXaddrows"<<endl;
					
				}

				free(inst->rmatval);
				free(inst->rmatind);
			}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	if (strncmp(inst->instance_type, "pmed", 4) == 0) {
		//add_pmed_constraint(inst);
	}

	add_singleflow_constraints_MPIF(inst);



#ifdef write_prob
	inst->status = CPXwriteprob(inst->env_MPIF, inst->lp_MPIF, "MPIF.lp", NULL);
	if (inst->status != 0) {
		cout<<"error in CPXwriteprob"<<endl;
		
	}
#endif
}

//multi commmodity with w
void build_model_MPIF_MW(instance* inst) {
	cout << "Multi commodity ..." << endl;
	//CPLEX environment
	inst->env_MPIF = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		cout<<"error in CPXopenCPLEX"<<endl;
		
	}

	//CPLEX problem
	inst->lp_MPIF = CPXcreateprob(inst->env_MPIF, &(inst->status), "MPIF");
	if (inst->status != 0)
	{
		cout<<"error in CPXcreateprob"<<endl;
		
	}

	//adding variables
	inst->ccnt = inst->dimension * inst->dimension + inst->dimension + inst->dimension * inst->dimension*inst->dimension - inst->dimension;
	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));
    inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(20, sizeof(char)); }

	inst->y = (int*)calloc(inst->dimension, sizeof(int));
	inst->w = (int**)calloc(inst->dimension, sizeof(int*));
	for (int i = 0; i < inst->dimension; i++) { inst->w[i] = (int*)calloc(inst->dimension, sizeof(int)); }
	inst->fm = (int***)calloc(inst->dimension, sizeof(int**));
	for (int i = 0; i < inst->dimension; i++) { inst->fm[i] = (int**)calloc(inst->dimension, sizeof(int*)); }
	for (int i = 0; i < inst->dimension; i++)
		for (int j = 0; j < inst->dimension; j++)
		{
			inst->fm[i][j] = (int*)calloc(inst->dimension, sizeof(int));
		}



	int counter = 0;
	for (int j = 0; j < inst->dimension; j++) {

		inst->obj[counter] = inst->fixed_cost;
		inst->lb[counter] = 0.0;
		inst->ub[counter] = 1.0;
#ifdef relax_y
		inst->c_type[counter] = 'C';
#else
		inst->c_type[counter] = 'B';
#endif
		sprintf(inst->colname[counter], "y.%d", j);
		inst->y[j] = counter;
		counter++;
	}


	// y_0 = 1
	inst->lb[inst->y[0]] = 1.0;
	inst->obj[inst->y[0]] = 0.0;

	for (int i = 0; i < inst->dimension; i++)
		for (int j = 0; j < inst->dimension; j++)
			if(i!=j)
			{
				//inst->obj[counter] = inst->demand[j];
				inst->obj[counter] = 1.0 * inst->c[i][j];
				inst->lb[counter] = 0.0;

				if (inst->c[i][j] > inst->R)
					inst->ub[counter] = 0.0;
				else
					inst->ub[counter] = 1.0;

				inst->c_type[counter] = 'B';
#ifdef relax_w
				inst->c_type[counter] = 'C';
#endif

				sprintf(inst->colname[counter], "w.%d.%d", i, j);
				inst->w[i][j] = counter;
				counter++;
			}

	for (int i = 0; i < inst->dimension; i++)
		for (int j = 0; j < inst->dimension; j++)
			for (int k = 0; k < inst->dimension; k++)
			if(j!=k && inst->c[j][k] <= inst->r) {
				inst->obj[counter] = 0.0;
				inst->lb[counter] = 0.0;

				if (inst->c[j][k] > inst->r)
					inst->ub[counter] = 0.0;
				else
					inst->ub[counter] = CPX_INFBOUND;

				inst->c_type[counter] = 'C';
				sprintf(inst->colname[counter], "f_m.%d.%d.%d", i, j, k);
				inst->fm[i][j][k] = counter;
				counter++;
			}
	inst->ccnt = counter;
#ifdef use_lp_solver
	inst->status = CPXnewcols(inst->env_MPIF, inst->lp_MPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, NULL, inst->colname);
#else
	inst->status = CPXnewcols(inst->env_MPIF, inst->lp_MPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
#endif // use_lp_solver


	if (inst->status != 0)
	{
		cout<<"error in CPXnewcols"<<endl;
		
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

    for (int i = 0; i < inst->ccnt; i++) { free(inst->colname[i]); }
	free(inst->colname);


	//optimization sense - minimization
	CPXchgobjsen(inst->env_MPIF, inst->lp_MPIF, CPX_MIN);

	//adding constraints
	cout << "Adding constraints: assignment" << endl;
	//y_{k}+\sum_{i \in I \setminus {k}} w_{k i} = 1
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 1.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension;

	for (int k = 0; k < inst->dimension; k++) {

		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->dimension; i++)
		{
			if (i!=k) {
				inst->rmatval[i] = 1.0;
				inst->rmatind[i] = inst->w[k][i];
			}
		}
		inst->rmatval[inst->dimension-1] = 1.0;
		inst->rmatind[inst->dimension-1] = inst->y[k];

		inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout<<"error in CPXaddrows"<<endl;
			
		}
		free(inst->rmatval);
		free(inst->rmatind);
	}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	cout << "Adding constraints: existing facilities" << endl;
	//w_{k i} \leq y_{i}
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = 2;

	for (int k = 0; k < inst->dimension; k++)
		for (int i = 0; i < inst->dimension; i++)
			if (k != i) {
				inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
				inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));
				inst->rmatval[0] = 1.0;
				inst->rmatval[1] = -1.0;
				inst->rmatind[0] = inst->w[k][i];
				inst->rmatind[1] = inst->y[i];
				inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
				if (inst->status != 0)
				{
					cout<<"error in CPXaddrows"<<endl;
					
				}
				free(inst->rmatval);
				free(inst->rmatind);
			}

	free(inst->rmatbeg);

	free(inst->rhs);
	free(inst->sense);

	if (strncmp(inst->instance_type, "pmed", 4) == 0) {
		//add_pmed_constraint(inst);
	}

	add_multiflow_constraints_MPIF(inst);

#ifdef write_prob
	inst->status = CPXwriteprob(inst->env_MPIF, inst->lp_MPIF, "MPIF.lp", NULL);
	if (inst->status != 0) {
		cout<<"error in CPXwriteprob"<<endl;
		
	}
#endif
}

void build_model_MPIF_NW(instance* inst)
{
	cout << "Node separator with W variables" << endl;
	//CPLEX environment
	inst->env_MPIF = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		cout<<"error in CPXopenCPLEX"<<endl;
		
	}

	//CPLEX problem
	inst->lp_MPIF = CPXcreateprob(inst->env_MPIF, &(inst->status), "MPIF");
	if (inst->status != 0)
	{
		cout<<"error in CPXcreateprob"<<endl;
		
	}

	//adding variables
	inst->ccnt = inst->dimension + inst->dimension * inst->dimension - inst->dimension;
	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(20, sizeof(char)); }

	inst->y = (int*)calloc(inst->dimension, sizeof(int));
	inst->w = (int**)calloc(inst->dimension, sizeof(int*));
	for (int i = 0; i < inst->dimension; i++) { inst->w[i] = (int*)calloc(inst->dimension, sizeof(int)); }

	int counter = 0;
	for (int j = 0; j < inst->dimension; j++) {

		inst->obj[counter] = inst->fixed_cost;
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

	for (int i = 0; i < inst->dimension; i++)
		for (int j = 0; j < inst->dimension; j++)
			if (i != j)
			{
				//inst->obj[counter] = inst->demand[j];
				inst->obj[counter] = 1.0 * inst->c[i][j];
				inst->lb[counter] = 0.0;

				if (inst->c[i][j] > inst->R)
					inst->ub[counter] = 0.0;
				else
					inst->ub[counter] = 1.0;


#ifdef relax_w
				inst->c_type[counter] = 'C';
#else
				inst->c_type[counter] = 'B';
#endif

				sprintf(inst->colname[counter], "w.%d.%d", i, j);
				inst->w[i][j] = counter;
				counter++;
			}

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

		for (auto neighbour:inst->fstar_facility[i])
		{
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

#ifdef use_lp_solver
	inst->status = CPXnewcols(inst->env_MPIF, inst->lp_MPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, NULL, inst->colname);
#else
	inst->status = CPXnewcols(inst->env_MPIF, inst->lp_MPIF, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
#endif // use_lp_solver

	if (inst->status != 0)
	{
		cout<<"error in CPXnewcols"<<endl;
		
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

    for (int i = 0; i < inst->ccnt; i++) { free(inst->colname[i]); }
	free(inst->colname);

	//optimization sense - minimization
	CPXchgobjsen(inst->env_MPIF, inst->lp_MPIF, CPX_MIN);

	//adding constraints

	//y_{k}+\sum_{i \in I} w_{k i} = 1
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 1.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'E';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = inst->dimension + 1;

	for (int k = 0; k < inst->dimension; k++) {

		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->dimension; i++)
		{
			if (i != k) {
				inst->rmatval[i] = 1.0;
				inst->rmatind[i] = inst->w[k][i];
			}

		}

		inst->rmatval[inst->dimension] = 1.0;
		inst->rmatind[inst->dimension] = inst->y[k];

		inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
		if (inst->status != 0)
		{
			cout<<"error in CPXaddrows"<<endl;
			
		}

		free(inst->rmatval);
		free(inst->rmatind);
	}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	//w_{k i} \leq y_{i}
	inst->rcnt = 1;
	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->rhs[0] = 0.0;
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));
	inst->sense[0] = 'L';
	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatbeg[0] = 0;
	inst->nzcnt = 2;

	for (int k = 0; k < inst->dimension; k++)
		for (int i = 0; i < inst->dimension; i++)
			if (k != i)
			{
				inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
				inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));
				inst->rmatval[0] = 1.0;
				inst->rmatval[1] = -1.0;
				inst->rmatind[0] = inst->w[k][i];
				inst->rmatind[1] = inst->y[i];
				inst->status = CPXaddrows(inst->env_MPIF, inst->lp_MPIF, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
				if (inst->status != 0)
				{
					cout<<"error in CPXaddrows"<<endl;
					
				}

				free(inst->rmatval);
				free(inst->rmatind);
			}

	free(inst->rmatbeg);
	free(inst->rhs);
	free(inst->sense);

	if (strncmp(inst->instance_type, "pmed", 4) == 0) {
		add_pmed_constraint(inst);
	}
	else {
		add_cardinal_constraint(inst);
	}


#ifdef write_prob
	inst->status = CPXwriteprob(inst->env_MPIF, inst->lp_MPIF, "MPIF.lp", NULL);
	if (inst->status != 0) {
		cout<<"error in CPXwriteprob"<<endl;
		
	}

#endif
}

//SOLVE

//solve by standard solver or automatic Benders
void solve_model_MPIF(instance* inst) {

	CPXsetintparam(inst->env_MPIF, CPX_PARAM_SCRIND, CPX_ON);
	//CPXsetintparam(inst->env_MPIF, CPXPARAM_Preprocessing_NumPass, 1);

	CPXsetintparam(inst->env_MPIF, CPX_PARAM_THREADS, 1);

	inst->status = CPXsetdblparam(inst->env_MPIF, CPX_PARAM_TILIM, inst->timelimit);
	if (inst->status)
	{
		cout<<"error setting CPX_PARAM_EPRHS"<<endl;
	}

	clock_t time_start = clock();

	//if (strncmp(inst->options, "MW", 2) == 0) {
	//	CPXsetintparam(inst->env_MPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
	//}

	if (strncmp(inst->benders, "autobenders", 4) == 0) {
		CPXsetintparam(inst->env_MPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
		inst->status = CPXsetintparam(inst->env_MPIF, CPXPARAM_Benders_Strategy, CPX_BENDERSSTRATEGY_FULL);

		inst->status = CPXwritebendersannotation(inst->env_MPIF, inst->lp_MPIF, "benders.ann");
		if (inst->status) {
			cout<<"Failed to write the annotation file."<<endl;
		}

		inst->status = CPXbendersopt(inst->env_MPIF, inst->lp_MPIF);
		if (inst->status) {
			cout<<"Failure in optimization."<<endl;
			
		}

		/*inst->status = CPXmipopt(inst->env_MPIF, inst->lp_MPIF);
		if (inst->status != 0)
		{
			cout<<"error in CPXmipopt"<<endl;
			
		}*/
	}
	else
	if (strncmp(inst->benders, "nobenders", 4) == 0){
		
#ifdef use_lp_solver
		CPXsetintparam(inst->env_MPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
		//CPXsetintparam(inst->env_MPIF, CPX_PARAM_PREDUAL, CPX_OFF);
		//CPXsetintparam(inst->env_MPIF, CPX_PARAM_PREIND, CPX_OFF);
		//CPXsetintparam(inst->env_MPIF, CPX_PARAM_RELAXPREIND, CPX_OFF);
		//CPXsetintparam(inst->env_MPIF, CPX_PARAM_PREPASS, CPX_OFF);


		inst->status = CPXlpopt(inst->env_MPIF, inst->lp_MPIF);
		if (inst->status != 0)
		{
			cout<<"error in CPXlpopt"<<endl;
			
	}
#else
		//dummy callback
		inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_MPIF, mydummycallback, inst);
		if (inst->status)
		{
			cout<<"error in CPXsetlazyconstraintcallbackfunc"<<endl;
		}

		inst->status = CPXmipopt(inst->env_MPIF, inst->lp_MPIF);
		if (inst->status != 0)
		{
			cout<<"error in CPXmipopt"<<endl;
			
		}
#endif


	}



	clock_t time_end = clock();


	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;
#ifdef use_lp_solver
	inst->status = CPXgetobjval(inst->env_MPIF, inst->lp_MPIF, &(inst->objval));
#else
	inst->status = CPXgetmipobjval(inst->env_MPIF, inst->lp_MPIF, &(inst->objval));
#endif
	if (inst->status != 0)
	{
		cout<<"error in CPXgetmipobjval"<<endl;
	}

	cout<<"CPLEX objval: " << inst->objval <<endl;

	cout << "Time elapsed: " <<solution_time;

	print_solution_MPIF(inst);

	write_log_MPIF(inst,solution_time);
}

void solve_model_MPIF_NW(instance* inst) {

	CPXsetintparam(inst->env_MPIF, CPX_PARAM_SCRIND, CPX_ON);
	CPXsetintparam(inst->env_MPIF, CPX_PARAM_THREADS, 1);
	inst->status = CPXsetdblparam(inst->env_MPIF, CPX_PARAM_TILIM, inst->timelimit);
	if (inst->status)
	{
		cout<<"error setting CPX_PARAM_EPRHS"<<endl;
	}

	CPXsetintparam(inst->env_MPIF, CPX_PARAM_MIPCBREDLP, CPX_OFF);        
	CPXsetintparam(inst->env_MPIF, CPX_PARAM_PRELINEAR, CPX_OFF);        
	CPXsetintparam(inst->env_MPIF, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);

	//CPXsetintparam(inst->env_MPIF, CPXPARAM_Preprocessing_Presolve, CPX_OFF);

	//with minimal separator only
	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_MPIF, mylazycallback_NW_min, inst);

	if (inst->status)
	{
		cout<<"error in CPXsetlazyconstraintcallbackfunc"<<endl;
	}

	inst->lazy_cuts_connectivity = 0;
	inst->lazy_cuts_cost = 0;
	inst->user_cuts_cost = 0;

	double solution_time;

	if (strncmp(inst->instance_type, "pmed", 4) != 0) {
		clock_t time_start = clock();

		double root_solution = calculate_only_root_case(inst);
		double two_solution = calculate_only_2(inst);

		double cut_up = min(root_solution, two_solution);

		CPXsetdblparam(inst->env_MPIF, CPX_PARAM_CUTUP, cut_up);

		inst->status = CPXmipopt(inst->env_MPIF, inst->lp_MPIF);
		if (inst->status != 0)
		{
			cout<<"error in CPXmipopt"<<endl;
		}

		clock_t time_end = clock();

		solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

		inst->status = CPXgetmipobjval(inst->env_MPIF, inst->lp_MPIF, &(inst->objval));
		if (inst->status != 0)
		{
			cout<<"error in CPXgetmipobjval"<<endl;
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
			for (int i = 0; i < inst->dimension; i++) {
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
	}
	else {
		clock_t time_start = clock();

		inst->status = CPXmipopt(inst->env_MPIF, inst->lp_MPIF);
		if (inst->status != 0)
		{
			cout<<"error in CPXmipopt"<<endl;
		}

		clock_t time_end = clock();

		solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

		inst->status = CPXgetmipobjval(inst->env_MPIF, inst->lp_MPIF, &(inst->objval));
		if (inst->status != 0)
		{
			cout<<"error in CPXgetmipobjval"<<endl;
		}
	}

	print_solution_MPIF(inst);
	write_log_MPIF(inst, solution_time);

	cout << "CPLEX objval: " << inst->objval<<endl;
	cout << "Time elapsed: " << solution_time <<endl;

	cout << "Lazy cuts connectivity: " << inst->lazy_cuts_connectivity << "  Lazy cuts cost: " << inst->lazy_cuts_cost << " User cuts cost: " << inst->user_cuts_cost << endl;

}

void clean_model_MPIF(instance* inst) {
	inst->status = CPXfreeprob(inst->env_MPIF, &(inst->lp_MPIF));
	if (inst->status != 0) { 
		cout<<"error in CPXfreeprob"<<endl;  
	}

	inst->status = CPXcloseCPLEX(&(inst->env_MPIF));
	if (inst->status != 0) { 
		cout<<"error in CPXcloseCPLEX"<<endl;  
	}
}
