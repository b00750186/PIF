#include "global_functions.h"

//#define print_star

//Floyd Warshall algorithm to treat some instances
/*
void floyd(instance* inst) {
	for (int k = 0; k < inst->dimension_f; k++)
	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {
			if (inst->c_f[i][j] > inst->c_f[i][k] + inst->c_f[k][j])
				inst->c_f[i][j] = inst->c_f[i][k] + inst->c_f[k][j];
		}



	//distance matrix
	cout << "Distance matrix" <<endl;
	for (int i = 0; i < inst->dimension; i++) {
		for (int j = 0; j < inst->dimension; j++) {
			cout << inst->c[i][j] << " ";
		}
		cout << endl;
	}

}
*/

double distance(double x1, double x2, double y1, double y2) {

	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

}

void read_input(instance* inst)
/*****************************************************************/
{
	cout << "Reading input file \t" << inst->input_file << " ... " << endl;

	inst->instance_name = (char*)calloc(100, sizeof(char));

	string s = inst->input_file;
	string delimiter = "\\";

	size_t pos = 0;

	while ((pos = s.find(delimiter)) != string::npos) {
		s.erase(0, pos + delimiter.length());
	}

	strcpy(inst->instance_name, s.c_str());

	ifstream in(inst->input_file);
	string line;

	if (!in)
	{
		ofstream err("Error.log", ios::app);
		cout << "Error while opening file. " << endl;
		exit(1);
	}

	in >> inst->dimension_f;
	in >> inst->dimension_c;

	inst->fixed_cost_ind = vector<double>(inst->dimension_f);
	inst->demand = vector<double>(inst->dimension_c);

	inst->pos_f = vector<vector<double>>(inst->dimension_f, vector<double>(2,0.0));
	inst->pos_c = vector<vector<double>>(inst->dimension_c, vector<double>(2,0.0));

	for (int j = 0; j < inst->dimension_f; j++)
	{
		char dummy_char;
		int dummy_int;

		in >> dummy_char;
		in >> dummy_int;

		in >> inst->pos_f[j][0];
		in >> inst->pos_f[j][1];
		in >> inst->fixed_cost_ind[j];
	}

	inst->total_demand = 0.0;

	for (int j = 0; j < inst->dimension_c; j++)
	{

		char dummy_char;
		int dummy_int;

		in >> dummy_char;
		in >> dummy_int;

		in >> inst->pos_c[j][0];
		in >> inst->pos_c[j][1];
		in >> inst->demand[j];
		inst->total_demand += inst->demand[j];

	}

	inst->c_fc = vector<vector<double>>(inst->dimension_f,vector<double>(inst->dimension_c));
	inst->c_f = vector<vector<double>>(inst->dimension_f, vector<double>(inst->dimension_f));
	//forwardstar with vectors

	//average distance
	double sum_dist_f = 0.0;
	double sum_dist_fc = 0.0;

	for (int i = 0; i < inst->dimension_c; i++) {
		inst->fstar_cover.push_back(vector<int>());
		for (int j = 0; j < inst->dimension_f; j++) {
			inst->c_fc[j][i] = distance(inst->pos_f[j][0], inst->pos_c[i][0], inst->pos_f[j][1], inst->pos_c[i][1]);
			sum_dist_fc += inst->c_fc[j][i];
			if (inst->c_fc[j][i] <= inst->R) {
				inst->fstar_cover[i].push_back(j);
			}
		}
	}

	for (int i = 0; i < inst->dimension_f; i++) {
		inst->fstar_facility.push_back(vector<int>());
		for (int j = 0; j < inst->dimension_f; j++) {
			inst->c_f[i][j] = distance(inst->pos_f[i][0],inst->pos_f[j][0],inst->pos_f[i][1],inst->pos_f[j][1]);
			sum_dist_f += inst->c_f[i][j];
			if (inst->c_f[i][j] <= inst->r) {
				inst->fstar_facility[i].push_back(j);
			}
		}
	}

	for (int i = 0; i < inst->dimension_f; i++) {
		inst->fstar_cover_by_facility.push_back(vector<int>());
		for (int j = 0; j < inst->dimension_c; j++) {
			if (inst->c_fc[i][j] <= inst->R) {
				inst->fstar_cover_by_facility[i].push_back(j);
			}
		}
	}


#ifdef print_star

	for (int i = 0; i < inst->dimension_f; i++) {
		cout << i << ":" << endl;
		for (auto j : inst->fstar_facility[i])
			cout << j << " ";
		cout << endl;
	}

#endif  print_star

	//Js clients that can be covered by one facility only
	inst->Js = (int*)calloc(inst->dimension_c, sizeof(int));
	cout << "Js: " <<endl;
	for (int i = 0; i < inst->dimension_c; i++)
		if (inst->fstar_cover[i].size() == 1) {
			inst->Js[i] = 1;
			cout << i << " ";
		}

    if (inst->fixed_cost_ind[0] == 1)
			inst->alpha = 0.0;

	cout << "Average distance inter-facility: "<<sum_dist_f/(inst->dimension_f*(inst->dimension_f-1))<<" average distance facility-customer: "<<sum_dist_fc/(inst->dimension_f*inst->dimension_c)<<endl;

}

void free_and_null(char** ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
}

void print_solution_MPIF(instance* inst) {
	//free(inst->x);
	inst->x = (double*)calloc(inst->ccnt, sizeof(double));
	inst->status = CPXgetmipx(inst->env_MPIF, inst->lp_MPIF, inst->x, 0, inst->ccnt - 1);
	if (inst->status != 0)
	{
		printf("error in CPXgetmipx\n");
	}

#ifdef print_solution

	//getting variables names
	char** cur_colname = NULL;
	char* cur_colnamestore = NULL;
	int           cur_colnamespace;
	int surplus;
	int cur_numrows, cur_numcols;
	cur_numcols = CPXgetnumcols(inst->env_MPIF, inst->lp_MPIF);
	inst->status = CPXgetcolname(inst->env_MPIF, inst->lp_MPIF, NULL, NULL, 0, &surplus, 0, cur_numcols - 1);
	cur_colnamespace = -surplus;
	cur_colname = (char**)malloc(sizeof(char*) * cur_numcols);
	cur_colnamestore = (char*)malloc(cur_colnamespace);

	inst->status = CPXgetcolname(inst->env_MPIF, inst->lp_MPIF, cur_colname, cur_colnamestore, cur_colnamespace, &surplus, 0, cur_numcols - 1);



	printf("\n\nAll variables:\n");
	for (int i = 0; i < cur_numcols; i++) {
		if (inst->x[i] > 0)
			cout << cur_colname[i] << " : " << inst->x[i] << "\t";
	}
#endif

#ifdef drawing
	draw_tex_MPIF(inst, inst->x);
#endif // drawing

}

void print_solution_CPIF(instance* inst) {

	inst->x = (double*)calloc(inst->ccnt, sizeof(double));
	inst->status = CPXgetmipx(inst->env_CPIF, inst->lp_CPIF, inst->x, 0, inst->ccnt - 1);
	if (inst->status != 0)
	{
		printf("error in CPXgetmipx\n");
	}

#ifdef print_solution
	//getting variables names
	char** cur_colname = NULL;
	char* cur_colnamestore = NULL;
	int           cur_colnamespace;
	int surplus;
	int cur_numrows, cur_numcols;
	cur_numcols = CPXgetnumcols(inst->env_CPIF, inst->lp_CPIF);
	inst->status = CPXgetcolname(inst->env_CPIF, inst->lp_CPIF, NULL, NULL, 0, &surplus, 0, cur_numcols - 1);
	cur_colnamespace = -surplus;
	cur_colname = (char**)malloc(sizeof(char*) * cur_numcols);
	cur_colnamestore = (char*)malloc(cur_colnamespace);

	inst->status = CPXgetcolname(inst->env_CPIF, inst->lp_CPIF, cur_colname, cur_colnamestore, cur_colnamespace, &surplus, 0, cur_numcols - 1);

	printf("\n\nAll variables:\n");
	for (int i = 0; i < cur_numcols; i++) {
		if (inst->x[i] > 0)
			cout << cur_colname[i] << " : " << inst->x[i] << "\t";
	}
#endif

#ifdef drawing
	//draw_tex_CPIF(inst, inst->x);
#endif // drawing
}

void calculate_sorted_c(instance* inst) {

	inst->index_sorted = vector<vector<int>>(inst->dimension_c, vector<int>(inst->dimension_f));

	for(int i = 0; i < inst->dimension_c; i++)
		iota(inst->index_sorted[i].begin(), inst->index_sorted[i].end(),0);

	for (int j = 0; j < inst->dimension_c; j++) {
		sort(inst->index_sorted[j].begin(), inst->index_sorted[j].end(),
			[&](int a, int b) { return inst->c_fc[a][j] < inst->c_fc[b][j]; });
	}
}

void write_log_CPIF(instance* inst,double solution_time, int cut_connect , int cut_cover) {
#ifdef to_file

	FILE* pFile;
	int n;
	char name[100];
	pFile = fopen("log.txt", "a");

	//get instance number
	string instance_number("NA");


	const regex base_regex("[\\D]+([\\d]+)[\\D]+");
	const string input = string(inst->instance_name);
	smatch base_match;
	if (regex_match(input, base_match, base_regex)) {
		if (base_match.size() == 2) {
			std::ssub_match base_sub_match = base_match[1];
			instance_number = base_sub_match.str();
		}
	}

	int lpstat = CPXgetstat(inst->env_CPIF, inst->lp_CPIF);

	CPXgetbestobjval(inst->env_CPIF, inst->lp_CPIF, &inst->best_lb);

	int nodecount = CPXgetnodecnt(inst->env_CPIF, inst->lp_CPIF);

	double solution;
	if (inst->use_custom_solution) {
		solution = inst->custom_solution_value;
		lpstat = 999;
	}
	else {
		solution = inst->objval;
	}

	fprintf(pFile, "%s\t %s \t %f \t %f \t %s \t %s \t %s \t %f \t %f \t %f \t %f \t %f \t %d \t %d \t %d \t %d \t %f \t lpstatus: %d\n", inst->instance_name, instance_number.c_str(), inst->r, inst->R, inst->problem_type, inst->options, inst->benders, solution, inst->objval, inst->best_lb, inst->root_bound, inst->opened_facilities, inst->lazy_cuts_connectivity, inst->lazy_cuts_coverage, inst->user_cuts_coverage, nodecount, solution_time, lpstat);
	fclose(pFile);

#endif
}

void write_log_MPIF(instance* inst, double solution_time, int cut_connect, int cut_cost, int user_cut_cost) {
#ifdef to_file
	FILE* pFile;
	int n;
	char name[100];
	pFile = fopen("log.txt", "a");

	//get instance number
	string instance_number("NA");


	const regex base_regex("[\\D]+([\\d]+)[\\D]+");
	const string input=string(inst->instance_name);
	smatch base_match;
	if (regex_match(input, base_match, base_regex)) {
		if (base_match.size() == 2) {
			std::ssub_match base_sub_match = base_match[1];
			instance_number = base_sub_match.str();
		}
	}

	int lpstat = CPXgetstat(inst->env_MPIF, inst->lp_MPIF);

	CPXgetbestobjval(inst->env_MPIF, inst->lp_MPIF, &inst->best_lb);


	double solution;
	if (inst->use_custom_solution) {
		solution = inst->custom_solution_value;
		lpstat = 999;
	}
	else {
		solution = inst->objval;
	}

	int nodecount = CPXgetnodecnt(inst->env_MPIF,inst->lp_MPIF);

	fprintf(pFile, "%s\t %s \t %f \t %s \t %s \t %s \t %f \t %f \t %f \t %f \t %f \t %d \t %d \t %d \t %d \t %f \t lpstatus: %d\n", inst->instance_name, instance_number.c_str(), inst->r, inst->problem_type, inst->options, inst->benders, solution, inst->objval, inst->best_lb, inst->root_bound,inst->opened_facilities, inst->lazy_cuts_connectivity, inst->lazy_cuts_cost, inst->user_cuts_cost, nodecount, solution_time, lpstat);
	fclose(pFile);
#endif
}

void debug(instance* inst) {
	cout <<"Debugging..."<<endl;
	inst->env_CPIF = CPXopenCPLEX(&(inst->status));
	printf("CPLEX version is %s \n", CPXversion(inst->env_CPIF));
	printf("%f \n",(double)CLOCKS_PER_SEC);


	/*
	strcpy(inst->input_file, "pmed_test.txt");
	inst->R = 9999;
	strcpy(inst->instance_type, "pmed");
	READ_NEW_FILE(inst);
	calculate_sorted_c(inst);
	for (int i = 0; i < inst->dimension; i++) {
		for (int j = 0; j < inst->dimension; j++) {
			cout<<inst->index_sorted[i][j]<<" ";
		}
		cout<<endl;
	}
	*/
}

void draw_tex_MPIF(instance* inst, double *x, vector<bool>* uncovered)
{

	vector<int> assign;
	for (int i = 0; i < inst->dimension_f; i++) {
		assign.push_back(-1);
	}

	//for (int i = 0; i < inst->dimension; i++) {
	//	if (x[inst->z[i]] > 1 - epsilon) {
	//		double min_dist = INFINITY;
	//		for (int j = 0; j < inst->dimension; j++) {
	//			if (x[inst->y[j]] > 1 - epsilon && inst->c[i][j] < min_dist) {
	//				assign[i] = j;
	//				min_dist = inst->c[i][j];
	//			}
	//		}
	//	}
	//}

	//eucl paper instances
	double size_grid = 100;

	double thresh;

#ifdef relax_y
	thresh = 0.01;
#else
	thresh = 1 - epsilon;
#endif //

	typedef struct colors_struct{ int red = 0; int green = 0; int blue = 0; } col;
	map<int, col> colormap;

	double scale_factor = size_grid / 10;

	double point = size_grid / scale_factor + 1;

	cout << "Generating LaTex picture ...\n";
	//cin.get();

	ofstream compact_file;
	if (uncovered != NULL) { compact_file.open("MPIF_instance_unreach.tex"); }
	else { compact_file.open("MPIF_instance.tex"); }

	compact_file
		<< "\\documentclass[]{article}" << "\n"
		<< "\\usepackage{tikz}" << "\n"
		<< "\\begin{document}" << "\n"
		<< "\\begin{tikzpicture}" << "\n"
		<< "\\usetikzlibrary{shapes}" << "\n"
		<< "\\tikzset{cross/.style={cross out, draw=black,inner sep=2pt, outer sep=2pt}}"
		<< endl;

	compact_file << "\\clip(0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;

	vector<int> colorcounter_red = { 0,   255, 255,   0,   0, 255,   0, 255, 192, 128, 128, 128,   0, 128,   0,   0 };
	vector<int> colorcounter_green = { 0, 255,   0, 255,   0, 255, 255,   0, 192, 128,   0, 128, 128,   0, 128,   0 };
	vector<int> colorcounter_blue = { 1,  255,   0,   0, 255,   0, 255, 255, 192, 128,   0,   0,   0, 128, 128, 128 };
	int color_counter = 0;

	for (int i = 0; i < inst->dimension_f; i++) {

		if (i == 0) {
			colormap.insert(pair<int, col>(i, { (colorcounter_red[color_counter % 16]) ,(colorcounter_green[color_counter % 16]) ,(colorcounter_blue[color_counter % 16]) }));
			color_counter++;

			compact_file <<
				"\\node[fill=black, regular polygon, regular polygon sides = 4, inner sep = 3pt] at (" <<
				round((inst->pos_f[i][0]) / scale_factor * 100) / 100 <<
				"," <<
				round((inst->pos_f[i][1]) / scale_factor * 100) / 100 <<
				"){};\n"
				<< endl;

		}
		else

			if (x[inst->y[i]] > thresh) {

				colormap.insert(pair<int, col>(i, { (colorcounter_red[color_counter % 16]) ,(colorcounter_green[color_counter % 16]) ,(colorcounter_blue[color_counter % 16]) }));
				color_counter++;

				compact_file <<
					"\\node[fill={rgb:red," << colormap.find(i)->second.red << " ;green," << colormap.find(i)->second.green << ";blue," << colormap.find(i)->second.blue << "}, regular polygon, regular polygon sides = 3, inner sep = 1.7pt] at (" <<
					round((inst->pos_f[i][0]) / scale_factor * 100) / 100 <<
					"," <<
					round((inst->pos_f[i][1]) / scale_factor * 100) / 100 <<
					"){};\n"
					<< endl;

			}

		if (strncmp(inst->options, "S", 1) == 0) {
			for (int j = 0; j < inst->dimension_f; j++) {
				if (x[inst->f[i][j]] > thresh) {
					compact_file << "\\draw[line width=0.01mm] (" << inst->pos_f[j][0] / scale_factor << "," << inst->pos_f[j][1] / scale_factor << ") -- (" << inst->pos_f[i][0] / scale_factor << "," << inst->pos_f[i][1] / scale_factor << ");\n" << endl;
				}
			}
		}
		else if (strncmp(inst->options, "M", 1) == 0) {
			for (int j = 0; j < inst->dimension_c; j++) {
				for (int k = 0; k < inst->dimension_f; k++) {
					if (x[inst->fm[k][i][j]] > thresh) {
						compact_file << "\\draw[line width=0.01mm] (" << inst->pos_f[j][0] / scale_factor << "," << inst->pos_f[j][1] / scale_factor << ") -- (" << inst->pos_f[i][0] / scale_factor << "," << inst->pos_f[i][1] / scale_factor << ");\n" << endl;
					}
				}
			}
		}

		if (uncovered != NULL) {
			if (uncovered->at(i)) {
				compact_file <<
					"\\draw[line width=0.01mm,fill={rgb:red," << 100 << " ;green," << 50 << ";blue," << 100 << "}] (" <<
					round((inst->pos_c[i][0]) / scale_factor * 100) / 100 <<
					"," <<
					round((inst->pos_c[i][1]) / scale_factor * 100) / 100 <<
					") circle (0.06);\n"
					<< endl;
			}
		}

	}



	for (int i = 0; i < inst->dimension_f; i++) {
		if (x[inst->y[i]] > thresh)
			for (int j = i + 1; j < inst->dimension_f; j++) {
				if (x[inst->y[j]] > thresh) {
					if (inst->c_f[i][j] <= inst->r) {
						compact_file << "\\draw[line width=0.01mm] (" << inst->pos_f[j][0] / scale_factor << "," << inst->pos_f[j][1] / scale_factor << ") -- (" << inst->pos_f[i][0] / scale_factor << "," << inst->pos_f[i][1] / scale_factor << ");\n" << endl;
					}
				}
			}
	}

	for (int i = 0; i < inst->dimension_c; i++) {

		compact_file << "\\node[] at(" << round((inst->pos_c[i][0]) / scale_factor * 100) / 100 + 0.1 << "," << round((inst->pos_c[i][1]) / scale_factor * 100) / 100 + 0.1 << ") {" << i << " };" << endl;

	}

	compact_file << "\\draw[ line width=0.7mm ] (0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;

	compact_file
		<< "\\end{tikzpicture}" << "\n"
		<< "\\end{document}" << "\n"
		<< endl;

	compact_file.close();




//	read_points(inst);
//
//	//eucl paper instances
//	int size_grid = 7000;
//
//	double thresh;
//
//#ifdef relax_y
//	thresh = 0.01;
//#else
//	thresh = 1 - epsilon;
//#endif //
//
//	typedef struct { int red = 0; int green = 0; int blue = 0; } col;
//	map<int,col> colormap;
//
//	double scale_factor = size_grid/10;
//
//	double point = size_grid / scale_factor;
//
//	cout << "Generating LaTex picture ...\n";
//	//cin.get();
//
//	ofstream compact_file;
//	if(uncovered!=NULL){ compact_file.open("MPIF_instance_unreach.tex"); }
//	else { compact_file.open("MPIF_instance.tex"); }
//
//	compact_file
//		<< "\\documentclass[]{article}" << "\n"
//		<< "\\usepackage{tikz}" << "\n"
//		<< "\\begin{document}" << "\n"
//		<< "\\begin{tikzpicture}" << "\n"
//		<< "\\usetikzlibrary{shapes}" << "\n"
//		<< endl;
//
//	compact_file << "\\clip(0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;
//
//	int colorcounter_red = 0;
//	int colorcounter_green = 0;
//	int colorcounter_blue = 0;
//
//	bool red_switch = true;
//	bool green_switch = false;
//	bool blue_switch = false;
//
//	for (int i = 0; i < inst->dimension; i++) {
//
//		if (i == 0) {
//			colormap.insert(pair<int, col>(i, { 256,0 ,0 }));
//
//			compact_file <<
//				"\\draw[line width=0.01mm,fill={rgb:red," << 256 << " ;green," << 0 << ";blue," << 0 << "}] (" <<
//				(inst->pos[i][0]) / scale_factor <<
//				"," <<
//				(inst->pos[i][1]) / scale_factor <<
//				") circle (0.10);\n"
//				<< endl;
//		}
//		else
//
//		if (x[inst->y[i]] > thresh) {
//
//			if (red_switch) {
//				red_switch = false;
//				blue_switch = true;
//				green_switch = false;
//				colorcounter_red += 80;
//			}
//			else if (blue_switch) {
//				red_switch = false;
//				blue_switch = false;
//				green_switch = true;
//				colorcounter_green += 80;
//			}
//			else if (green_switch) {
//				red_switch = true;
//				blue_switch = false;
//				green_switch = false;
//				colorcounter_blue += 80;
//			}
//
//			colormap.insert(pair<int,col>(i,{(colorcounter_red)%256,(colorcounter_green) % 256 ,(colorcounter_blue) % 256 }));
//
//			compact_file <<
//				"\\draw[line width=0.01mm,fill={rgb:red," << colormap.find(i)->second.red << " ;green," << colormap.find(i)->second.green << ";blue," << colormap.find(i)->second.blue << "},opacity =0.0001,fill opacity=0.0001] (" <<
//				round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
//				"," <<
//				round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
//				") circle (" << inst->r / scale_factor << ");\n"
//				<< endl;
//
//			compact_file <<
//				"\\node[fill={rgb:red," << colormap.find(i)->second.red << " ;green," << colormap.find(i)->second.green << ";blue," << colormap.find(i)->second.blue << "}, regular polygon, regular polygon sides = 3, inner sep = 1.7pt] at (" <<
//				round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
//				"," <<
//				round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
//				"){};\n"
//				<< endl;
//
//
//		}
//
//		if (strncmp(inst->options, "S", 1) == 0) {
//			for (int j = 0; j < inst->dimension; j++) {
//				if (x[inst->f[i][j]] > thresh) {
//					compact_file << "\\draw[line width=0.01mm] (" << inst->pos[j][0] / scale_factor << "," << inst->pos[j][1] / scale_factor << ") -- (" << inst->pos[i][0] / scale_factor << "," << inst->pos[i][1] / scale_factor << ");\n" << endl;
//				}
//			}
//		}
//		else if (strncmp(inst->options, "M", 1) == 0) {
//			for (int j = 0; j < inst->dimension; j++) {
//				for (int k = 0; k < inst->dimension; k++) {
//					if (x[inst->fm[k][i][j]] > thresh) {
//						compact_file << "\\draw[line width=0.01mm] (" << inst->pos[j][0] / scale_factor << "," << inst->pos[j][1] / scale_factor << ") -- (" << inst->pos[i][0] / scale_factor << "," << inst->pos[i][1] / scale_factor << ");\n" << endl;
//					}
//				}
//			}
//		}
//
//		if (uncovered != NULL) {
//			if (uncovered->at(i)) {
//				compact_file <<
//					"\\draw[line width=0.01mm,fill={rgb:red," << 100 << " ;green," << 50 << ";blue," << 100 << "}] (" <<
//					(inst->pos[i][0]) / scale_factor <<
//					"," <<
//					(inst->pos[i][1]) / scale_factor <<
//					") circle (0.06);\n"
//					<< endl;
//			}
//		}
//
//	}
//
//	for (int i = 0; i < inst->dimension; i++) {
//		if (x[inst->y[i]] < thresh) {
//			//color coded assignment
//			if (inst->options[1] == 'W') {
//				for (int j = 0; j < inst->dimension; j++) {
//					if (x[inst->w[i][j]] > thresh) {
//
//						compact_file <<
//							"\\draw[line width=0.01mm,fill={rgb:red," << colormap.find(j)->second.red << " ;green," << colormap.find(j)->second.green << ";blue," << colormap.find(j)->second.blue << "}] (" <<
//							round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
//							"," <<
//							round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
//							") circle (0.06);\n"
//							<< endl;
//					}
//				}
//			}
//			else {
//
//				compact_file <<
//					"\\draw[line width=0.01mm,fill={rgb:red," << 50 << " ;green," << 100 << ";blue," << 0 << "}] (" <<
//					round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
//					"," <<
//					round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
//					") circle (0.06);\n"
//					<< endl;
//			}
//		}
//		compact_file << "\\node[] at(" << round((inst->pos[i][0]) / scale_factor * 100) / 100 + 0.1 << "," << round((inst->pos[i][1]) / scale_factor * 100) / 100 + 0.1 << ") {" << i << " };" << endl;
//	}
//	compact_file << "\\draw[ line width=0.7mm ] (0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;
//
//	compact_file
//		<< "\\end{tikzpicture}" << "\n"
//		<< "\\end{document}" << "\n"
//		<< endl;
//
//	compact_file.close();

}

/*
void draw_tex_CPIF(instance* inst, double* x, vector<bool>* uncovered)
{

	vector<int> assign;
	for (int i = 0; i < inst->dimension; i++) {
		assign.push_back(-1);
	}

	for (int i = 0; i < inst->dimension; i++) {
		if (x[inst->z[i]] > 1 - epsilon) {
			double min_dist = INFINITY;
			for (int j = 0; j < inst->dimension; j++) {
				if (x[inst->y[j]] > 1 - epsilon && inst->c[i][j] < min_dist) {
					assign[i] = j;
					min_dist = inst->c[i][j];
				}
			}
		}
	}

	read_points(inst);

	//eucl paper instances
	double size_grid = 7000;

	for (int i = 0; i < inst->dimension; i++) {
		if (inst->pos[i][0] > size_grid) {
			size_grid = inst->pos[i][0];
		}
		if (inst->pos[i][1] > size_grid) {
			size_grid = inst->pos[i][1];
		}
	}

	double thresh;

#ifdef relax_y
	thresh = 0.01;
#else
	thresh = 1 - epsilon;
#endif //

	typedef struct color_struct2{ int red = 0; int green = 0; int blue = 0; } col;
	map<int, col> colormap;

	double scale_factor = size_grid / 10;

	double point = size_grid / scale_factor + 1;

	cout << "Generating LaTex picture ...\n";
	//cin.get();

	ofstream compact_file;
	if (uncovered != NULL) { compact_file.open("CPIF_instance_unreach.tex"); }
	else { compact_file.open("CPIF_instance.tex"); }

	compact_file
		<< "\\documentclass[]{article}" << "\n"
		<< "\\usepackage{tikz}" << "\n"
		<< "\\begin{document}" << "\n"
		<< "\\begin{tikzpicture}" << "\n"
		<< "\\usetikzlibrary{shapes}" << "\n"
		<<"\\tikzset{cross/.style={cross out, draw=black,inner sep=2pt, outer sep=2pt}}"
		<< endl;

	compact_file << "\\clip(0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;

 	vector<int> colorcounter_red = { 0,   255, 255,   0,   0, 255,   0, 255, 192, 128, 128, 128,   0, 128,   0,   0 };
	vector<int> colorcounter_green = { 0, 255,   0, 255,   0, 255, 255,   0, 192, 128,   0, 128, 128,   0, 128,   0 };
	vector<int> colorcounter_blue = { 1,  255,   0,   0, 255,   0, 255, 255, 192, 128,   0,   0,   0, 128, 128, 128 };
	int color_counter = 0;

	for (int i = 0; i < inst->dimension; i++) {

		if (i == 0) {
			colormap.insert(pair<int, col>(i, { (colorcounter_red[color_counter % 16]) ,(colorcounter_green[color_counter % 16]) ,(colorcounter_blue[color_counter % 16]) }));
			color_counter++;

			compact_file <<
				"\\node[fill=black, regular polygon, regular polygon sides = 4, inner sep = 3pt] at (" <<
				round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
				"," <<
				round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
				"){};\n"
				<< endl;

		}
		else

			if (x[inst->y[i]] > thresh) {

				colormap.insert(pair<int, col>(i, { (colorcounter_red[color_counter % 16]) ,(colorcounter_green[color_counter % 16]) ,(colorcounter_blue[color_counter % 16]) }));
				color_counter++;

				compact_file <<
					"\\node[fill={rgb:red," << colormap.find(i)->second.red << " ;green," << colormap.find(i)->second.green << ";blue," << colormap.find(i)->second.blue << "}, regular polygon, regular polygon sides = 3, inner sep = 1.7pt] at (" <<
					round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
					"," <<
					round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
					"){};\n"
					<< endl;

			}

		if (strncmp(inst->options, "S", 1) == 0) {
			for (int j = 0; j < inst->dimension; j++) {
				if (x[inst->f[i][j]] > thresh) {
					compact_file << "\\draw[line width=0.01mm] (" << inst->pos[j][0] / scale_factor << "," << inst->pos[j][1] / scale_factor << ") -- (" << inst->pos[i][0] / scale_factor << "," << inst->pos[i][1] / scale_factor << ");\n" << endl;
				}
			}
		}
		else if (strncmp(inst->options, "M", 1) == 0) {
			for (int j = 0; j < inst->dimension; j++) {
				for (int k = 0; k < inst->dimension; k++) {
					if (x[inst->fm[k][i][j]] > thresh) {
						compact_file << "\\draw[line width=0.01mm] (" << inst->pos[j][0] / scale_factor << "," << inst->pos[j][1] / scale_factor << ") -- (" << inst->pos[i][0] / scale_factor << "," << inst->pos[i][1] / scale_factor << ");\n" << endl;
					}
				}
			}
		}

		if (uncovered != NULL) {
			if (uncovered->at(i)) {
				compact_file <<
					"\\draw[line width=0.01mm,fill={rgb:red," << 100 << " ;green," << 50 << ";blue," << 100 << "}] (" <<
					round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
					"," <<
					round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
					") circle (0.06);\n"
					<< endl;
			}
		}

	}

	for (int i = 0; i < inst->dimension; i++) {
		if (x[inst->y[i]] < thresh) {
			//color coded assignment
			if (inst->options[1] == 'Z') {
				{
					if (assign[i] != -1) {

						if(assign[i] == 0){
							compact_file <<
								"\\draw[line width=0.01mm,fill=black] (" <<
								round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
								"," <<
								round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
								") circle (0.06);\n"
								<< endl;
						}
						else {
							compact_file <<
								"\\draw[line width=0.01mm,fill={rgb:red," << colormap.find(assign[i])->second.red << " ;green," << colormap.find(assign[i])->second.green << ";blue," << colormap.find(assign[i])->second.blue << "}] (" <<
								round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
								"," <<
								round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
								") circle (0.06);\n"
								<< endl;
													}
					}
					else {
						compact_file <<
							"\\node[cross, rotate = 90] at (" <<
							round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
							"," <<
							round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
							"){};\n"
							<< endl;
					}
				}
			}
			else {

				compact_file <<
					"\\draw[line width=0.01mm,fill={rgb:red," << 50 << " ;green," << 100 << ";blue," << 0 << "}] (" <<
					round((inst->pos[i][0]) / scale_factor * 100) / 100 <<
					"," <<
					round((inst->pos[i][1]) / scale_factor * 100) / 100 <<
					") circle (0.06);\n"
					<< endl;
			}
		}

		//node labels
		//compact_file << "\\node[] at(" << (inst->pos[i][0]) / scale_factor + 0.1 << "," << (inst->pos[i][1]) / scale_factor + 0.1 << ") {" << i << " };" << endl;
	}

	for (int i = 0; i < inst->dimension; i++) {
		if(x[inst->y[i]] > thresh)
		for (int j = i + 1; j < inst->dimension; j++) {
			if (x[inst->y[j]] > thresh) {
				if (inst->c[i][j] <= inst->r) {
					compact_file << "\\draw[line width=0.01mm] (" << inst->pos[j][0] / scale_factor << "," << inst->pos[j][1] / scale_factor << ") -- (" << inst->pos[i][0] / scale_factor << "," << inst->pos[i][1] / scale_factor << ");\n" << endl;
				}
			}
		}
	}

	compact_file << "\\draw[ line width=0.7mm ] (0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;

	compact_file
		<< "\\end{tikzpicture}" << "\n"
		<< "\\end{document}" << "\n"
		<< endl;

	compact_file.close();

}
*/
