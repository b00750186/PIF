#include "global_variables.h"
#include "global_functions.h"
#include "MPIF.h"
#include "MPIF2.h"
#include "CPIF.h"
#include "CPIF2.h"

#ifndef UNICODE
#define UNICODE
#endif 

//#define print_c


int main(int argc, char** argv) {
	
	instance inst;
	inst.input_file = (char*)calloc(150, sizeof(char));
	inst.problem_type = (char*)calloc(10, sizeof(char));
	inst.instance_type = (char*)calloc(10, sizeof(char));
	inst.benders = (char*)calloc(20, sizeof(char));
	inst.options = (char*)calloc(30, sizeof(char));
	inst.posfile = (char*)calloc(150, sizeof(char));
	
	if (argc == 9)
	{
		strcpy(inst.input_file, argv[1]);
		strcpy(inst.instance_type, argv[2]);
		strcpy(inst.problem_type, argv[3]);
		strcpy(inst.options, argv[4]);
		inst.timelimit = atoi(argv[5]);
		inst.r = atof(argv[6]);
		inst.R = atof(argv[7]);
		strcpy(inst.benders,argv[8]);
	}
	else if(argc == 10){
		strcpy(inst.input_file, argv[1]);
		strcpy(inst.instance_type, argv[2]);
		strcpy(inst.problem_type, argv[3]);
		strcpy(inst.options, argv[4]);
		inst.timelimit = atoi(argv[5]);
		inst.r = atof(argv[6]);
		inst.R = atof(argv[7]);
		strcpy(inst.benders, argv[8]);
		strcpy(inst.posfile, argv[9]);
	}
	else{
		cout << "This program is intended to solve the MPIF and the CPIF for euclidean and p-median instances:" << endl;
		cout << "List of parameters:" << endl;
		cout << "Parameter 1:\t path to file with facility and customers data\n";
		cout << "Parameter 2:\t instance type: eucl - euclidean, pmed - p-median\n";
		cout << "Parameter 3:\t MPIF - median, CPIF - covering\n";
		cout << "Parameter 4:\t Model options: SZ,MZ,NZ,NT;SW,MW,NW,NV\n";
		cout << "Parameter 5:\t time limit in seconds\n";
		cout << "Parameter 6:\t radius r\n";
		cout << "Parameter 7:\t radius R\n";
		cout << "Parameter 8:\t Solution strategy: autobenders, nobenders (has no effect for the NV, NW, NZ, NT models)\n";
		cout << "Parameter 9 (optional)\t File with coordinates for LaTex output generation\n";
		exit(-1);
	}

	read_input(&inst);

#ifdef print_c
	for (int i = 0; i < inst.dimension; i++) {
		for (int j = 0; j < inst.dimension; j++) {
			cout << inst.c[i][j] << " ";
		}
		cout << endl;
	}
#endif // print_c

	cout << "Solving instance: " << inst.instance_name<<" "<<inst.problem_type <<" "<< inst.options << " " << inst.benders << endl;

	if (strncmp(inst.problem_type, "MPIF", 4) == 0) {

		//model type selection
		if (strncmp(inst.options, "SW", 2) == 0) {
			build_model_MPIF_SW(&inst);
		}
		
		//multicommodity
		if (strncmp(inst.options, "MW", 2) == 0) {
			build_model_MPIF_MW(&inst);
		}

		//node separator and v variables model
		if (strncmp(inst.options, "NV", 2) == 0) {
			build_model_MPIF_NV(&inst);
		}
		
		//node separator
		if (strncmp(inst.options, "NW", 2) == 0) {
			build_model_MPIF_NW(&inst);
		}

		//solver algorithm selection
		if (strncmp(inst.options, "NV", 2) == 0) {
			solve_model_MPIF_NV(&inst);
		}
		else
		if (strncmp(inst.options, "NW", 2) == 0) {
			solve_model_MPIF_NW(&inst);

		}
		else {			
			if (strncmp(inst.benders, "nobenders", 4) == 0) {
				solve_model_MPIF(&inst);
			}
			if (strncmp(inst.benders, "autobenders", 4) == 0) {
				solve_model_MPIF(&inst);
			}
			if (strncmp(inst.benders, "silly", 4) == 0) {
				cout << inst.instance_name << "  r = "<<inst.r << " cost to connect all to the root = " << calculate_only_root_case(&inst);
			}
		}		
		
	}

	if (strncmp(inst.problem_type, "CPIF", 4) == 0) {
		if (strncmp(inst.options, "SZ", 2) == 0) {
			build_model_CPIF_SZ(&inst);
		}
		if (strncmp(inst.options, "MZ", 2) == 0) {
			build_model_CPIF_MZ(&inst);
		}
		if (strncmp(inst.options, "NZ", 2) == 0) {
			build_model_CPIF_NZ(&inst);
		}
		
		if (strncmp(inst.options, "NT", 2) == 0) {
			build_model_CPIF_NT(&inst);
		}

		if (strncmp(inst.options, "NT", 2) == 0) {
			solve_model_CPIF_NT(&inst);
		}else
		if (strncmp(inst.options, "NZ", 2) == 0) {
			solve_model_CPIF_NZ(&inst);
		}else{

			if (strncmp(inst.benders, "nobenders", 4) == 0) {
				solve_model_CPIF(&inst);
			}

			if (strncmp(inst.benders, "autobenders", 4) == 0) {
				solve_model_CPIF(&inst);
			}
		}

		clean_model_CPIF(&inst);
	}
	
	return 0;
}