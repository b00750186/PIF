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
	inst.input_file = (char*)calloc(1000, sizeof(char));
	inst.problem_type = (char*)calloc(10, sizeof(char));
	inst.benders = (char*)calloc(20, sizeof(char));
	inst.options = (char*)calloc(30, sizeof(char));
	inst.posfile = (char*)calloc(1000, sizeof(char));

	if (argc == 8)
	{
		strcpy(inst.input_file, argv[1]);
		strcpy(inst.problem_type, argv[2]);
		strcpy(inst.options, argv[3]);
		inst.timelimit = atoi(argv[4]);
		inst.r = atof(argv[5]);
		inst.R = atof(argv[6]);
		strcpy(inst.benders,argv[7]);
	}
	else{
		cout << "This program is intended to solve the MPIF and the CPIF for euclidean and p-median instances:" << endl;
		cout << "List of parameters:" << endl;
		cout << "Parameter 1:\t path to file with facility and customers data" << endl;
		cout << "Parameter 2:\t MPIF - median, CPIF - covering" << endl;
		cout << "Parameter 3:\t Model options: SZ,MZ,NZ,NT;SW,MW,NW,NV" << endl;
		cout << "Parameter 4:\t time limit in seconds" << endl;
		cout << "Parameter 5:\t radius r" << endl;
		cout << "Parameter 6:\t radius R" << endl;
		cout << "Parameter 7:\t Solution strategy: autobenders, nobenders (has no effect for the NV, NW, NZ, NT models)" << endl;
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

	cout << "Solving instance: " << inst.instance_name << " " << inst.problem_type << " " << inst.options << " " << inst.benders << endl;

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
		}

	}

	if (strncmp(inst.problem_type, "CPIF", 4) == 0) {
		//inst.alpha = 1.0;
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
