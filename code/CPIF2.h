#ifndef CPIF_HEADER_2
#define CPIF_HEADER_2

#include "global_variables.h"
#include "global_functions.h"

double calculate_only_root_case_CPIF(instance* inst);
double calculate_2_facility_sol_CPIF(instance* inst);
void add_cardinal_constraint_CPIF(instance* inst);
void solve_model_CPIF_NT(instance* inst);
void build_model_CPIF_NT(instance* inst);

#endif
