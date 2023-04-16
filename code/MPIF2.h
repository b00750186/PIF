#ifndef MPIF_HEADER_2
#define MPIF_HEADER_2

#include "MPIF.h"
#include "global_variables.h"
#include "global_functions.h"

double calculate_only_root_case(instance *inst);
double calculate_only_2(instance* inst);
void add_cardinal_constraint(instance* inst);

void solve_model_MPIF_NV(instance* inst);
void build_model_MPIF_NV(instance* inst);

#endif
