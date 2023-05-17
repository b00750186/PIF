#ifndef CPIF_HEADER
#define CPIF_HEADER

#include "global_variables.h"
#include "global_functions.h"

using namespace std;
void add_pmed_constraint_CPIF(instance* inst);

void build_model_CPIF_SZ(instance* inst);
void build_model_CPIF_MZ(instance* inst);
void build_model_CPIF_NZ(instance* inst);

void solve_model_CPIF(instance* inst);
void solve_model_CPIF_NZ(instance* inst);

void clean_model_CPIF(instance* inst);
#endif