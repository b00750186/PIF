#ifndef MPIF_HEADER
#define MPIF_HEADER


#include "global_variables.h"
#include "global_functions.h"
#include "MPIF2.h"

using namespace std;

void add_pmed_constraint(instance* inst);
void add_singleflow_constraints_MPIF(instance* inst);
void add_multiflow_constraints_MPIF(instance* inst);
void initialize_v(instance* inst);

//in the paper
void build_model_MPIF_SW(instance* inst); //* scf with compact assignment

void build_model_MPIF_MW(instance* inst); //* mcf with compact assignment 

void build_model_MPIF_NW(instance* inst); // non-compact network with compact assignment
void solve_model_MPIF_NW(instance* inst);

void solve_model_MPIF(instance* inst); // nobenders and autobenders

void clean_model_MPIF(instance* inst);
#endif
