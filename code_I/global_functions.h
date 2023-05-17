#ifndef FUNCTIONS_local_HEADER
#define FUNCTIONS_local_HEADER

#include "global_variables.h"

void read_input(instance* inst);

void floyd(instance* inst);

void print_solution_CPIF(instance* inst);

void print_solution_MPIF(instance* inst);

void calculate_sorted_c(instance* inst);

void write_log_CPIF(instance* inst, double solution_time, int cut_connect = 0, int cut_cover = 0);

void write_log_MPIF(instance* inst, double solution_time, int cut_connect = 0, int cut_cost = 0, int user_cut_cost = 0);

void free_and_null(char** ptr);

void draw_tex_MPIF(instance* inst, double* x, vector<bool> *uncovered = NULL);

void draw_tex_CPIF(instance* inst, double* x,vector<bool> *uncovered = NULL);

void read_points(instance* inst);


#endif
