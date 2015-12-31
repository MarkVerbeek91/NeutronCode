#ifndef AUXILIARY_H_INCLUDED
#define AUXILIARY_H_INCLUDED

void output_data(void);
void plot_function_1D( double (*funcPtr)(double), double Start, double End, double step, const char name[], const char input_file_name[]);
void plot_function_2D(double (*funcPtr)(double, double), double Start1, double End1, double Start2, double End2, double step1, double step2, const char name[], const char input_file_name[]);
void print_data_dd( double (*funcPtr)(double), double Start, double End, double step, char name[], int ID);
void print_data_ddd(double (*funcPtr)(double,double), double Start, double End, double step, double sec, char name[], int ID);
void plot_table_1D(double *table, const char name[], const char input_file_name[]);
void plot_table_2D(double (*table)[N_TABLE], const char name[], const char input_file_name[]);
void print_table(int choice, char name[]);
void readfile(FILE * input);
bool readline_ini_file_bools(FILE* input);
void print_program_parameters(void);

#endif // AUXILIARY_H_INCLUDED
