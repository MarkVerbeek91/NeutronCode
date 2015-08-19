#ifndef AUXILIARY_H_INCLUDED
#define AUXILIARY_H_INCLUDED

void plot_function_dd( double (*funcPtr)(double), double Start, double End, double step, const char name[], const char input_file_name[]);
void plot_function_ddd(double (*funcPtr)(double, double), double var, double Start, double End, double step, const char name[], const char input_file_name[]);
void print_data_dd( double (*funcPtr)(double), double Start, double End, double step, char name[], int ID);
void print_data_ddd(double (*funcPtr)(double,double), double Start, double End, double step, double sec, char name[], int ID);
void plot_table_1D(const char name[], const char input_file_name[]);
void plot_table_2D();
void print_table(int choice, char name[]);
void readfile(FILE * input);

#endif // AUXILIARY_H_INCLUDED
