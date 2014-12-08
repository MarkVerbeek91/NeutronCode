#ifndef AUXILIARY_H_INCLUDED
#define AUXILIARY_H_INCLUDED

void print_data_dd(double (*funcPtr)(double), double Start, double End, double step, char name[]);
void print_data_ddd(double (*funcPtr)(double,double), double Start, double End, double step, double sec, char name[]);
void print_data_ddd(double (*funcPtr)(double,double), double Start, double End, double step, double sec, char name[]);
void print_table(int choice, char name[]);
void readfile(FILE ** input);

#endif // AUXILIARY_H_INCLUDED
