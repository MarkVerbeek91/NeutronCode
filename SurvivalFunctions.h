#ifndef SURVIVALFUNCTIONS_H_INCLUDED
#define SURVIVALFUNCTIONS_H_INCLUDED

double f_inte(double r);
double f(double r);
double Intergrant(double r);
double gamma(double r);
double A(double r);
double g_inte(double r, double dr);
double g(double r, double r1);
double kernel(double r, double r1);
void kernel_to_table(void);


#endif // SURVIVALFUNCTIONS_H_INCLUDED
