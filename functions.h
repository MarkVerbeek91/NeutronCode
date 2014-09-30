#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED


/**
    function declarations
*/
double Potential_Phi(double);
double SIIEE(double);
double ParticleEnergy1(float);
double ParticleEnergy2(float, float);
double CrosssecCX(float);
double CrosssecIon(float);
double CrosssecTot(float);
double CrosssecFusion(float);
double f(float);
double Intergrant(float);
double A(int);
double g(int, int);
double gamma(int);
double kernel(int, int);
void print_data(void);
void print_SIIEE(void);
void print_cross_section(void);
void print_kernel(void);




#endif // FUNCTIONS_H_INCLUDED
