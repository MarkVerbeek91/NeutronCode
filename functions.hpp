
/**
    function declarations
*/
double Potential_Phi(double);
double SIIEE(double);
double ParticleEnergy1(double);
double ParticleEnergy2(double, double);
double CrosssecCX(double);
double CrosssecIon(double);
double CrosssecTot(double);
double CrosssecFusion(double);
double f(double);
double Intergrant(double);
double A(double);
double g(double, double);
double gamma(double);
double kernel(double, double);
double S_1(double);
double S_2(double);
double S_3(double);
void print_data_dd(double (*funcPtr)(double), double, double, double, char*);
void print_data_ddd(double (*funcPtr)(double,double), double, double, double, double, char*);
void print_SIIEE(void);
void print_cross_section(void);
void print_kernel(void);


