#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

struct Fusor{
    int a =  500; // cathode radius in um
    int b = 2500; // anode radius in um
    int V0 = -55000; // voltage
    double wire_diameter = 0.5; // in mm
    double Tc;
};

Fusor fusor;

// the precision of the functions.
#define N_pres 2500

double q = 1.602e-19;
double pressure = 0.5;  // Pa
double Tgas = 400; // K
double ngas = 9.05401e19; //6.022e23 * pressure / (8.314 * Tgas);
double E0 = 0.0001;          // reducing errors
double Itot = 0.1;

/**
    each millimeter has it's own data point.
*/
struct data_struct{
    double phi[N_pres];
    double ParticleEnergy[N_pres];
    double ParticleEnergy2[N_pres][N_pres];
    double SIIEE[100000];
    double Crosssec_CX[500000];
    double Crosssec_Ion[500000];
    double Crosssec_Tot[500000];
    double Crosssec_Fus[500000];
    double f[N_pres];
    double A[N_pres];
    double g[N_pres][N_pres];
    double Kernel[N_pres][N_pres];
} data;

/**
    Coefficients for the CX cross section of deuterium
*/
struct CX_cross_sections{
    double A1cx =   3.245     ;
    double A2cx = 235.88      ;
    double A3cx =   0.038371  ;
    double A4cx =   3.8068e-6 ;
    double A5cx =   1.1832e-10;
    double A6cx =   2.3713    ;
} CS_cx;

struct Ion_impact_cross_sections{
    double A1Ion = 12.899    ;
    double A2Ion = 61.897    ;
    double A3Ion =  9.27e3   ;
    double A4Ion =  4.9749e-4;
    double A5Ion =  3.989e-2 ;
    double A6Ion = -1.59     ;
    double A7Ion =  3.1834   ;
    double A8Ion = -3.7154   ;
} CS_Ion;

struct Fusion_cross_section{
    double A1Fusion =  47.88     ;
    double A2Fusion = 482        ;
    double A3Fusion =   3.08e-4  ;
    double A4Fusion =   1.177    ;
    double A5Fusion =   0        ;
} CS_Fusion;


/**
    function declarations
*/
double Potential_Phi(double);
double SIIEE(double);
double ParticleEnergy1(int);
double ParticleEnergy2(int, int);
double CrosssecCX(int);
double CrosssecIon(int);
double CrosssecTot(int);
double CrosssecFusion(int);
double f(int);
double Intergrant(int);
double A(int);
double g(int, int);
double gamma(int);
double kernel(int, int);
void print_data(void);
void print_SIIEE(void);
void print_cross_section(void);
void print_kernel(void);

int main()
{
    fusor.Tc = 1 - (2.5*fusor.wire_diameter *(fusor.a/pow(fusor.a,2)));

    // filling the potential array and particle energy
    cout << "-- Start of program --" << endl;
    cout << "Potential and particle energy calculation" << endl;
    for (int r=0; r < N_pres; r++)
    {
        data.phi[r] = Potential_Phi(r);
        data.ParticleEnergy[r] = ParticleEnergy1(r);
        if ( r%250 == 0)
            cout << ".";
    }
    cout << endl;

    FILE * abc;
    abc = fopen("ParticalEnergy2.csv","w");
    for (int r1=0; r1 < N_pres; r1++)
    {
        for (int r=0; r < N_pres; r++)
        {
            data.ParticleEnergy2[r][r1] = ParticleEnergy2(r,r1);
            fprintf(abc,"%E",data.ParticleEnergy2[r][r1]);
        }
        fprintf(abc,"\n");
        if ( r1%250 == 0)
            cout << ".";
    }
    fclose(abc);


    // filling the SIIEE data array
    cout << "\nSIIEE calculation" << endl;
    for (int E=0; E < -fusor.V0; E++)
    {
        data.SIIEE[E] = SIIEE(E);

        if ( E%(fusor.V0/50) == 0)
            cout << ".";
    }

  //  for (int E=0; E< 500; E++)
  //      data.Crosssec_Fus[E] = CrosssecFusion(E);

    // filling the crosssections data arrays
    cout << "\nCrosssection calculation" << endl;
    for (int E=0; E < 500000; E++)
    {
        data.Crosssec_CX[E]  = CrosssecCX(E);
        data.Crosssec_Ion[E] = CrosssecIon(E);
        data.Crosssec_Tot[E] = CrosssecTot(E);
    //    data.Crosssec_Fus[E] = CrosssecFusion(E);
        if (E < 50000)
            data.Crosssec_Fus[E] = 0;
        else
            data.Crosssec_Fus[E] = CrosssecFusion(E);

        if ( E%10000 == 0)
            cout << ".";
    }

    FILE * inte;
    inte = fopen("integrant1.csv","w");
    for (int r=0; r <= N_pres; r++)
    {
        Intergrant(r);
        fprintf(inte,"%i,%E\n", r, Intergrant(r));
    }

    fclose(inte);


    // start of survival function calculation
    cout << "\nSurvival function calculation" << endl;
    for (int r=0; r <= N_pres; r++)
    {
        data.f[r] = f(r);
 //       cout << "f(" << r << ") = " << data.f[r] << endl;
        data.A[r] = -1; // A(r);

        for (int r1=r; r1<= N_pres; r1++)
            data.g[r][r1] = -1; // g(0, r1);

        if ( r%250 == 0)
            cout << ".";
    }
/*
    // filling of the kernel
    cout << "\nFilling of Kernel" << endl;
    for (int r=0; r < 250; r++)
    {
        for (int r1=0; r1 < 250; r1++)
          //  data.Kernel[r][r1] = kernel(r,r1);

        if ( r%25 == 0)
            cout << ".";
    }
*/
    // printing output
    cout << "\nPrinting output to files" << endl;
    print_data();           cout << ".";
    print_SIIEE();          cout << ".";
    print_cross_section();  cout << ".";
    print_kernel();         cout << ".";

    cout << "\n-- Done --" << endl;
    return 0;
}

/** \brief This function returns the potential on position r. Input is a double r in
   meter and return is a double phi in KiloVolt.
 *
 * \param r = radius
 * \param
 * \return potential in V
 *
 */
double Potential_Phi(double r)
{
    double phi;

    if ( r <= fusor.a)
        phi = fusor.V0;
    else
        phi = 1000 * (fusor.a * (fusor.b - r) * -1 * fusor.V0/1000.) / (r * (fusor.a - fusor.b));

    return phi;
}

/**
    In E in eV
*/
double SIIEE(double energy)
{
    double tmp;
    tmp = 1.5*pow((1.15*(energy/97861)),-0.667)*(1-exp(-1.8*pow(energy/97891,1.2)));
    return tmp;
}

/** \brief calculates the energy of particle at radius r coming from outer edge
 *
 * \param r = radius where particle is
 * \return energie in eV
 *
 */
double ParticleEnergy1(int r)
{
    double energy;
    energy = - data.phi[r] + 1;
    return energy;
}

/** \brief Calculates the energy of particle at radius r coming from born radius r1
 *
 * \param r = radius where particle is
 * \param r1= radius where particle is born
 * \return
 *
 */
double ParticleEnergy2(int r, int r1)
{
    double energy;
    energy = 0.5 * (data.phi[r1] - data.phi[r]) + 1;
    return energy;
}

/** \brief Calculates the CX crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecCX(int E)
{
    double crosssection;
    double energy = E/1000.;
    crosssection = 1e-20 * CS_cx.A1cx;
    crosssection = crosssection * log((CS_cx.A2cx/energy) + CS_cx.A6cx);
    crosssection = crosssection / (1 + CS_cx.A3cx * energy + CS_cx.A4cx * pow(energy,3.5) + CS_cx.A5cx * pow(energy,5.4));

    return crosssection;
}

/** \brief Calculates the Ion crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecIon(int E)
{
    double crosssection = 0;

    // because the int E is in eV and the formula in keV, the factor 1/1000 is introduced.
    double energy = E/1000.;
    crosssection = (exp(-CS_Ion.A2Ion/energy) * log(1 + CS_Ion.A3Ion * energy)) / energy;
    crosssection = crosssection + CS_Ion.A4Ion * exp(-CS_Ion.A5Ion * energy) / (exp(CS_Ion.A6Ion) + CS_Ion.A7Ion*exp(CS_Ion.A8Ion));
    crosssection = crosssection * 1e-20 * CS_Ion.A1Ion;
    return crosssection;
}

/** \brief Calculates the Tot crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecTot(int energy)
{
    double crosssection = 0;
    crosssection = data.Crosssec_CX[energy] + data.Crosssec_Ion[energy];
    return crosssection;
}

/** \brief Calculates the fusion crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecFusion(int E)
{
    double crosssection, energy = E/1000.;
    crosssection = 1e-28 * (CS_Fusion.A5Fusion + (CS_Fusion.A2Fusion / (pow(CS_Fusion.A4Fusion - CS_Fusion.A3Fusion * energy,2)+1)));
    crosssection = crosssection /(energy * (exp(CS_Fusion.A1Fusion / sqrt(energy))-1));
    return crosssection;
}

/**
    The next function compute the chance that a particle survise to that radius
*/
double f(int r)
{
    double tmp = 0;
    int energy;

    // very simple intergration.
    for (r; r < fusor.b; r++ )
    {
        energy = data.ParticleEnergy[r];
        tmp = tmp + ngas * data.Crosssec_CX[energy] * 0.0001 ;
   //     cout << "r: " << r << " E:" << energy << " tmp: " << tmp << endl;
    }

 //   cout << -1 * tmp << endl;

    tmp = exp(-1 * tmp);

    return tmp;
}

double Intergrant(int r)
{
    double tmp;
    int energy;

    energy = data.ParticleEnergy[r];
    tmp = ngas * data.Crosssec_CX[energy];
 //   cout << "r = " << r << " E = " << energy << " Int: " << tmp << endl;

    return tmp;

}

double A(int r)
{
    double tmp = -1;
    int energy = data.ParticleEnergy[r];
    tmp = ngas * gamma(r) * data.Crosssec_Tot[energy];

    return -1;
}

double g(int r, int r1)
{
    double tmp = -1.0;
    int energy;

    if ( r > r1)
    {
        cout << "error: r >= r1" << endl;
        return -2;
    }

    for (int r2=r; r2<r1; r2++)
    {
        energy = data.ParticleEnergy2[r2][r];
        tmp = tmp + data.Crosssec_CX[energy];
    }

    tmp = ngas * tmp;
    tmp = exp(-tmp);
    return tmp;
}

double gamma(int r)
{
    double Gamma = -1;
    Gamma = pow(fusor.b,2)/pow(r,2) * (data.f[r] + pow(fusor.Tc * data.f[0],2)/data.f[r]);
    return Gamma;
}

double kernel(int r, int r1)
{
    double tmp;
     if ( r < r1 )
    {
        int energy;
        energy = data.ParticleEnergy2[r][r1];
        cout << energy << endl;
        tmp = pow(r1/r,2) * (data.g[r][r1] + pow(fusor.Tc,2)*data.g[0][r1]/data.g[r][r1])/(1-pow(fusor.Tc,2)*data.g[0][r1]);
        cout << tmp << endl;
        tmp = tmp * ngas * data.Crosssec_Tot[energy];
        cout << tmp << endl;
    }
    else
        tmp = 0;



    return tmp;
}

/**
    These functions write the data to files that are matlab readable.
*/
void print_data(void)
{
    FILE * output;
    output = fopen("DATA.csv","w");
    for (int r=1; r<N_pres; r++)
        fprintf (output, "%d,%E,%E,%E,%E,%E\n",r,data.phi[r],data.ParticleEnergy[r],data.f[r],data.g[0][r],data.A[r]);

    fclose(output);
    return;
}

void print_SIIEE(void)
{
    FILE * output;
    output = fopen("DATA_SIIEE.csv","w");
    for (int E=1; E < -fusor.V0; E++)
        fprintf (output, "%d,%E\n",E,data.SIIEE[E]);

    fclose(output);
    return;
}

void print_cross_section(void)
{
    FILE * output;
    output = fopen("DATA_cross_sections.csv","w");
    for (int E=1; E < 500000; E++)
        fprintf (output, "%d,%E,%E,%E,%E\n",E,data.Crosssec_CX[E],data.Crosssec_Ion[E],data.Crosssec_Tot[E],data.Crosssec_Fus[E]);

    fclose(output);
    return;
}

void print_kernel(void)
{
    FILE * output;
    output = fopen("Kernel.csv","w");
    for (int r=0; r <= N_pres; r++)
    {
        for (int r1=r; r1<= N_pres; r1++)
            fprintf(output,"%E,",data.Kernel[r][r1]);

        fprintf(output,"\n");
    }
    return;
}
