#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

struct Fusor{
    int a =  50; // cathode radius in mm
    int b = 250; // anode radius in mm
    int V0 = -55000; // voltage
    float wire_diameter = 0.5; // in mm
    float Tc;
};

Fusor fusor;

float q = 1.602e-19;
float pressure = 0.5;  // Pa
float Tgas = 300; // K
float ngas = 6.022e23 * pressure / (8.314 * Tgas);   //2.35404e19;
float E0 = 0.0001;          // reducing errors
float Itot = 0.1;

/**
    each millimeter has it's own data point.
*/
struct data_struct{
    float phi[250];
    float ParticleEnergy[250];
    float ParticleEnergy2[250][250];
    float SIIEE[100000];
    double Crosssec_CX[500000];
    double Crosssec_Ion[500000];
    double Crosssec_Tot[500000];
    double Crosssec_Fus[500000];
    float f[250];
    float A[250];
    float g[250][250];
    float Kernel[250][250];
} data;

/**
    Coefficients for the CX cross section of deuterium
*/
struct CX_cross_sections{
    float A1cx =   3.245     ;
    float A2cx = 235.88      ;
    float A3cx =   0.03871   ;
    float A4cx =   3.8068e-6 ;
    float A5cx =   1.1832e-10;
    float A6cx =   2.3713    ;
} CS_cx;

struct Ion_impact_cross_sections{
    float A1Ion = 12.899    ;
    float A2Ion = 61.897    ;
    float A3Ion =  9.27e3   ;
    float A4Ion =  4.749e-4 ;
    float A5Ion =  3.989e-2 ;
    float A6Ion = -1.59     ;
    float A7Ion =  3.1838   ;
    float A8Ion = -3.7154   ;
} CS_Ion;

struct Fusion_cross_section{
    float A1Fusion =  47.88     ;
    float A2Fusion = 482        ;
    float A3Fusion =   3.08e-4  ;
    float A4Fusion =   1.177    ;
    float A5Fusion =   0        ;
} CS_Fusion;


/**
    function declarations
*/
float Potential_Phi(float);
float SIIEE(float);
float ParticleEnergy1(int);
float ParticleEnergy2(int, int);
float CrosssecCX(int);
float CrosssecIon(int);
float CrosssecTot(int);
float CrosssecFusion(int);
float f(int);
float A(int);
float g(int, int);
float gamma(int);
float kernel(int, int);
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
    for (int r=0; r < 250; r++)
    {
        data.phi[r] = Potential_Phi(r);
        data.ParticleEnergy[r] = ParticleEnergy1(r);
        if ( r%25 == 0)
            cout << ".";
    }
    cout << endl;

    FILE * abc;
    abc = fopen("test.csv","w");
    for (int r1=0; r1 < 250; r1++)
    {
        for (int r=0; r < 250; r++)
        {
            data.ParticleEnergy2[r][r1] = ParticleEnergy2(r,r1);
            fprintf(abc,"%E",data.ParticleEnergy2[r][r1]);
        }
        fprintf(abc,"\n");
        if ( r1%25 == 0)
            cout << ".";
    }

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

    // start of survival function calculation
    cout << "\nSurvival function calculation" << endl;
    for (int r=0; r <= 250; r++)
    {
        data.f[r] = f(r);
   //     cout << data.f[r] << endl;
        data.A[r] = -1; // A(r);

        for (int r1=r; r1<= 250; r1++)
            data.g[r][r1] = g(0, r1);

        if ( r%25 == 0)
            cout << ".";
    }

    // filling of the kernel
    cout << "\nFilling of Kernel" << endl;
    for (int r=0; r < 250; r++)
    {
        for (int r1=0; r1 < 250; r1++)
          //  data.Kernel[r][r1] = kernel(r,r1);

        if ( r%25 == 0)
            cout << ".";
    }

    // printing output
    cout << "\nPrinting output to files" << endl;
    print_data();           cout << ".";
    print_SIIEE();          cout << ".";
    print_cross_section();  cout << ".";
    print_kernel();         cout << ".";

    cout << "\n-- Done --" << endl;
    return 0;
}


/**
   This function returns the potential on position r. Input is a float r in
   meter and return is a float phi in volt.
*/
float Potential_Phi(float r)
{
    float phi;

    if ( r <= fusor.a)
        phi = fusor.V0/1000;
    else
        phi = (fusor.a * (fusor.b - r) * -fusor.V0) / (r * (fusor.a - fusor.b) * 1000);

    return phi;
}

float SIIEE(float energy)
{
    float tmp;
    tmp = 1.5*pow((1.15*(energy/97861)),-0.667)*(1-exp(-1.8*pow(energy/97891,1.2)));
    return tmp;
}

float ParticleEnergy1(int r)
{
    float energy;
    energy = -0.5 * data.phi[r];
    return energy;
}

float ParticleEnergy2(int r, int r1)
{
    float energy;
//    energy = 0.5 * (Potential_Phi(r1) - Potential_Phi(r));
    energy = 0.5 * (data.phi[r1] - data.phi[r]);
    return energy;
}

/**
    The next three functions compute the cross section at given energy
*/
float CrosssecCX(int E)
{
    float crosssection;
    float energy = E/1000.;
    crosssection = 1e-20 * CS_cx.A1cx; // * log((CS_cx.A2cx/energy) + CS_cx.A6cx)) / (1 + CS_cx.A3cx * energy + CS_cx.A4cx * pow(energy,3.5) + CS_cx.A5cx * pow(energy,5.4));
//    cout << crosssection << endl;
    crosssection = crosssection * log((CS_cx.A2cx/energy) + CS_cx.A6cx);
//    cout << crosssection << endl;
    crosssection = crosssection / (1 + CS_cx.A3cx * energy + CS_cx.A4cx * pow(energy,3.5) + CS_cx.A5cx * pow(energy,5.4));
//    cout << "CX: " << crosssection << endl;

    return crosssection;
}

float CrosssecIon(int E)
{
    float crosssection = 0;
    float energy = E/1000.;
 //   crosssection = 1e-20 * CS_Ion.A1Ion * (exp(-CS_Ion.A2Ion/energy) * log(1 + CS_Ion.A3Ion * energy) / energy + CS_Ion.A4Ion * exp(-CS_Ion.A5Ion * energy) / (exp(CS_Ion.A6Ion) + CS_Ion.A7Ion*exp(CS_Ion.A8Ion)));
    crosssection = (exp(-CS_Ion.A2Ion/energy) * log(1 + CS_Ion.A3Ion * energy)) / energy;
//    cout << crosssection << endl;
    crosssection = crosssection + CS_Ion.A4Ion * exp(-CS_Ion.A5Ion * energy) / (exp(CS_Ion.A6Ion) + CS_Ion.A7Ion*exp(CS_Ion.A8Ion));
//    cout << crosssection << endl;
    crosssection = crosssection * 1e-20 * CS_Ion.A1Ion;
//    cout << "Ion: " <<crosssection << endl;
    return crosssection;
}

float CrosssecTot(int energy)
{
    float crosssection = 0;
    crosssection = data.Crosssec_CX[energy] + data.Crosssec_Ion[energy];
//    cout << "Tot: " << crosssection << endl;
    return crosssection;
}

float CrosssecFusion(int E)
{
    double crosssection, energy = E/1000., tmp;

    crosssection = 1e-28 * (CS_Fusion.A5Fusion + (CS_Fusion.A2Fusion / (pow(CS_Fusion.A4Fusion - CS_Fusion.A3Fusion * energy,2)+1)));
 //   cout << "1: " << crosssection << endl;
    crosssection = crosssection /(energy * (exp(CS_Fusion.A1Fusion / sqrt(energy))-1));
 //   cout << "2: " << crosssection << endl;
    return crosssection;
}

/**
    The next function compute the chance that a particle survise to that radius
*/
float f(int r)
{
    float tmp = 0;
    for (r; r < fusor.b; r++ )
    {
        tmp = tmp + ngas * data.Crosssec_CX[r];
        //cout << tmp << endl;
    }

    tmp = exp(-tmp);

    return tmp;
}

float A(int r)
{
    float tmp = -1;
    int energy = data.ParticleEnergy[r];
    tmp = ngas * gamma(r) * data.Crosssec_Tot[energy];

    return -1;
}

float g(int r, int r1)
{
    float tmp = -1.0;
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

float gamma(int r)
{
    float Gamma = -1;
    Gamma = pow(fusor.b,2)/pow(r,2) * (data.f[r] + pow(fusor.Tc * data.f[0],2)/data.f[r]);
    return Gamma;
}

float kernel(int r, int r1)
{
    float tmp;
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
    This function writes the data to a file that's matlab readeble.
*/
void print_data(void)
{
    FILE * output;
    output = fopen("DATA.csv","w");
    for (int r=1; r<249; r++)
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
    for (int r=0; r <= 250; r++)
    {
        for (int r1=r; r1<= 250; r1++)
            fprintf(output,"%E,",data.Kernel[r][r1]);

        fprintf(output,"\n");
    }
    return;
}
