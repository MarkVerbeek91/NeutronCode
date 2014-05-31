#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

struct Fusor{
    int a =  50; // cathode radius in mm
    int b = 250; // anode radius in mm
    int V0 = -55000; // voltage
    float wire_diameter = 0.5; // in mm

};

float q = 1.602e-19;
float ngas = 2.35404e19;

Fusor fusor;

/**
    each millimeter has it's own data point.
*/
struct data_struct{
    float phi[250];
    float ParticleEnergy[250];
    float SIIEE[100000];
    float Crosssec_CX[100000];
    float Crosssec_Ion[100000];
    float Crosssec_Tot[100000];
    float f[250];
    float A[250];
    float g[250][250];
} data;

/**
    Coefficients for the CX cross section of deuterium
*/
struct CX_cross_sections{
    float A_cx[6] = {3.245,235.88,0.03871,3.8068e-6,1.1832e-10,2.3713};
} CX_cs;

struct Ion_impact_cross_sections{
    float A_Ion[8] = {12.899,61.897,9.27e3,4.749e-4,3.989e-2,-1.59,3.1838,-3.7154};
} Ion_cs;

struct Fusion_cross_section{
    float A_fusion[5] = {47.88,482,3.08e-4,1.177,0};
} Fusion_cs;


/**
    function declarations
*/
float Potential_Phi(int);
float SIIEE(float);
float ParticleEnergy1(int);
float ParticleEnergy2(int, int);
float CrosssecCX(int);
float CrosssecIon(int);
float CrosssecTot(int);
float f(int);
void print_data(void);

int main()
{
    // filling the potential array and particle energy
    cout << " Start of program" << endl;
    for (int r=0; r <= 250; r++)
    {
        data.phi[r] = Potential_Phi(r);
        data.ParticleEnergy[r] = ParticleEnergy1(r);
    }

    // filling the SIIEE and crosssection data arrays
    cout << " Start of crosssection calculation " << endl;
    for (int E=0; E < -fusor.V0; E++)
    {
        data.SIIEE[E] = SIIEE(E);
        data.Crosssec_CX[E] = CrosssecCX(E);
        data.Crosssec_Ion[E] = CrosssecIon(E);
        data.Crosssec_Tot[E] = CrosssecTot(E);

        if ( E%1000 == 0)
            cout << ".";
    }

    // start of survival function calculation
    cout << "\n Start of crosssection calculation \n" << endl;
    for (int r=0; r <= 250; r++)
    {
        data.f[r] = f(r);
        cout << "." << endl;
    }

    int tmp;
    for (int r1=0; r1<250; r1++)
        for (int r=0; r < r1; r++)
            tmp = ParticleEnergy2(r,r1);

    print_data();

    cout << "Done" << endl;
    return 0;
}


/**
   This function returns the potential on position r. Input is a float r in
   meter and return is a float phi in volt.
*/
float Potential_Phi(int r)
{
    float phi;

    if ( r <= fusor.a)
        phi = fusor.V0;
    else
        phi = (fusor.a * (fusor.b - r) * -fusor.V0) / (r * (fusor.a - fusor.b));

    return phi;
}

float SIIEE(float energy)
{
    float tmp;
    tmp = 1.5*pow((1.15*(energy/97.861)),-0.667)*(1-exp(-1.8*pow(energy/97.891,1.2)));
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
float CrosssecCX(int energy)
{
    float crosssection;
    crosssection = (1e-20 * CX_cs.A_cx[0] * log( CX_cs.A_cx[2]/energy + CX_cs.A_cx[6])) / (1 + CX_cs.A_cx[3] * energy + CX_cs.A_cx[4] * pow(energy,3.5) + CX_cs.A_cx[5] * pow(energy,5.4));
    return crosssection;
}

float CrosssecIon(int energy)
{
    float crosssection = 0;

    return crosssection;
}

float CrosssecTot(int energy)
{
    float crosssection = 0;
    crosssection = data.Crosssec_CX[energy] + data.Crosssec_Ion[energy];
    return crosssection;
}

/**
    The next function compute the chance that a particle survise to that radius
*/
float f(int r)
{
    float tmp = 0;
    for ( int i=r; r = fusor.b; r++ )
        tmp = tmp + ngas * data.Crosssec_CX[r];

    tmp = exp(-tmp);

    return tmp;
}

float A(int r)
{
    float tmp = -1;

    return tmp;
}

float g(int r, int r1)
{
    float tmp = -1;

    return tmp;
}

/**
    This function writes the data to a file that's matlab readeble.
*/
void print_data(void)
{
    FILE * output;
    output = fopen("DATA.txt","w");
    for (int r=0; r<250; r++)
        fprintf (output, "%d,%E,%E,%E\n",r,data.phi[r],data.ParticleEnergy[r],f[r]);

    fclose(output);
    return;
}
