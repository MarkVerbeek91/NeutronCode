#include <iostream>
#include <math.h>
#include <stdio.h>

#include "constants.hpp"
#include "functions.hpp"


using namespace std;

Fusor fusor;

int main()
{
    // filling the potential array and particle energy
    printf("-- Start of program -- \n");

    // writing the potential to a file for plotting

    if ( false )
    {
    printf("Potential calculation\n");

    double (*Potential_PhiPtr)(double);
    Potential_PhiPtr = &Potential_Phi;

    double (*ParticleEnergy2Ptr)(double,double);
    ParticleEnergy2Ptr = &ParticleEnergy2;

    print_data_ddd(*ParticleEnergy2Ptr, 0.0, fusor.b+0.001, 0.01, 0.0, "Particle2.csv");

    print_data_dd(*Potential_PhiPtr, 0.0, 0.25, 0.001, "Potential.csv");


    }

    // writing the SIIEE to a file for plotting
    if ( false )
    {
        printf("SIIEE calculation\n");

        double (*SIIEEPtr)(double);
        SIIEEPtr = &SIIEE;

        print_data_dd(*SIIEEPtr, 1.0, -fusor.V0, 1, "SIIEE.csv");
    }

    // writing the Cross sections for Charge Exchange, Iononisation and the sum
    // of those to files to a file for plotting
    if ( false )
    {
        printf("Cross section calculation\n");

        double (*CrosssecCXPtr)(double);
        CrosssecCXPtr = &CrosssecCX;

        print_data_dd(*CrosssecCXPtr, 1.0, 500000, 10, "CrosssecCX.csv");

        double (*CrosssecIonPtr)(double);
        CrosssecIonPtr = &CrosssecIon;

        print_data_dd(*CrosssecIonPtr, 1.0, 500000, 10, "CrosssecIon.csv");

        double (*CrosssecTotPtr)(double);
        CrosssecTotPtr = &CrosssecTot;

        print_data_dd(*CrosssecTotPtr, 1.0, 500000, 10, "CrosssecTot.csv");
    }


    // Writing the survival functions to a file for plotting.
    if ( false )
    {
        printf("Survival function calculation\n");

        double (*fPtr)(double);
        fPtr = &f;

        print_data_dd(*fPtr, 0.0, fusor.b+0.001, 0.001, "f.csv");

        double (*gPtr)(double,double);
        gPtr = &g;

        print_data_ddd(*gPtr, 0.0, fusor.b+0.001, 0.001, 0, "g.csv");
    }


    // writing A to a file for plotting
    if ( false )
    {
        printf("doing some thing with A\n");

        double (*APtr)(double);
        APtr = &A;

        print_data_dd(*APtr, 0.05, fusor.b+0.01, 0.001, "A.csv");
    }

    // building the "Kernel"
    if ( false )
    {
        printf("Kernel\n");
        double (*KPtr)(double,double);
        KPtr = &kernel;

        print_data_ddd(*KPtr, 0.0, fusor.b+0.001, 0.01, 0.0001, "K.csv");
    }

    // source rate for first generation of Class II ions.

    kernel_to_table();
    S();

    // print the S tables to screen
    if ( true )
    {
        printf("printing table to screen\n");

        print_table(1, "Atable.csv");
//      print_table(&Table.K, "Ktable.csv");
        print_table(3, "S1.csv");
        print_table(4, "S2.csv");
        print_table(5, "S3.csv");

    }

    // program is done
    printf("\n-- Done --");
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
        phi = (fusor.a * (fusor.b - r) * -1 * fusor.V0) / (r * (fusor.a - fusor.b));

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
double ParticleEnergy1(double r)
{
    double energy;
    energy = - Potential_Phi(r) + 0.001;
    return energy;
}

/** \brief Calculates the energy of particle at radius r coming from born radius r1
 *
 * \param r = radius where particle is
 * \param r1= radius where particle is born
 * \return energy in eV.
 *
 */
double ParticleEnergy2(double r, double r1)
{
    double energy;
    energy = (Potential_Phi(r1) - Potential_Phi(r)) + 4;
    return energy;
}

/** \brief Calculates the CX crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecCX(double E)
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
double CrosssecIon(double E)
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
double CrosssecTot(double energy)
{
    double crosssection = 0;
    crosssection = CrosssecCX(energy) + CrosssecIon(energy);
    return crosssection;
}

/** \brief Calculates the fusion crosssection for particle with energy E
 *
 * \param E = energy of particle
 * \param
 * \return crosssection in m2
 *
 */
double CrosssecFusion(double E)
{
    double crosssection, energy = E/1000.;
    crosssection = 1e-28 * (CS_Fusion.A5Fusion + (CS_Fusion.A2Fusion / (pow(CS_Fusion.A4Fusion - CS_Fusion.A3Fusion * energy,2)+1)));
    crosssection = crosssection /(energy * (exp(CS_Fusion.A1Fusion / sqrt(energy))-1));
    return crosssection;
}

/**
    The next function compute the chance that a particle survise to that radius
*/
double f(double r)
{
    double sum = CrosssecCX(ParticleEnergy1(0)) + CrosssecCX(ParticleEnergy1(r));
    double step = (fusor.b - r)/N_pres;

    // very simple intergration.
    for (r; r < fusor.b; r+=step )
    {
        sum += 2.0*CrosssecCX(ParticleEnergy1(r));
    }

    sum = sum * step / 2.0;

    sum *= ngas;

    sum = exp(-1 * sum);

    return sum;
}

double Intergrant(double r)
{
    double tmp;

    tmp = ngas * CrosssecCX(ParticleEnergy1(r));

    return tmp;

}

double A(double r)
{
    double tmp = -1;
    tmp = ngas * CrosssecTot(ParticleEnergy1(r)) * gamma(r);

    if (isnan(tmp))
    {
        printf("went to the end\n");
        tmp = 0;
    }
    return tmp;
}

double gamma(double r)
{
    double Gamma = -1;
    Gamma = pow(fusor.b/r,2) * (f(r) + pow(fusor.Tc * f(0),2)/f(r));
    return Gamma;
}

double g(double r, double r1)
{
    if ( r > r1)
    {
        printf("error: r >= r1");
        return -2;
    }

    double sum = CrosssecCX(ParticleEnergy2(r,r1))+CrosssecCX(ParticleEnergy2(r1,r1)), step = (r1 - r)/N_pres;

    for (double r2=r; r2<r1; r2+= step)
    {
        sum += 2.0 * CrosssecCX(ParticleEnergy2(r2,r1));
 //       printf("r: %E, r1: %E, r2: %E tmp: %E\n",r,r1,r2,tmp);
    }

    sum *= ngas;
    sum = sum * step / 2.0;

 //   printf("integration done\n");

    sum = exp(-sum);
    return sum;
}

double kernel(double r, double r1)
{
    double tmp;
    if ( r < r1 )
    {
        tmp = pow(r1/r,2) * ((g(r,r1) + (pow(fusor.Tc*g(0,r1),2)/g(r,r1)))/(1.0-pow(fusor.Tc*g(0,r1),2)));
        tmp = tmp * ngas * CrosssecTot(ParticleEnergy2(r,r1));
    }
    else
        tmp = 0;

    return tmp;
}

void kernel_to_table(void)
{
    double step = (fusor.b - fusor.a)/(N_TABLE-1);

    for ( int i = 0; i < N_TABLE; i++)
    {
        for ( int j = 0; j < N_TABLE; j++)
        {
            Table.K[i][j] = kernel(i*step+fusor.a,j*step+fusor.a);
        }


        Table.A[i] = A(i*step+fusor.a);

        if ((N_TABLE / (i+1)) % 5 == 0)
            printf(".");
    }

    printf("\n");
/*
    for ( int i = 0; i < N_TABLE; i++)
    {
        printf("r: %f, A: %E, K: %E\n",i*step+fusor.a, Table.A[i],Table.K[0][i]);
    }
*/
    return;
}

void S_1(int r)
{
    // integral over Kernel times A from zero to fusor.b

    double sum = Table.K[r][0] + Table.K[r][N_TABLE-1];
    double step = (fusor.b-fusor.a)/N_pres;

    for (int i = 0; i<N_TABLE; i++)
    {
        sum += 2.0 * Table.K[r][i] * Table.A[i];
    }

    Table.S_1[r] = sum * (fusor.b-fusor.a) / (2.0 * N_TABLE);

    return;
}

void S_2(int r)
{
    double sum = Table.K[r][0] + Table.K[r][N_TABLE-1];
    double step = (fusor.b-fusor.a)/N_pres;

    for (int i = 0; i<N_TABLE; i++)
    {
        sum += 2.0 * Table.K[r][i] * Table.S_1[i];
    }

    Table.S_2[r] = sum * (fusor.b-fusor.a) / (2.0 * N_TABLE);

    return;
}

void S_3(int r)
{
    double sum = Table.K[r][0] + Table.K[r][N_TABLE-1];
    double step = (fusor.b-fusor.a)/N_pres;

    for (int i = 0; i<N_TABLE; i++)
    {
        sum += 2.0 * Table.K[r][i] * Table.A[i];
    }

    Table.S_1[r] = sum * (fusor.b-fusor.a) / (2.0 * N_TABLE);

    return;
}


void S(void)
{
    for ( int i = 0; i < N_TABLE; i++)
        S_1(i);

    for ( int i = 0; i < N_TABLE; i++)
        S_2(i);

    for ( int i = 0; i < N_TABLE; i++)
        S_3(i);

    return;
}


/**
    This function writes a given function to a file.
*/
void print_data_dd(double (*funcPtr)(double), double Start, double End, double step, char name[])
{
    FILE * output;
    output = fopen(name,"w");
    for (double r = Start; r<End; r+=step)
    {
        fprintf (output, "%E,%E\n",r, (*funcPtr)(r));
    }

    printf("%s: writing done\n",name);

    fclose(output);
    return;
}

/**
    This function writes a given function to a file.
*/
void print_data_ddd(double (*funcPtr)(double,double), double Start, double End, double step, double sec, char name[])
{
    FILE * output;
    output = fopen(name,"w");

    double counter = 0;

    for (double r = Start; r<End; r+=step)
    {
        fprintf (output, "%E,%E\n",r, (*funcPtr)(sec, r));


    }

    printf("%s: writing done\n",name);

    fclose(output);
    return;
}

void print_table(int choice, char name[])
{
    FILE * output;
    output = fopen(name,"w");

    switch (choice)
    {
        case 1:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%d, %E\n",i,Table.A[i]);
        }
        break;
        case 2:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%d, %E\n",i,Table.K[0][i]);
        }
        break;
        case 3:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%d, %E\n",i,Table.S_1[i]);
        }
        break;
        case 4:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%d, %E\n",i,Table.S_2[i]);
        }
        break;
        case 5:
        for ( int i = 0; i < N_TABLE; i++)
        {
            fprintf(output,"%d, %E\n",i,Table.S_3[i]);
        }

        break;


    }

    fclose(output);

    return;
}
