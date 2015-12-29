
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#include "constants.h"
#include "MathFunctions.h"
#include "PotentialFunctions.h"
#include "SurvivalFunctions.h"
#include "Crosssections.h"

#include "NeutralsClassISpectrum.h"

/**
        INWARD TRAVELING NEUTRALS FROM CLASS II IONS
*/

// intregral of neutral from claas II ions _outside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte1(double E, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());
/*
    if ( dr < r )
    {
        printf("NeutralsClassIISpectrumInwards_Inte1 error: dr < r \n");
        return NAN;
    }
*/
	double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double integrant;

    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// intregral of neutral from claas II ions _inside_ the cathode _inwards_
double NeutralsClassIISpectrumInwards_Inte2(double E, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    if ( giveCathodeRadius() < dr )
    {
        printf("NeutralsClassIISPectrumInwards error: d < dr \n");
        return NAN;
    }

	double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double integrant;

    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(dr, ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

double NeutralsClassIISpectrumInwards (double r, double E)
{
    double flux;
    double term1 = 0, term2 = 0;
	
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 /= pow(r,2);
    term1 *= ngas * CrosssecCX(E);

    if ( giveCathodeRadius() < r)
    {
        // outside the cathode region
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumInwards_Inte1;

        term1 *= NIntegration_2(*FunctPtr, E, r, giveAnodeRadius());

        flux = term1;
    }
    else
    {
        // inside the cathode region
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumInwards_Inte2;

        term1 *= NIntegration_2(FunctPtr, E, giveCathodeRadius(), giveAnodeRadius());

		double ddr = Potential_Phi_Inv(Potential_Phi(giveCathodeRadius() - E / giveq()));

        double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());
        
        term2  = pow(ddr/r, 2);
        term2 *= g(giveCathodeRadius(), ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr),2);
        term2 *= 1 - exp(ngas * CrosssecCX(E) * ( r - giveCathodeRadius()) );
        term2 *= interpolation(Table->S, ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = giveTransparency() * ( term1 + term2);
    }

    return flux;
}

/**
        OUTWARD TRAVELING NEUTRALS FROM CLASS II IONS
*/

// intregral of neutrals from claas II ions _outside_ the cathode _outward_
double NeutralsClassIISpectrumOutwards_Inte1(double E, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());
/*
    if ( dr < r )
    {
        printf("NeutralsClassIISpectrumOutwards_Inte1 error: dr < r \n");
        return NAN;
    }
*/
	double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    double integrant;

    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant *= g(giveCathodeRadius(), ddr);
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

// intregral of neutral from claas II ions _inside_ the cathode _outwards_
double NeutralsClassIISpectrumOutwards_Inte2(double E, double r, double ddr)
{
    double dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

    double integrant = 0;
    if ( giveCathodeRadius() < dr )
        integrant = g(dr, ddr);

    if ( dr < r )
        integrant += pow(g(0,ddr),2) / g(dr, ddr);

    if ( integrant == 0)
        return 0;

	double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;
	
    integrant  = interpolation(Table->S, ddr);
    integrant /= abs(differentiat(*PhiPtr, dr));
    integrant /= 1 - pow(giveTransparency() * g(0,ddr),2);
    integrant *= pow(ddr,2);

    return integrant;
}

double NeutralsClassIISpectrumOutwards (double r, double E)
{
    double flux, dr, ddr;
    double term1 = 0, term2 = 0;
	
    double (*PhiPtr)(double);
    PhiPtr = &Potential_Phi;

    term1  = 1 / giveq();
    term1 /= pow(r,2);
    term1 *= ngas * CrosssecCX(E);

    if ( r < giveCathodeRadius() )
    {
        // inside the cathode region

        double (*FunctPtr)(double, double, double);
        FunctPtr = &NeutralsClassIISpectrumOutwards_Inte2;

        term1 *= NIntegration_3(*FunctPtr, E, r, r, giveAnodeRadius());

		ddr = Potential_Phi_Inv(Potential_Phi(giveCathodeRadius() - E / giveq()));

        dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());

        term2  = pow(ddr/r, 2);
        term2 *= g(giveCathodeRadius(), ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr),2);
        term2 *= 1 - exp(- ngas * CrosssecCX(E) * ( r + giveCathodeRadius()) );
        term2 *= interpolation(Table->S, ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = giveTransparency() * ( term1 + term2);
    }
    else
    {
        // outside the cathode region
        double (*FunctPtr)(double, double);
        FunctPtr = &NeutralsClassIISpectrumOutwards_Inte1;

        term1 *= NIntegration_2(FunctPtr, E, giveCathodeRadius(), giveAnodeRadius());

        ddr = Potential_Phi_Inv(Potential_Phi(giveCathodeRadius() - E / giveq()));
		dr = Potential_Phi_Inv(Potential_Phi(ddr) - E/giveq());
        
        term2  = pow(ddr/r, 2);
        term2 *= g(giveCathodeRadius(), ddr);
        term2 /= 1 - pow(giveTransparency() * g(0, ddr),2);
        term2 *= 1 - exp(-2 * ngas * CrosssecCX(E) * giveCathodeRadius());
        term2 *= interpolation(Table->S, ddr) / ( giveq() * abs(differentiat(*PhiPtr, dr)));

        flux = pow(giveTransparency(),2) * ( term1 + term2);
    }

    return flux;
}
