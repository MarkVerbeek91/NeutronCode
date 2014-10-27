/**
    This file contains all prototypes and contands needed for cross section
    calculation.

*/



double CrosssecCX(double);
double CrosssecIon(double);
double CrosssecTot(double);
double CrosssecFusion(double);

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