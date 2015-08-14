/**
    This file contains all prototypes and contands needed for cross section
    calculation.

*/

#ifndef CROSSSECTIONSCONT_H_INCLUDED
#define CROSSSECTIONSCONT_H_INCLUDED

/**
    Coefficients for the CX cross section of deuterium
*/
struct CX_cross_sections{
    const double A1cx =   3.245     ;
    const double A2cx = 235.88      ;
    const double A3cx =   0.038371  ;
    const double A4cx =   3.8068e-6 ;
    const double A5cx =   1.1832e-10;
    const double A6cx =   2.3713    ;
} CS_cx;

struct Ion_impact_cross_sections{
    const double A1Ion = 12.899    ;
    const double A2Ion = 61.897    ;
    const double A3Ion =  9.27e3   ;
    const double A4Ion =  4.9749e-4;
    const double A5Ion =  3.989e-2 ;
    const double A6Ion = -1.59     ;
    const double A7Ion =  3.1834   ;
    const double A8Ion = -3.7154   ;
} CS_Ion;

struct Fusion_cross_section{
    const double A1Fusion =  47.88     ;
    const double A2Fusion = 482        ;
    const double A3Fusion =   3.08e-4  ;
    const double A4Fusion =   1.177    ;
    const double A5Fusion =   0        ;
} CS_Fusion;



#endif // CROSSSECTIONSCONT_HPP_INCLUDED
