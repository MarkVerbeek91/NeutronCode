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
    double A1cx;
    double A2cx;
    double A3cx;
    double A4cx;
    double A5cx;
    double A6cx;
} CS_cx;

struct Ion_impact_cross_sections{
    double A1Ion;
    double A2Ion;
    double A3Ion;
    double A4Ion;
    double A5Ion;
    double A6Ion;
    double A7Ion;
    double A8Ion;
} CS_Ion;

struct Fusion_cross_section{
    double A1Fusion;
    double A2Fusion;
    double A3Fusion;
    double A4Fusion;
    double A5Fusion;
} CS_Fusion;



#endif // CROSSSECTIONSCONT_HPP_INCLUDED
