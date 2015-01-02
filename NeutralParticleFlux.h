#ifndef NEUTRALPARTICLEFLUX_H_INCLUDED
#define NEUTRALPARTICLEFLUX_H_INCLUDED

// neutrals from Class I ions

double Sfn1_InMinInte(double r, double dr );
double Sfn1_InMin( double r );
double Sfn1_OutMinInte(double r, double dr );
double Sfn1_OutMin( double r );
double Sfn1_InPlusInte(double r, double dr);
double Sfn1_InPlus(double r);
double Sfn1_OutPlusInte(double r, double dr);
double Sfn1_OutPlus(double r);

// neutrals from Class II ions
double Sfn2_InMinInte(double r, double dr );
double Sfn2_InMin( double r );
double Sfn2_OutMinInte(double r, double dr );
double Sfn2_OutMin( double r );
double Sfn2_InPlusInte(double r, double dr);
double Sfn2_InPlus(double r);
double Sfn2_OutPlusInte(double r, double dr);
double Sfn2_OutPlus(double r);

#endif // NEUTRALPARTICLEFLUX_H_INCLUDED
