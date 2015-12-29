/**
    The calculation of the source rate of the source rate for the
    class II ions.
*/


#include "constants.h"
#include "SourceFunction.h"

void S(void)
{
    // a volterra solver
    int i, j;
    double h, sum;

    h = ( giveAnodeRadius() - giveCathodeRadius()) / N_TABLE;

    // case i = N
    Table->S[N_TABLE-1] = Table->A[N_TABLE-1];

    for ( i = N_TABLE-2; i >= 0; i--)
    {
        sum = 0;
        for ( j = N_TABLE-2; j > i; j--)
            sum += Table->K[i][j] * Table->S[j];

        Table->S[i] = h * ( 0.5 * Table->K[0][i] * Table->S[N_TABLE-1] + sum);
        Table->S[i] += Table->A[i];
        Table->S[i] /= ( 1 -  0.5 * h * Table->K[i][i] );
    }

    return;
}


