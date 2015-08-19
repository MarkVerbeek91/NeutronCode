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

    // case i = 0
    Table.S[0] = Table.A[0];

    for ( int r = 0; r < N_TABLE; r++)
    {

    }
    for ( i = 1; i < N_TABLE; i++)
    {

        sum = 0;
        for ( j = 1; j < i - 1; j++)
            sum += Table.K[i][j] * Table.S[j];

        Table.S[i] = h * ( 0.5 * Table.K[i][0] * Table.S[0] + sum);
        Table.S[i] += Table.A[i];
        Table.S[i] /= ( 1 - (h/2.0) * Table.K[i][i] );
    }

    return;
}


