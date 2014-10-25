
// yeay, i know, this is most likely not the correct way to split up code but for now I'm lazy to read the correct way
// i will regret that

double interpolation(double r)
{
    double value = -1;
    int low = N_TABLE - 2, upper = N_TABLE - 1;

    for ( int i = 0; i < N_TABLE; i++)
    {
        if ( r < Table.R[i])
        {
            low = i - 1;
            upper = i;
        }
    }

    double slope = ( Table.S_4[upper] - Table.S_4[low])/(Table.R[upper] - Table.R[low]);

    value = slope*r + Table.S_4[low] - Table.R[low] * slope;

    if (value < 0)
        value *= -1;

    return value;
}

double NIntegration( double (*funcPtr)(double), double Start, double End)
{
    double sum = (*funcPtr)(Start) + (*funcPtr)(End), step = (End - Start)/N_pres;

    for (double r=Start; r<End; r += step)
    {
        sum += 2.0 * (*funcPtr)(r);
    }

    sum = sum * step / 2.0;

    return sum;
}

double NIntegration_2( double (*funcPtr)(double, double), double Bar, double Start, double End)
{
    double sum = (*funcPtr)(Bar,Start) + (*funcPtr)(Bar,End), step = (End - Start)/N_pres;

    for (double r=Start; r<End; r += step)
    {
        sum += 2.0 * (*funcPtr)(Bar, r);
    }

    sum = sum * step / 2.0;

    return sum;
}

// special integrator for the Class II ions.
double NIntegration_3( double (*funcPtr)(double, double), double Start, double End)
{
    double sum = (*funcPtr)(Start,End) + (*funcPtr)(End,End), step = (End - Start)/N_pres;

    for (double r=Start; r<End; r += step)
    {
        sum += 2.0 * (*funcPtr)(r, End);
    }

    sum = sum * step / 2.0;

    return sum;
}

