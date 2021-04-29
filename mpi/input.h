#include <stdio.h>
#include <stdlib.h> // for exit(0)

// PDE parameters

const double xmin = 0.0, xmax = 1.0;
const double coeff = 1.0; // Not working for different coefficient.

double init_condn(double x)
{
  double value = sin(2.0 * M_PI * x);
  return value;
}

const int N = 100;                                // number of points
const double cfl = 0.95;
const double Tf = 0.2;
