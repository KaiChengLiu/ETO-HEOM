#include "constant.h"

gsl_complex ONE = gsl_complex_rect(1.0, 0.0);
gsl_complex ZERO = gsl_complex_rect(0.0, 0.0);
gsl_complex MINUS_ONE = gsl_complex_rect(-1.0, 0.0);
gsl_complex MINUS_ICNT = gsl_complex_rect(0.0, -1.0);
gsl_complex ICNT = gsl_complex_rect(0.0, 1.0);
double boltz_k = 0.695; // cm^-1 / K
double h_bar = 5.308e3; // cm^-1 * fs

double pi = 3.141592653589793238462643383279502884L;