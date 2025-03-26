#ifndef HEOM_CONSTANT_H
#define HEOM_CONSTANT_H
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

extern gsl_complex ONE;
extern gsl_complex ZERO;
extern gsl_complex MINUS_ONE;
extern gsl_complex MINUS_ICNT;
extern gsl_complex ICNT;
extern double boltz_k; // cm^-1 / K
extern double h_bar; // cm^-1 * fs

extern double pi;
#endif