#include "HEOM_constant.h"

const int ThreadPerBlock = 16;
const float pi = 3.1415926535897932;
const data_type ONE = make_cuComplex(1.0, 0.0);
const data_type ZERO = make_cuComplex(0.0, 0.0);
const data_type MINUS_ONE = make_cuComplex(-1.0, 0.0);
const data_type MINUS_ICNT = make_cuComplex(0.0, -1.0);
const data_type ICNT = make_cuComplex(0.0, 1.0);
const data_type boltz_k = make_cuComplex(0.695, 0.0); // cm^-1 / K
const data_type h_bar = make_cuComplex(5.308e+3, 0.0); // cm^-1 * fs