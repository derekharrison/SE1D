/*
 * user_types.hpp
 *
 *  Created on: Dec 20, 2020
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

#include "complex.hpp"

typedef struct domain_data {
    int Nx;
    int Nt;
    double L;
} d_data;

typedef struct time_data {
    double to;
    double tf;
} t_data;

typedef struct physical_params {
    double h;
    double m;
} p_params;

typedef struct boundaries_and_source {
    Complex q;
    Complex psi_west;
    Complex psi_east;
} bound_and_source;

typedef struct solver_data {
    Complex* psi;
    double* x_c;
    double error_real;
    double error_im;
} s_data;

#endif /* USER_TYPES_HPP_ */
