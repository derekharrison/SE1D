/*
 * main.cpp
 *
 *  Created on: Dec 20, 2020
 *      Author: d-w-h
 *
 *      Numerical solution of the Schrodinger equation:
 *
 *      ih*d(psi)/dt = -h*h/2m*d^2(psi)/dx^2 + q
 */

#include <math.h>
#include <stdio.h>
#include "complex.hpp"
#include "solver.hpp"
#include "user_types.hpp"

int main(int argc, char* argv[]) {
    d_data domain_data;
    t_data time_data;
    p_params physical_params;
    bound_and_source boundaries_and_source;
    s_data solver_data;

    /* Parameters */
    domain_data.Nx = 10;                           //Number of nodes along x axis
    domain_data.Nt = 20;                           //Number of timesteps
    domain_data.L = 1.0;                           //Length of domain

    time_data.to = 0.0;                            //Initial time
    time_data.tf = 1.5;                            //Final time

    physical_params.h = 1.0;                       //Constant
    physical_params.m = 1.0;                       //Mass of particle

    boundaries_and_source.q.a = 0.5;               //Real part of source term in SE equation
    boundaries_and_source.q.b = 0.5;               //Imaginary part of source term in SE equation
    boundaries_and_source.psi_west.a = 0.0;        //Real part of west boundary
    boundaries_and_source.psi_west.b = 0.0;        //Imaginary part of west boundary
    boundaries_and_source.psi_east.a = 0.0;        //Real part of east boundary
    boundaries_and_source.psi_east.b = 0.0;        //Imaginary part of east boundary

    /* Allocate memory for solver results */
    solver_data.psi = new Complex[domain_data.Nx];
    solver_data.x_c = new double[domain_data.Nx];

    /* Execute solver */
    solver(domain_data,
           time_data,
           physical_params,
           boundaries_and_source,
           &solver_data);

    /* Print results */
    for(int j = 0; j < domain_data.Nx; ++j) {
        printf("x_c[%i]: %f, real(psi[%i]): %f, Im(psi[%i]): %f\n",
                j, solver_data.x_c[j], j, solver_data.psi[j].a, j, solver_data.psi[j].b);
    }

    printf("error_real: %E\n", solver_data.error_real);
    printf("error_im: %E\n", solver_data.error_im);

    /* Deallocate memory for solver results */
    delete [] solver_data.psi;
    delete [] solver_data.x_c;

    return 0;
}
