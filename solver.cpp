/*
 * solver.cpp
 *
 *  Created on: Dec 20, 2020
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include "complex.hpp"
#include "user_types.hpp"

void solver(d_data domain_data,
            t_data time_data,
            p_params physical_params,
            bound_and_source boundaries_and_source,
            s_data* solver_data) {

    int Nx, Nt;
    double to, t, tf, m, h, L;
    Complex q(0,0);
    Complex psi_west(0,0);
    Complex psi_east(0,0);

    /* Parameters */
    Nx = domain_data.Nx;
    Nt = domain_data.Nt;
    L = domain_data.L;

    to = time_data.to;
    tf = time_data.tf;

    h = physical_params.h;
    m = physical_params.m;

    q.a = boundaries_and_source.q.a;
    q.b = boundaries_and_source.q.b;
    psi_west.a = boundaries_and_source.psi_west.a;
    psi_west.b = boundaries_and_source.psi_west.b;
    psi_east.a = boundaries_and_source.psi_east.a;
    psi_east.b = boundaries_and_source.psi_east.b;

    /* Allocate memory for psi data */
    Complex* psip = new Complex[Nx];
    Complex* psi_prev = new Complex[Nx];
    Complex* psipo = new Complex[Nx];

    /* Start calculations */
    Complex i(0,0);
    int max_iter = 300;
    double d_x = L/Nx;
    double d_t = (tf - to)/Nt;
    t = 0.0;

    /* Initialize psi and x_c data */
    for(int j = 0; j < Nx; ++j) {
        psip[j].a = 0.0;
        psip[j].b = 0.0;
        psipo[j].a = 0.0;
        psipo[j].b = 0.0;
        solver_data->x_c[j] = j*d_x + 0.5*d_x;
    }

    while(t < tf){
        /* Start Gauss-Seidel iterations */
        int it = 0;
        while(it < max_iter) {
            /* Check values previous iteration */
            for(int j = 0; j < Nx; ++j) {
                psi_prev[j].a = psip[j].a;
                psi_prev[j].b = psip[j].b;
            }

            /* First node */
            Complex psi_w(0,0);
            psi_w.a = -h*h/(m*d_x)*psi_west.a;
            psi_w.b = -h*h/(m*d_x)*psi_west.b;

            Complex psi_e(0,0);
            psi_e.a = -h*h/(2*m*d_x)*psip[1].a;
            psi_e.b = -h*h/(2*m*d_x)*psip[1].b;

            Complex qdx(0,0);
            qdx.a = q.a*d_x;
            qdx.b = q.b*d_x;

            Complex ihdx_psipo_dt(0,0);
            ihdx_psipo_dt = i * psipo[0];
            ihdx_psipo_dt.a = ihdx_psipo_dt.a*h*d_x/d_t;
            ihdx_psipo_dt.b = ihdx_psipo_dt.b*h*d_x/d_t;

            Complex Fp = psi_w + psi_e + qdx + ihdx_psipo_dt;

            double a = 3/2*h*h/(m*d_x);
            double b = h*d_x/d_t;
            Complex ib_plus_a(a,b);

            Complex min_Fp(-Fp.a, -Fp.b);

            psip[0] = min_Fp * ib_plus_a;
            psip[0].a = psip[0].a/(a*a + b*b);
            psip[0].b = psip[0].b/(a*a + b*b);

            /* Central nodes */
            for(int j = 1; j < Nx - 1; ++j) {
                psi_w.a = -h*h/(2*m*d_x)*psip[j-1].a;
                psi_w.b = -h*h/(2*m*d_x)*psip[j-1].b;

                psi_e.a = -h*h/(2*m*d_x)*psip[j+1].a;
                psi_e.b = -h*h/(2*m*d_x)*psip[j+1].b;

                ihdx_psipo_dt = i * psipo[j];
                ihdx_psipo_dt.a = ihdx_psipo_dt.a*h*d_x/d_t;
                ihdx_psipo_dt.b = ihdx_psipo_dt.b*h*d_x/d_t;

                Fp = psi_w + psi_e + qdx + ihdx_psipo_dt;

                a = h*h/(m*d_x);
                b = h*d_x/d_t;
                ib_plus_a.a = a;
                ib_plus_a.b = b;

                min_Fp.a = Fp.a*-1;
                min_Fp.b = Fp.b*-1;

                psip[j] = min_Fp * ib_plus_a;
                psip[j].a = psip[j].a/(a*a + b*b);
                psip[j].b = psip[j].b/(a*a + b*b);
            }

            /* Final node */
            psi_w.a = -h*h/(2*m*d_x)*psip[Nx-1-1].a;
            psi_w.b = -h*h/(2*m*d_x)*psip[Nx-1-1].b;

            psi_e.a = -h*h/(m*d_x)*psi_east.a;
            psi_e.b = -h*h/(m*d_x)*psi_east.b;

            ihdx_psipo_dt = i * psipo[Nx-1];
            ihdx_psipo_dt.a = ihdx_psipo_dt.a*h*d_x/d_t;
            ihdx_psipo_dt.b = ihdx_psipo_dt.b*h*d_x/d_t;

            Fp = psi_w + psi_e + qdx + ihdx_psipo_dt;

            a = 3/2*h*h/(m*d_x);
            b = h*d_x/d_t;
            ib_plus_a.a = a;
            ib_plus_a.b = b;

            min_Fp.a = Fp.a*-1;
            min_Fp.b = Fp.b*-1;

            psip[Nx-1] = min_Fp * ib_plus_a;
            psip[Nx-1].a = psip[Nx-1].a/(a*a + b*b);
            psip[Nx-1].b = psip[Nx-1].b/(a*a + b*b);

            it++;
        }

        /* Check convergence */
        solver_data->error_real = 0.0;
        solver_data->error_im = 0.0;
        for(int j = 0; j < Nx; ++j) {
            solver_data->error_real += fabs(psip[j].a - psi_prev[j].a);
            solver_data->error_im += fabs(psip[j].b - psi_prev[j].b);
        }

        solver_data->error_real /= Nx;
        solver_data->error_im /= Nx;

        /* Update old timestep values */
        for(int j = 0; j < Nx; ++j) {
            psipo[j].a = psip[j].a;
            psipo[j].b = psip[j].b;
        }

        t = t + d_t;
    }

    /* Set solver results */
    for(int j = 0; j < Nx; ++j) {
        solver_data->psi[j].a = psip[j].a;
        solver_data->psi[j].b = psip[j].b;
    }

    /* Deallocate memory for psi data */
    delete [] psip;
    delete [] psi_prev;
    delete [] psipo;
}
