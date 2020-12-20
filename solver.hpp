/*
 * solver.hpp
 *
 *  Created on: Dec 20, 2020
 *      Author: d-w-h
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include "user_types.hpp"

void solver(d_data domain_data,
            t_data time_data,
            p_params physical_params,
            bound_and_source boundaries_and_source,
            s_data* solver_data);

#endif /* SOLVER_HPP_ */
