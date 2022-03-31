//
//  lib.hpp
//  two-bulb-imp
//
//  Created by dwh on 14/03/2022.
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>

#include "user_types.h"

void compute_bulb_compositions(e_params_t e_params,
                               p_params_t p_params,
                               t_params_t t_params,
                               int ng,
                               b_data_t & bulb_data);

void testing();

void compute_bulb_compositions_NNNN(e_params_t e_params,
                                    p_params_t_NNNN p_params,
                                    p_params_t p_params_ref,
                                    t_params_t t_params,
                                    int ng,
                                    int n,
                                    b_data_t_NNNN & bulb_data,
                                    b_data_t & bulb_data_ref);
#endif /* lib_hpp */
