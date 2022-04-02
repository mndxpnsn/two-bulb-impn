//
//  main.cpp
//  two-bulb-imp
//
//  Created by dwh on 12/03/2022.
//  Three-component two-bulb diffusion
//  using implicit time discretization
//

#include <iostream>

#include "lib.hpp"
#include "user_types.h"

int main(int argc, const char * argv[]) {
    
    // Number of grid nodes
    int ng = 5;
    int n = 4;
    
    // Experimental setup parameters
    e_params_t e_params;
    e_params.V = 5e-4; // Volume of compartments (m3)
    e_params.d = 2e-3; // Diameter of tube connecting compartments (m)
    e_params.len = 1e-2; // Length of tube connecting compartments (m)
    e_params.A = 0.25 * 3.14 * e_params.d * e_params.d; // Cross-section area of tube (m2)
    e_params.dz = e_params.len / ng; // Tube resolution (m)
    
    // Initial composition bulb 1
    b_data_t bulb_data;
    bulb_data.mol_fracs_bulb1.x1 = 0.0;//0.501; // Bulb 1 H2 fraction
    bulb_data.mol_fracs_bulb1.x2 = 0.501;//0.499; // Bulb 1 N2 fraction
    bulb_data.mol_fracs_bulb1.x3 = 1.0 - 0.501 - 0.0; // Bulb 1 CO2 fraction
    
    // Initial composition bulb 2
    bulb_data.mol_fracs_bulb2.x1 = 0.501;//0.0; // Bulb 2 H2 fraction
    bulb_data.mol_fracs_bulb2.x2 = 0.499; // Bulb 2 N2 fraction
    bulb_data.mol_fracs_bulb2.x3 = 1.0 - 0.499 - 0.501; // Bulb 2 CO2 fraction
    
    // Initial composition bulb 1
    b_data_t_NNNN bulb_data_N;
    bulb_data_N.mol_fracs_bulb1.n = n;
    bulb_data_N.mol_fracs_bulb1.x = new double[n];
    bulb_data_N.mol_fracs_bulb1.x[0] = 0.399;//0.201;
    bulb_data_N.mol_fracs_bulb1.x[1] = 0.2;//0.0;
    bulb_data_N.mol_fracs_bulb1.x[2] = 0.25;//0.15;
    bulb_data_N.mol_fracs_bulb1.x[3] = 1.0 - 0.399 - 0.2 - 0.25;
    
    // Initial composition bulb 2
    bulb_data_N.mol_fracs_bulb2.n = n;
    bulb_data_N.mol_fracs_bulb2.x = new double[n];
    bulb_data_N.mol_fracs_bulb2.x[0] = 0.201;//0.399;
    bulb_data_N.mol_fracs_bulb2.x[1] = 0.0;//0.2;
    bulb_data_N.mol_fracs_bulb2.x[2] = 0.15;//0.25;
    bulb_data_N.mol_fracs_bulb2.x[3] = 1.0 - 0.201 - 0.15;
    
    // Total concentration
    p_params_t p_params;
    p_params.ct = 1.0; // Total concentration (mol / m3)
 
    // Total concentration
    p_params_t_NNNN p_params_N;
    p_params_N.ct = 1.0; // Total concentration (mol / m3)
    
    // Time parameters
    t_params_t t_params;
    t_params.to = 0.0; // Initial time (h)
    t_params.tf = 10.0; // Final time (h)
    t_params.nt = 10; // Number of time steps.
    t_params.dt = (t_params.tf - t_params.to) / t_params.nt; // Time sampling
    
    // Diffusivities
    p_params.D12 = 8.33e-5 * 3600; // units (m2 / h)
    p_params.D13 = 6.8e-5 * 3600; // units (m2 / h)
    p_params.D23 = 1.68e-5 * 3600; // units (m2 / h)

    // Diffusivities NNNN
    p_params_N.D = new double * [n];
    for(int j = 0; j < n; ++j) {
        p_params_N.D[j] = new double[n];
        for(int c = 0; c < n; ++c) {
            p_params_N.D[j][c] = 0.0;
        }
    }
    
    double D12 = 8.33e-5 * 3600; // units are (m2 / h)
    double D13 = 6.8e-5 * 3600; // units are (m2 / h)
    double D14 = 3.8e-5 * 3600; // units are (m2 / h)
    double D23 = 1.68e-5 * 3600; // units are (m2 / h)
    double D24 = 4.68e-5 * 3600; // units are (m2 / h)
    double D34 = 5.68e-5 * 3600; // units are (m2 / h)
    
    p_params_N.D[0][1] = D12; // units (m2 / h)
    p_params_N.D[0][2] = D13; // units (m2 / h)
    p_params_N.D[0][3] = D14; // units (m2 / h)
    p_params_N.D[1][2] = D23; // units (m2 / h)
    p_params_N.D[1][3] = D24; // units (m2 / h)
    p_params_N.D[2][3] = D34; // units (m2 / h)
    
    p_params_N.D[1][0] = D12; // units (m2 / h)
    p_params_N.D[2][0] = D13; // units (m2 / h)
    p_params_N.D[3][0] = D14; // units (m2 / h)
    p_params_N.D[2][1] = D23; // units (m2 / h)
    p_params_N.D[3][1] = D24; // units (m2 / h)
    p_params_N.D[3][2] = D34; // units (m2 / h)
    
    p_params_N.D12 = 8.33e-5 * 3600; // units (m2 / h)
    p_params_N.D13 = 6.8e-5 * 3600; // units (m2 / h)
    p_params_N.D23 = 1.68e-5 * 3600; // units (m2 / h)
    
    // Perform two-bulb diffusion experiment
//    compute_bulb_compositions(e_params, p_params, t_params, ng, bulb_data);
    
    // Perform two-bulb diffusion experiment
    compute_bulb_compositions_NNNN(e_params, p_params_N, p_params, t_params, ng, n, bulb_data_N, bulb_data);

    // Print results
    std::cout << "ref data: " << std::endl;
    std::cout << "bulb 1 frac 1: " << bulb_data.mol_fracs_bulb1.x1 << std::endl;
    std::cout << "bulb 1 frac 2: " << bulb_data.mol_fracs_bulb1.x2 << std::endl;
    std::cout << "bulb 1 frac 3: " << bulb_data.mol_fracs_bulb1.x3 << std::endl;
    
    std::cout << "bulb 2 frac 1: " << bulb_data.mol_fracs_bulb2.x1 << std::endl;
    std::cout << "bulb 2 frac 2: " << bulb_data.mol_fracs_bulb2.x2 << std::endl;
    std::cout << "bulb 2 frac 3: " << bulb_data.mol_fracs_bulb2.x3 << std::endl;
    
    
    // Print results
    std::cout << "calc data: " << std::endl;
    std::cout << "bulb 1 frac 1: " << bulb_data_N.mol_fracs_bulb1.x[0] << std::endl;
    std::cout << "bulb 1 frac 2: " << bulb_data_N.mol_fracs_bulb1.x[1] << std::endl;
    std::cout << "bulb 1 frac 3: " << bulb_data_N.mol_fracs_bulb1.x[2] << std::endl;
    std::cout << "bulb 1 frac 4: " << bulb_data_N.mol_fracs_bulb1.x[3] << std::endl;
    
    std::cout << "bulb 2 frac 1: " << bulb_data_N.mol_fracs_bulb2.x[0] << std::endl;
    std::cout << "bulb 2 frac 2: " << bulb_data_N.mol_fracs_bulb2.x[1] << std::endl;
    std::cout << "bulb 2 frac 3: " << bulb_data_N.mol_fracs_bulb2.x[2] << std::endl;
    std::cout << "bulb 2 frac 4: " << bulb_data_N.mol_fracs_bulb2.x[3] << std::endl;
    
//    testing();
//    
    return 0;
}
