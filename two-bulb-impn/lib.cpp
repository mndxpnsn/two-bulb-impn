//
//  lib.cpp
//  two-bulb-imp
//
//  Created by dwh on 14/03/2022.
//

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "lib.hpp"
#include "user_types.h"

double a1(p_params_t & p_params, double x2) {
    double v1 = 1.0 / p_params.D13;
    double v2 = 1.0 / p_params.D12;
    
    return (v2 - v1) * x2 + v1;
}

double a2(p_params_t & p_params, double x1) {
    double v1 = 1.0 / p_params.D13;
    double v2 = 1.0 / p_params.D12;
    
    return x1 * (v1 - v2);
}

double a1_NNNN(p_params_t_NNNN & p_params, double x2) {
    double v1 = 1.0 / p_params.D13;
    double v2 = 1.0 / p_params.D12;
    
    return (v2 - v1) * x2 + v1;
}

double a2_NNNN(p_params_t_NNNN & p_params, double x1) {
    double v1 = 1.0 / p_params.D13;
    double v2 = 1.0 / p_params.D12;
    
    return x1 * (v1 - v2);
}

double b1(p_params_t & p_params, double x2) {
    double v2 = 1.0 / p_params.D12;
    double v3 = 1.0 / p_params.D23;
    
    return x2 * (v3 - v2);
}

double b2(p_params_t & p_params, double x1) {
    double v2 = 1.0 / p_params.D12;
    double v3 = 1.0 / p_params.D23;
    
    return (v2 - v3) * x1 + v3;
}

double b1_NNNN(p_params_t_NNNN & p_params, double x2) {
    double v2 = 1.0 / p_params.D12;
    double v3 = 1.0 / p_params.D23;
    
    return x2 * (v3 - v2);
}

double b2_NNNN(p_params_t_NNNN & p_params, double x1) {
    double v2 = 1.0 / p_params.D12;
    double v3 = 1.0 / p_params.D23;
    
    return (v2 - v3) * x1 + v3;
}

double alpha1(p_params_t & p_params, double x1, double x2) {
    double ct = p_params.ct;
    
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double b1_loc = b1(p_params, x2);
    double b2_loc = b2(p_params, x1);
    
    return -ct / (a2_loc - a1_loc * b2_loc / b1_loc);
}

double alpha2(p_params_t & p_params, double x1, double x2) {
    double ct = p_params.ct;
    
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double b1_loc = b1(p_params, x2);
    double b2_loc = b2(p_params, x1);
    
    return (a1_loc / b1_loc) * ct / (a2_loc - a1_loc * b2_loc / b1_loc);
}

double alpha1_NNNN(p_params_t_NNNN & p_params, double x1, double x2) {
    double ct = p_params.ct;
    
    double a1_loc = a1_NNNN(p_params, x2);
    double a2_loc = a2_NNNN(p_params, x1);
    double b1_loc = b1_NNNN(p_params, x2);
    double b2_loc = b2_NNNN(p_params, x1);
    
    return -ct / (a2_loc - a1_loc * b2_loc / b1_loc);
}

double alpha2_NNNN(p_params_t_NNNN & p_params, double x1, double x2) {
    double ct = p_params.ct;
    
    double a1_loc = a1_NNNN(p_params, x2);
    double a2_loc = a2_NNNN(p_params, x1);
    double b1_loc = b1_NNNN(p_params, x2);
    double b2_loc = b2_NNNN(p_params, x1);
    
    return (a1_loc / b1_loc) * ct / (a2_loc - a1_loc * b2_loc / b1_loc);
}

double beta1(p_params_t & p_params, double x1, double x2) {
    double ct = p_params.ct;
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double alpha1_loc = alpha1(p_params, x1, x2);
    
    return -ct / a1_loc - a2_loc * alpha1_loc / a1_loc;
}

double beta2(p_params_t & p_params, double x1, double x2) {
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double alpha2_loc = alpha2(p_params, x1, x2);
    
    return - a2_loc * alpha2_loc / a1_loc;
}

double beta1_NNNN(p_params_t_NNNN & p_params, double x1, double x2) {
    double ct = p_params.ct;
    double a1_loc = a1_NNNN(p_params, x2);
    double a2_loc = a2_NNNN(p_params, x1);
    double alpha1_loc = alpha1_NNNN(p_params, x1, x2);
    
    return -ct / a1_loc - a2_loc * alpha1_loc / a1_loc;
}

double beta2_NNNN(p_params_t_NNNN & p_params, double x1, double x2) {
    double a1_loc = a1_NNNN(p_params, x2);
    double a2_loc = a2_NNNN(p_params, x1);
    double alpha2_loc = alpha2_NNNN(p_params, x1, x2);
    
    return - a2_loc * alpha2_loc / a1_loc;
}

void bulb1_c1(c_data_t & comp_data) {

    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb1.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb1.x2;

    double beta1_loc = beta1(comp_data.p_params, x1, x2);
    double beta2_loc = beta2(comp_data.p_params, x1, x2);
    
    // Store coefficients
    comp_data.coeff_mat_boundary[0].A_inv[0][0] = beta1_loc;
    comp_data.coeff_mat_boundary[0].A_inv[0][1] = beta2_loc;

    double ap1 = V * comp_data.p_params.ct / dt - 2 * beta1_loc * A / dz;

    double x11 = comp_data.tube_fracs[0].x1;
    double x21 = comp_data.tube_fracs[0].x2;

    double old_term1 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb1.x1 / dt;
    
    double x1_b1 = -2 * beta1_loc * x11 * A / dz - 2 * beta2_loc * (x21 - x2) * A / dz + old_term1;

    x1_b1 = x1_b1 / ap1;

    comp_data.bulb_data.mol_fracs_bulb1.x1 = x1_b1;
    
//    std::cout << "A_inv_b1 ref:" << std::endl;
//    std::cout << beta1_loc << ", " << beta2_loc << std::endl;
}

void bulb1_c2(c_data_t & comp_data) {

    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb1.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb1.x2;

    double alpha1_loc = alpha1(comp_data.p_params, x1, x2);
    double alpha2_loc = alpha2(comp_data.p_params, x1, x2);

    // Store coefficients
    comp_data.coeff_mat_boundary[0].A_inv[1][0] = alpha1_loc;
    comp_data.coeff_mat_boundary[0].A_inv[1][1] = alpha2_loc;
    
    double ap2 = V * comp_data.p_params.ct / dt - 2 * alpha2_loc * A / dz;

    double x11 = comp_data.tube_fracs[0].x1;
    double x21 = comp_data.tube_fracs[0].x2;

    double old_term2 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb1.x2 / dt;

    double x2_b1 = -2 * alpha1_loc / dz * (x11 - x1) * A  - 2 * alpha2_loc * x21 * A / dz + old_term2;

    x2_b1 = x2_b1 / ap2;

    comp_data.bulb_data.mol_fracs_bulb1.x2 = x2_b1;
    comp_data.bulb_data.mol_fracs_bulb1.x3 = 1.0 - comp_data.bulb_data.mol_fracs_bulb1.x1 - x2_b1;
    
//    std::cout << alpha1_loc << ", " << alpha2_loc << std::endl;
}

void bulb2_c1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb2.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb2.x2;

    double beta1_loc = beta1(comp_data.p_params, x1, x2);
    double beta2_loc = beta2(comp_data.p_params, x1, x2);

    // Store coefficients
    comp_data.coeff_mat_boundary[ng].A_inv[0][0] = beta1_loc;
    comp_data.coeff_mat_boundary[ng].A_inv[0][1] = beta2_loc;
    
    double ap1 = V * comp_data.p_params.ct / dt - 2 * beta1_loc * A / dz;

    double x1n = comp_data.tube_fracs[ng - 1].x1;
    double x2n = comp_data.tube_fracs[ng - 1].x2;

    double old_term1 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb2.x1 / dt;
    
    double x1_b2 = -2 * beta1_loc * x1n * A / dz + 2 * beta2_loc * (x2 - x2n) * A / dz + old_term1;

    x1_b2 = x1_b2 / ap1;

    comp_data.bulb_data.mol_fracs_bulb2.x1 = x1_b2;
    
//    std::cout << "A_inv_b2 ref:" << std::endl;
//    std::cout << beta1_loc << ", " << beta2_loc << std::endl;
}

void bulb2_c2(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb2.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb2.x2;

    double alpha1_loc = alpha1(comp_data.p_params, x1, x2);
    double alpha2_loc = alpha2(comp_data.p_params, x1, x2);

    // Store coefficients
    comp_data.coeff_mat_boundary[ng].A_inv[1][0] = alpha1_loc;
    comp_data.coeff_mat_boundary[ng].A_inv[1][1] = alpha2_loc;
    
    double ap2 = V * comp_data.p_params.ct / dt - 2 * alpha2_loc * A / dz;

    double x1n = comp_data.tube_fracs[ng - 1].x1;
    double x2n = comp_data.tube_fracs[ng - 1].x2;

    double old_term2 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb2.x2 / dt;

    double x2_b2 = 2 * alpha1_loc / dz * (x1 - x1n) * A  - 2 * alpha2_loc * x2n * A / dz + old_term2;

    x2_b2 = x2_b2 / ap2;

    comp_data.bulb_data.mol_fracs_bulb2.x2 = x2_b2;
    comp_data.bulb_data.mol_fracs_bulb2.x3 = 1.0 - comp_data.bulb_data.mol_fracs_bulb2.x1 - x2_b2;
    
//    std::cout << alpha1_loc << ", " << alpha2_loc << std::endl;
}

void tube0_c1(c_data_t & comp_data) {

    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = comp_data.bulb_data_inter.mol_fracs_bulb1.x1;
    double x2w = comp_data.bulb_data_inter.mol_fracs_bulb1.x2;
    double x1e = 0.5 * (comp_data.tube_fracs_inter[0].x1 + comp_data.tube_fracs_inter[1].x1);
    double x2e = 0.5 * (comp_data.tube_fracs_inter[0].x2 + comp_data.tube_fracs_inter[1].x2);
    
    double x2P = comp_data.tube_fracs[0].x2;
    double x1E = comp_data.tube_fracs[1].x1;
    double x2E = comp_data.tube_fracs[1].x2;
    
    double beta1w = beta1(comp_data.p_params, x1w, x2w);
    double beta2w = beta2(comp_data.p_params, x1w, x2w);
    double beta1e = beta1(comp_data.p_params, x1e, x2e);
    double beta2e = beta2(comp_data.p_params, x1e, x2e);
    
    // Store coefficients
    comp_data.coeff_mat_boundary[0].A_inv[0][0] = beta1w;
    comp_data.coeff_mat_boundary[0].A_inv[0][1] = beta2w;
    comp_data.coeff_mat_boundary[1].A_inv[0][0] = beta1e;
    comp_data.coeff_mat_boundary[1].A_inv[0][1] = beta2e;
    
    double ap1_tube0 = comp_data.p_params.ct * dz / dt - beta1w / (0.5 * dz) - beta1e / dz;
    
    double old_term_tube1 = comp_data.p_params.ct * dz * comp_data.tube_fracs_old[0].x1 / dt;
    
    double x1_tube0 = - beta1w / (0.5 * dz) * x1w + beta2w / (0.5 * dz) * (x2P - x2w) - beta1e / dz * x1E - beta2e / dz * (x2E - x2P) + old_term_tube1;
    
    x1_tube0 = x1_tube0 / ap1_tube0;
    
    comp_data.tube_fracs[0].x1 = x1_tube0;
}

void tube0_c2(c_data_t & comp_data) {
    
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = comp_data.bulb_data_inter.mol_fracs_bulb1.x1;
    double x2w = comp_data.bulb_data_inter.mol_fracs_bulb1.x2;
    double x1e = 0.5 * (comp_data.tube_fracs_inter[0].x1 + comp_data.tube_fracs_inter[1].x1);
    double x2e = 0.5 * (comp_data.tube_fracs_inter[0].x2 + comp_data.tube_fracs_inter[1].x2);
    
    double x1E = comp_data.tube_fracs[1].x1;
    double x2E = comp_data.tube_fracs[1].x2;
    double x1P = comp_data.tube_fracs[0].x1;
    
    double alpha1w = alpha1(comp_data.p_params, x1w, x2w);
    double alpha2w = alpha2(comp_data.p_params, x1w, x2w);
    double alpha1e = alpha1(comp_data.p_params, x1e, x2e);
    double alpha2e = alpha2(comp_data.p_params, x1e, x2e);
    
    // Store coefficients
    comp_data.coeff_mat_boundary[0].A_inv[1][0] = alpha1w;
    comp_data.coeff_mat_boundary[0].A_inv[1][1] = alpha2w;
    comp_data.coeff_mat_boundary[1].A_inv[1][0] = alpha1e;
    comp_data.coeff_mat_boundary[1].A_inv[1][1] = alpha2e;
    
    double ap2_tube0 = comp_data.p_params.ct * dz / dt - alpha2w / (0.5 * dz) - alpha2e / dz;
    
    double old_term_tube2 = comp_data.p_params.ct * dz * comp_data.tube_fracs_old[0].x2 / dt;
    
    double x2_tube0 = alpha1w / (0.5 * dz) * (x1P - x1w) - alpha2w / (0.5 * dz) * x2w - alpha1e / dz * (x1E - x1P) - alpha2e / dz * x2E + old_term_tube2;
    
    x2_tube0 = x2_tube0 / ap2_tube0;
    
    comp_data.tube_fracs[0].x2 = x2_tube0;
    comp_data.tube_fracs[0].x3 = 1.0 - comp_data.tube_fracs[0].x1 - x2_tube0;
}

void mid_nodes1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    for(int node = 1; node < ng - 1; ++node) {
        double x1w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x1 + comp_data.tube_fracs_inter[node].x1);
        double x2w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x2 + comp_data.tube_fracs_inter[node].x2);
        double x1e = 0.5 * (comp_data.tube_fracs_inter[node].x1 + comp_data.tube_fracs_inter[node + 1].x1);
        double x2e = 0.5 * (comp_data.tube_fracs_inter[node].x2 + comp_data.tube_fracs_inter[node + 1].x2);

        double x1W = comp_data.tube_fracs[node - 1].x1;
        double x2W = comp_data.tube_fracs[node - 1].x2;
        double x1E = comp_data.tube_fracs[node + 1].x1;
        double x2E = comp_data.tube_fracs[node + 1].x2;
        double x2P = comp_data.tube_fracs[node].x2;

        double beta1w = beta1(comp_data.p_params, x1w, x2w);
        double beta2w = beta2(comp_data.p_params, x1w, x2w);
        double beta1e = beta1(comp_data.p_params, x1e, x2e);
        double beta2e = beta2(comp_data.p_params, x1e, x2e);

        // Store coefficients
        comp_data.coeff_mat_boundary[node].A_inv[0][0] = beta1w;
        comp_data.coeff_mat_boundary[node].A_inv[0][1] = beta2w;
        comp_data.coeff_mat_boundary[node + 1].A_inv[0][0] = beta1e;
        comp_data.coeff_mat_boundary[node + 1].A_inv[0][1] = beta2e;
        
        double ap1_tube_node = comp_data.p_params.ct * dz / dt - beta1w / dz - beta1e / dz;
        
        double old_term_tube1 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[node].x1;

        double x1_tube_node = - beta1w / dz * x1W + beta2w / dz * (x2P - x2W) - beta1e / dz * x1E - beta2e / dz * (x2E - x2P) + old_term_tube1;

        x1_tube_node = x1_tube_node / ap1_tube_node;

        comp_data.tube_fracs[node].x1 = x1_tube_node;
    }
}

void mid_nodes2(c_data_t & comp_data)  {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    for(int node = 1; node < ng - 1; ++node) {
        double x1w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x1 + comp_data.tube_fracs_inter[node].x1);
        double x2w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x2 + comp_data.tube_fracs_inter[node].x2);
        double x1e = 0.5 * (comp_data.tube_fracs_inter[node].x1 + comp_data.tube_fracs_inter[node + 1].x1);
        double x2e = 0.5 * (comp_data.tube_fracs_inter[node].x2 + comp_data.tube_fracs_inter[node + 1].x2);

        double x1W = comp_data.tube_fracs[node - 1].x1;
        double x2W = comp_data.tube_fracs[node - 1].x2;
        double x1E = comp_data.tube_fracs[node + 1].x1;
        double x2E = comp_data.tube_fracs[node + 1].x2;
        double x1P = comp_data.tube_fracs[node].x1;

        double alpha1w = alpha1(comp_data.p_params, x1w, x2w);
        double alpha2w = alpha2(comp_data.p_params, x1w, x2w);
        double alpha1e = alpha1(comp_data.p_params, x1e, x2e);
        double alpha2e = alpha2(comp_data.p_params, x1e, x2e);

        // Store coefficients
        comp_data.coeff_mat_boundary[node].A_inv[1][0] = alpha1w;
        comp_data.coeff_mat_boundary[node].A_inv[1][1] = alpha2w;
        comp_data.coeff_mat_boundary[node + 1].A_inv[1][0] = alpha1e;
        comp_data.coeff_mat_boundary[node + 1].A_inv[1][1] = alpha2e;
        
        double ap2_tube_node = comp_data.p_params.ct * dz / dt - alpha2w / dz - alpha2e / dz;
        
        double old_term_tube2 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[node].x2;

        double x2_tube_node = alpha1w / dz * (x1P - x1W) - alpha2w / dz * x2W - alpha1e * (x1E - x1P) / dz - alpha2e / dz * x2E + old_term_tube2;

        x2_tube_node = x2_tube_node / ap2_tube_node;

        comp_data.tube_fracs[node].x2 = x2_tube_node;
    }
}

void tube_n_c1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x1 + comp_data.tube_fracs_inter[ng - 1].x1);
    double x2w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x2 + comp_data.tube_fracs_inter[ng - 1].x2);
    double x1e = comp_data.bulb_data_inter.mol_fracs_bulb2.x1;
    double x2e = comp_data.bulb_data_inter.mol_fracs_bulb2.x2;
    
    double x1W = comp_data.tube_fracs[ng - 2].x1;
    double x2W = comp_data.tube_fracs[ng - 2].x2;
    double x1E = x1e;
    double x2E = x2e;
    double x2P = comp_data.tube_fracs[ng - 1].x2;
    
    double beta1w = beta1(comp_data.p_params, x1w, x2w);
    double beta2w = beta2(comp_data.p_params, x1w, x2w);
    double beta1e = beta1(comp_data.p_params, x1e, x2e);
    double beta2e = beta2(comp_data.p_params, x1e, x2e);
    
    // Store coefficients
    comp_data.coeff_mat_boundary[ng - 1].A_inv[0][0] = beta1w;
    comp_data.coeff_mat_boundary[ng - 1].A_inv[0][1] = beta2w;
    comp_data.coeff_mat_boundary[ng].A_inv[0][0] = beta1e;
    comp_data.coeff_mat_boundary[ng].A_inv[0][1] = beta2e;
    
    double ap1_tube_n = comp_data.p_params.ct * dz / dt - beta1w / dz - beta1e / (0.5 * dz);
    
    double old_term_tube_n_1 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[ng - 1].x1;
    
    double x1_tube_n = - beta1w / dz * x1W + beta2w / dz * (x2P - x2W) - beta1e / (0.5 * dz) * x1E - beta2e / dz * (x2E - x2P) + old_term_tube_n_1;
    
    x1_tube_n = x1_tube_n / ap1_tube_n;
    
    comp_data.tube_fracs[ng - 1].x1 = x1_tube_n;
}

void tube_n_c2(c_data_t & comp_data)  {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x1 + comp_data.tube_fracs_inter[ng - 1].x1);
    double x2w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x2 + comp_data.tube_fracs_inter[ng - 1].x2);
    double x1e = comp_data.bulb_data_inter.mol_fracs_bulb2.x1;
    double x2e = comp_data.bulb_data_inter.mol_fracs_bulb2.x2;
    
    double x1W = comp_data.tube_fracs[ng - 2].x1;
    double x2W = comp_data.tube_fracs[ng - 2].x2;
    double x1E = x1e;
    double x2E = x2e;
    double x1P = comp_data.tube_fracs[ng - 1].x1;
    
    double alpha1w = alpha1(comp_data.p_params, x1w, x2w);
    double alpha2w = alpha2(comp_data.p_params, x1w, x2w);
    double alpha1e = alpha1(comp_data.p_params, x1e, x2e);
    double alpha2e = alpha2(comp_data.p_params, x1e, x2e);
    
    // Store coefficients
    comp_data.coeff_mat_boundary[ng - 1].A_inv[1][0] = alpha1w;
    comp_data.coeff_mat_boundary[ng - 1].A_inv[1][1] = alpha2w;
    comp_data.coeff_mat_boundary[ng].A_inv[1][0] = alpha1e;
    comp_data.coeff_mat_boundary[ng].A_inv[1][1] = alpha2e;
    
    double ap2_tube_n = comp_data.p_params.ct * dz / dt - alpha2w / dz - alpha2e / (0.5 * dz);
    
    double old_term_tube_n_2 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[ng - 1].x2;
    
    double x2_tube_n = alpha1w / dz * (x1P - x1W) - alpha2w / dz * x2W - alpha1e / (0.5 * dz) * (x1E - x1P) - alpha2e / (0.5 * dz) * x2E + old_term_tube_n_2;
    
    x2_tube_n = x2_tube_n / ap2_tube_n;
    
    comp_data.tube_fracs[ng - 1].x2 = x2_tube_n;
}

void set_frac_comp3(c_data_t & comp_data) {
    int ng = comp_data.ng;
    
    for(int node = 0; node < ng; ++node) {
        comp_data.tube_fracs[node].x3 = 1.0 - comp_data.tube_fracs[node].x1 - comp_data.tube_fracs[node].x2;
    }
}

void set_frac_comp3_NNNN(c_data_t_NNNN & comp_data) {
    int ng = comp_data.ng;
    int n = comp_data.n;
    
    for(int node = 0; node < ng; ++node) {
        double sum_loc = 0.0;
        for(int c = 0; c < n - 1; ++c) {
            sum_loc = sum_loc + comp_data.tube_fracs[node].x[c];
        }
        comp_data.tube_fracs[node].x[n - 1] = 1.0 - sum_loc;
    }
    
    double sum_loc1 = 0.0;
    for(int c = 0; c < n - 1; ++c) {
        sum_loc1 = sum_loc1 + comp_data.bulb_data.mol_fracs_bulb1.x[c];
    }
    
    double sum_loc2 = 0.0;
    for(int c = 0; c < n - 1; ++c) {
        sum_loc2 = sum_loc2 + comp_data.bulb_data.mol_fracs_bulb2.x[c];
    }
    
    comp_data.bulb_data.mol_fracs_bulb1.x[n - 1] = 1.0 - sum_loc1;
    comp_data.bulb_data.mol_fracs_bulb2.x[n - 1] = 1.0 - sum_loc2;
}

void update_composition_estimate(c_data_t & comp_data) {
    
    // Bulb 1, component 1
    bulb1_c1(comp_data);
    
    // Bulb 1, component 2
    bulb1_c2(comp_data);

    // Tube node 0, component 1
    tube0_c1(comp_data);
    
    // Tube node 0, component 2
    tube0_c2(comp_data);
    
    // Tube mid nodes, component 1
    mid_nodes1(comp_data);
    
    // Tube mid nodes, component 2
    mid_nodes2(comp_data);
    
    // Tube node n, component 1
    tube_n_c1(comp_data);
    
    // Tube node n, component 2
    tube_n_c2(comp_data);
    
    // Set mole fraction component 3
    set_frac_comp3(comp_data);
    
    // Bulb 2, component 1
    bulb2_c1(comp_data);
    
    // Bulb 2, component 2
    bulb2_c2(comp_data);
}

void compute_bulb_compositions(e_params_t e_params,
                               p_params_t p_params,
                               t_params_t t_params,
                               int ng,
                               b_data_t & bulb_data) {
    
    // Organize data
    c_data_t comp_data;
    comp_data.ng = ng;

    comp_data.p_params = p_params;
    comp_data.e_params = e_params;
    comp_data.t_params = t_params;

    comp_data.bulb_data = bulb_data;
    comp_data.bulb_data_old = bulb_data;
    comp_data.bulb_data_inter = bulb_data;
    
    // Allocate data for tube composition
    comp_data.tube_fracs = new node_t[ng];
    comp_data.tube_fracs_inter = new node_t[ng];
    comp_data.tube_fracs_old = new node_t[ng];
    
    // Initialize tube composition
    for(int node = 0; node < ng; ++node) {
        comp_data.tube_fracs[node].x1 = 1.0 / 3;
        comp_data.tube_fracs[node].x2 = 1.0 / 3;
        comp_data.tube_fracs[node].x3 = 1.0 / 3;

        comp_data.tube_fracs_inter[node].x1 = 1.0 / 3;
        comp_data.tube_fracs_inter[node].x2 = 1.0 / 3;
        comp_data.tube_fracs_inter[node].x3 = 1.0 / 3;

        comp_data.tube_fracs_old[node].x1 = 1.0 / 3;
        comp_data.tube_fracs_old[node].x2 = 1.0 / 3;
        comp_data.tube_fracs_old[node].x3 = 1.0 / 3;
    }
    
    // Compute composition
    int max_out_it = MAX_OUT;
    int max_in_it = MAX_IN;
    
    double t = t_params.to;
    double dt = (t_params.tf - t_params.to) / t_params.nt;
    
    // Loop to time t = tf
    while(t < t_params.tf) {
        
        // Update bulb composition of previous time step
        comp_data.bulb_data_old = comp_data.bulb_data;
        
        // Update tube composition of previous time step
        for(int node = 0; node < ng; ++node)
            comp_data.tube_fracs_old[node] = comp_data.tube_fracs[node];
        
        // Outer Gauss-Seidel iterations
        int out_it = 0;
        while(out_it < max_out_it) {
            
            // Update intermediate bulb composition
            comp_data.bulb_data_inter = comp_data.bulb_data;
            
            // Update intermediate tube composition
            for(int node = 0; node < ng; ++node)
                comp_data.tube_fracs_inter[node] = comp_data.tube_fracs[node];
            
            // Inner Gauss-Seidel iterations
            int in_it = 0;
            while(in_it < max_in_it) {

                // Update estimates bulb and tube composition
                update_composition_estimate(comp_data);
                
                in_it++;
            }
            
            // Print current tube composition
//            std::cout << "ref 3 comps" << std::endl;
//            for(int node = 0; node < ng; ++node) {
//                printf("%f, ", comp_data.tube_fracs[node].x1);
//                printf("%f, ", comp_data.tube_fracs[node].x2);
//                printf("%f, ", comp_data.tube_fracs[node].x3);
//                printf("\n");
//            }
            
            out_it++;
        }
        
        t = t + dt;
    }
    
    // Set bulb data
    bulb_data = comp_data.bulb_data;
    
    delete [] comp_data.tube_fracs;
    delete [] comp_data.tube_fracs_inter;
    delete [] comp_data.tube_fracs_old;
}

double Determinant(double ** a, int n) {
   int i, j, j1, j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) {
       exit(5);
   }
   else if (n == 1) {
       det = a[0][0];
   }
   else if (n == 2) {
       det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   }
   else {
       det = 0;
       for (j1 = 0; j1 < n; j1++)
       {
           m = new double * [n - 1];
           for (i = 0; i < n - 1; i++)
               m[i] = new double[n - 1];

           for (i = 1; i < n; i++)
           {
               j2 = 0;

               for (j = 0; j < n; j++)
               {
                   if (j == j1)
                       continue;

                   m[i-1][j2] = a[i][j];
                   j2++;
               }
           }

           det += pow(-1.0, j1 + 2.0) * a[0][j1] * Determinant(m, n - 1);

           for (i = 0; i < n - 1; i++)
               delete [] m[i];

           delete [] m;
      }
   }

   return det;

}


// Find the cofactor matrix of a square matrix
void CoFactor(double ** a, int n, double ** b)
{
    
    int i, j, ii, jj, i1, j1;
    double det;
    double ** c;

    c = new double * [n - 1];
    for (i = 0; i < n - 1; i++)
        c[i] = new double[n - 1];

    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {

            // Form the adjoint a_ij
            i1 = 0;

            for (ii = 0; ii < n; ii++) {
                if (ii == i)
                    continue;

                j1 = 0;

                for (jj = 0; jj < n; jj++) {
                    if (jj == j)
                        continue;

                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }

            // Calculate the determinate
            det = Determinant(c, n - 1);

            // Fill in the elements of the cofactor
            b[i][j] = pow(-1.0, i + j + 2.0) * det;
        }
    }

    for (i = 0; i < n - 1; i++)
        delete [] c[i];

    delete [] c;
}


// Transpose of a square matrix, do it in place
void Transpose(double ** a, int n) {
    int i, j;
    double tmp;

    for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}

void compute_mat_inv(double ** mat, int n, double ** mat_inv) {
    int i, j;
    
    // Allocate data for adj matrix
    double ** mat_adj = new double * [n];
    
    for(int i = 0; i < n; ++i) {
        mat_adj[i] = new double[n];
    }
    
    // Calculating the determinant of M
    double det = Determinant(mat, n);


    // If matrix is singular, exit
    if(!det) {
        printf("\nTransformation matrix is singular\n");
        exit(4);
    }


    // Calculating co-factor M_adj of M
    CoFactor(mat, n, mat_adj);


    // Transposing M_adj
    Transpose(mat_adj, n);


    // Calculating inverse of transformation matrix
    for(i = 0; i < n; ++i) {
        for(j = 0; j < n; ++j) {
            mat_inv[i][j] = 1.0 / det * mat_adj[i][j];
        }
    }
}

void mat_prod(double ** mat1, double ** mat2, int n, double ** prod) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double sum_loc = 0.0;
            for(int k = 0; k < n; ++k) {
                sum_loc = sum_loc + mat1[i][k] * mat2[k][j];
            }
            prod[i][j] = sum_loc;
        }
    }
}

void print_mat(double ** mat, int n) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }
}

double ** mat2D(int n) {
    double ** r = new double * [n];
    for(int i = 0; i < n; ++i) {
        r[i] = new double[n];
    }
    
    return r;
}

void compute_coeff_matrix(int flux_node,
                          c_data_t_NNNN & comp_data,
                          double * x) {
    int n = comp_data.n - 0;
    
    for(int i = 0; i < n - 1; ++i) {
        for(int j = 0; j < n - 1; ++j) {
            if(i == j) {
                double sum_loc = x[i] / (comp_data.p_params.D[i][n - 1] + 1e-13);
                double sum_frac = 0.0;
                for(int k = 0; k < n - 1; ++k) {
                    sum_frac = sum_frac + x[k];
                }
                x[n - 1] = 1.0 - sum_frac;
                for(int k = 0; k < n; ++k) {
                    if(k != i) {
                        sum_loc = sum_loc + x[k] / (comp_data.p_params.D[i][k] + 1e-13);
                    }
                }
                comp_data.coeff_mat_boundary[flux_node].A[i][j] = sum_loc;
            }
            if(i != j) {
                comp_data.coeff_mat_boundary[flux_node].A[i][j] = -x[i] * (1.0 / (comp_data.p_params.D[i][j] + 1e-13) - 1.0 / (comp_data.p_params.D[i][n - 1] + 1e-13));
            }
        }
    }
    
    compute_mat_inv(comp_data.coeff_mat_boundary[flux_node].A,
                    n - 1,
                    comp_data.coeff_mat_boundary[flux_node].A_inv);
}

void compute_coefficients(c_data_t_NNNN & comp_data) {
    int flux_node = 0;
    int ng = comp_data.ng;
    int n = comp_data.n - 1;
    
    // Bulb 1
    compute_coeff_matrix(flux_node, comp_data, comp_data.bulb_data_inter.mol_fracs_bulb1.x);
    
    // Tube
    for(flux_node = 1; flux_node < ng; ++flux_node) {
        double * x_loc = new double [n + 1];
        double sum_l_c = 0.0;
        for(int c = 0; c < n; ++c) {
            double w_term = comp_data.tube_fracs[flux_node - 1].x[c];
            double e_term = comp_data.tube_fracs[flux_node].x[c];
            x_loc[c] = 0.5 * (w_term + e_term);
            sum_l_c = sum_l_c + x_loc[c];
        }
        
        x_loc[n] = 1.0 - sum_l_c;
        
        compute_coeff_matrix(flux_node, comp_data, x_loc);
    }
    
    // Bulb 2
    compute_coeff_matrix(ng, comp_data, comp_data.bulb_data_inter.mol_fracs_bulb2.x);
    
    
//    std::cout << "ANALYZE THERE" << std::endl;
//    for(int node_l = 0; node_l < ng + 1; node_l++) {
//        std::cout << "A_inv node " << node_l << std::endl;
//        print_mat(comp_data.coeff_mat_boundary[node_l].A_inv, n);
//    }
}

void reset_coefficients(c_data_t_NNNN & comp_data) {
    int ng = comp_data.ng;
    int n = comp_data.n;
    
    for(int flux_node = 0; flux_node < ng + 1; ++flux_node) {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                comp_data.coeff_mat_boundary[flux_node].A[i][j] = 0.0;
                comp_data.coeff_mat_boundary[flux_node].A_inv[i][j] = 0.0;
            }
        }
        
    }
}

void bulb1_NNNN(c_data_t_NNNN & comp_data_N, c_data_t & comp_data_ref) {

    int flux_node = 0;
    int n = comp_data_N.n - 1;
    double V = comp_data_N.e_params.V;
    double A = comp_data_N.e_params.A;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
    double ** A_inv = comp_data_N.coeff_mat_boundary[flux_node].A_inv;
//    double ** A_ = comp_data_N.coeff_mat_boundary[flux_node].A;
//    double ** prod = mat2D(n);
    
    // Check mat
//    printf("A\n");
//    print_mat(A_, n);
//    printf("A_inv\n");
//    print_mat(A_inv, n);
//    printf("prod\n");
//    mat_prod(A_, A_inv, n, prod);
//    print_mat(prod, n);
    
    for(int i = 0; i < n; ++i) {
        double sum_loc = 0.0;
        for(int c = 0; c < n; ++c) {
            if(c != i) {
                double dxj = comp_data_N.tube_fracs[0].x[c] - comp_data_N.bulb_data.mol_fracs_bulb1.x[c];
                double dxj_dz = dxj / (0.5 * dz);
                sum_loc = sum_loc + A * A_inv[i][c] * dxj_dz;
            }
        }

        sum_loc = sum_loc + A * A_inv[i][i] * comp_data_N.tube_fracs[0].x[i] / (0.5 * dz);
        sum_loc = sum_loc + V / dt * comp_data_N.bulb_data_old.mol_fracs_bulb1.x[i];
        double ap = V / dt + A * A_inv[i][i] / (0.5 * dz);

        comp_data_N.bulb_data.mol_fracs_bulb1.x[i] = sum_loc / ap;
    }
    
    // Set last frac values
    double sum_3 = 0.0;
    for(int c = 0; c < n; ++c) {
        sum_3 = sum_3 + comp_data_N.bulb_data.mol_fracs_bulb1.x[c];
    }
    
    comp_data_N.bulb_data.mol_fracs_bulb1.x[n] = 1.0 - sum_3;
    
    // Set ref values for testing

    bulb1_c1(comp_data_ref);

    bulb1_c2(comp_data_ref);
    
//    comp_data_N.bulb_data.mol_fracs_bulb1.x[0] = comp_data_ref.bulb_data.mol_fracs_bulb1.x1;
//    comp_data_N.bulb_data.mol_fracs_bulb1.x[1] = comp_data_ref.bulb_data.mol_fracs_bulb1.x2;
    
}

void bulb2_NNNN(c_data_t_NNNN & comp_data_N, c_data_t & comp_data) {
    int ng = comp_data_N.ng;
    
    int flux_node = ng;
    int n = comp_data_N.n - 1;
    double V = comp_data_N.e_params.V;
    double A = comp_data_N.e_params.A;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
    double ** A_inv = comp_data_N.coeff_mat_boundary[flux_node].A_inv;
    
    for(int i = 0; i < n; ++i) {
        double sum_loc = 0.0;
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                double dxj_dz = (comp_data_N.bulb_data.mol_fracs_bulb2.x[j] - comp_data_N.tube_fracs[ng - 1].x[j]) / (0.5 * dz);
                sum_loc = sum_loc + A_inv[i][j] * dxj_dz;
            }
        }
        sum_loc = - A * sum_loc + A * A_inv[i][i] * comp_data_N.tube_fracs[ng - 1].x[i] / (0.5 * dz) + V / dt * comp_data_N.bulb_data_old.mol_fracs_bulb2.x[i];

        comp_data_N.bulb_data.mol_fracs_bulb2.x[i] = sum_loc / (V / dt + A * A_inv[i][i] / (0.5 * dz));
    }
    
    // Testing code below
    bulb2_c1(comp_data);
    
    bulb2_c2(comp_data);
    
//    comp_data_N.bulb_data.mol_fracs_bulb1.x[0] = comp_data.bulb_data.mol_fracs_bulb1.x1;
//    comp_data_N.bulb_data.mol_fracs_bulb1.x[1] = comp_data.bulb_data.mol_fracs_bulb1.x2;
}

void tube0_NNNN(c_data_t_NNNN & comp_data_N, c_data_t & comp_data) {
    int flux_node_w = 0;
    int flux_node_e = 1;
    int n = comp_data_N.n - 1;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
    double ** A_inv_w = comp_data_N.coeff_mat_boundary[flux_node_w].A_inv;
    double ** A_inv_e = comp_data_N.coeff_mat_boundary[flux_node_e].A_inv;
    
//    std::cout << "A_inv_w:" << std::endl;
//    print_mat(A_inv_w, n);
//
//    std::cout << "A_inv_e:" << std::endl;
//    print_mat(A_inv_e, n);
    
    for(int i = 0; i < n; ++i) {
        double sum_loc_w = 0.0;
        double sum_loc_e = 0.0;
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                double dxj_dz_w = (comp_data_N.tube_fracs[0].x[j] - comp_data_N.bulb_data_inter.mol_fracs_bulb1.x[j]) / (0.5 * dz);
                sum_loc_w = sum_loc_w + A_inv_w[i][j] * dxj_dz_w;

                double dxj_dz_e = (comp_data_N.tube_fracs[1].x[j] - comp_data_N.tube_fracs[0].x[j]) / dz;
                sum_loc_e = sum_loc_e + A_inv_e[i][j] * dxj_dz_e;
            }
        }

        double sum_loc_tot = - sum_loc_w + sum_loc_e + A_inv_w[i][i] * comp_data_N.bulb_data_inter.mol_fracs_bulb1.x[i] / (0.5 * dz) + A_inv_e[i][i] * comp_data_N.tube_fracs[1].x[i] / dz + dz / dt * comp_data_N.tube_fracs_old[0].x[i];

        double ap = dz / dt + A_inv_w[i][i] / (0.5 * dz) + A_inv_e[i][i] / dz;

        comp_data_N.tube_fracs[0].x[i] = sum_loc_tot / ap;
    }
        
    // Testing code below
    tube0_c1(comp_data);
    
    tube0_c2(comp_data);
    
//    comp_data_N.tube_fracs[0].x[0] = comp_data.tube_fracs[0].x1;
//    comp_data_N.tube_fracs[0].x[1] = comp_data.tube_fracs[0].x2;
    
//    std::cout << "A_inv_w ref:" << std::endl;
//    std::cout << beta1w << ", " << beta2w << std::endl;
//    std::cout << alpha1w << ", " << alpha2w << std::endl;
//
//    std::cout << "A_inv_e ref:" << std::endl;
//    std::cout << beta1e << ", " << beta2e << std::endl;
//    std::cout << alpha1e << ", " << alpha2e << std::endl;
}

void tube_mid_NNNN(c_data_t_NNNN & comp_data_N, c_data_t & comp_data) {
    
    int ng = comp_data_N.ng;
    
    for(int node = 1; node < ng - 1; ++node) {
        int flux_node_w = node;
        int flux_node_e = node + 1;
        int n = comp_data_N.n - 1;
        double dz = comp_data_N.e_params.dz;
        double dt = comp_data_N.t_params.dt;

        double ** A_inv_w = comp_data_N.coeff_mat_boundary[flux_node_w].A_inv;
        double ** A_inv_e = comp_data_N.coeff_mat_boundary[flux_node_e].A_inv;

        for(int i = 0; i < n; ++i) {
            double sum_loc_w = 0.0;
            double sum_loc_e = 0.0;
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    double dxj_dz_w = (comp_data_N.tube_fracs[node].x[j] - comp_data_N.tube_fracs[node - 1].x[j]) / dz;
                    sum_loc_w = sum_loc_w + A_inv_w[i][j] * dxj_dz_w;

                    double dxj_dz_e = (comp_data_N.tube_fracs[node + 1].x[j] - comp_data_N.tube_fracs[node].x[j]) / dz;
                    sum_loc_e = sum_loc_e + A_inv_e[i][j] * dxj_dz_e;
                }
            }
            double xiW = comp_data_N.tube_fracs[node - 1].x[i];
            double xiE = comp_data_N.tube_fracs[node + 1].x[i];
            double sum_loc_tot = -sum_loc_w + sum_loc_e + A_inv_w[i][i] / dz * xiW + A_inv_e[i][i] / dz * xiE + dz / dt * comp_data_N.tube_fracs_old[node].x[i];

            comp_data_N.tube_fracs[node].x[i] = sum_loc_tot / (dz / dt + A_inv_w[i][i] / dz + A_inv_e[i][i] / dz);
        }
        
//        comp_data_N.tube_fracs[node].x[0] = comp_data.tube_fracs[node].x1;
//        comp_data_N.tube_fracs[node].x[1] = comp_data.tube_fracs[node].x2;
    }
    
    // Testing code
    mid_nodes1(comp_data);
    
    mid_nodes2(comp_data);
}

void tube_n_NNNN(c_data_t_NNNN & comp_data_N, c_data_t & comp_data) {
    int ng = comp_data.ng;
    
    int flux_node_w = ng - 1;
    int flux_node_e = ng;
    int n = comp_data_N.n - 1;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
//    double * x_loc_w = new double [n + 1];
//    for(int c = 0; c < n + 1; ++c) {
//        double w_term = comp_data_N.tube_fracs_inter[ng - 2].x[c];
//        double e_term = comp_data_N.tube_fracs_inter[ng - 1].x[c];
//        x_loc_w[c] = 0.5 * (w_term + e_term);
//    }
    
    double ** A_inv_w = comp_data_N.coeff_mat_boundary[flux_node_w].A_inv;
    double ** A_inv_e = comp_data_N.coeff_mat_boundary[flux_node_e].A_inv;
    
    for(int i = 0; i < n; ++i) {
        double sum_loc_w = 0.0;
        double sum_loc_e = 0.0;
        for(int j = 0; j < n; ++j) {
            if(i != j) {
                double dxj_dz_w = (comp_data_N.tube_fracs[ng - 1].x[j] - comp_data_N.tube_fracs[ng - 2].x[j]) / dz;
                sum_loc_w = sum_loc_w + A_inv_w[i][j] * dxj_dz_w;

                double dxj_dz_e = (comp_data_N.bulb_data.mol_fracs_bulb2.x[j] - comp_data_N.tube_fracs[ng - 1].x[j]) / (0.5 * dz);
                sum_loc_e = sum_loc_e + A_inv_e[i][j] * dxj_dz_e;
            }
        }
        double xiW = comp_data_N.tube_fracs[ng - 2].x[i];
        double xiE = comp_data_N.bulb_data.mol_fracs_bulb2.x[i];
        double sum_loc_tot = -sum_loc_w + sum_loc_e + A_inv_w[i][i] / dz * xiW + A_inv_e[i][i] / (0.5 * dz) * xiE + dz / dt * comp_data_N.tube_fracs_old[ng - 1].x[i];

        comp_data_N.tube_fracs[ng - 1].x[i] = sum_loc_tot / (dz / dt + A_inv_w[i][i] / dz + A_inv_e[i][i] / (0.5 * dz));
    }
    
    
    // Testing
    tube_n_c1(comp_data);
    
    tube_n_c2(comp_data);
    
//    comp_data_N.tube_fracs[ng - 1].x[0] = comp_data.tube_fracs[ng - 1].x1;
//    comp_data_N.tube_fracs[ng - 1].x[1] = comp_data.tube_fracs[ng - 1].x2;
}

void update_composition_estimate_NNNN(c_data_t_NNNN & comp_data, c_data_t & comp_data_ref) {

    bulb1_NNNN(comp_data, comp_data_ref);
    
    tube0_NNNN(comp_data, comp_data_ref);
    
    tube_mid_NNNN(comp_data, comp_data_ref);
    
    tube_n_NNNN(comp_data, comp_data_ref);
    
    bulb2_NNNN(comp_data, comp_data_ref);
    
    set_frac_comp3_NNNN(comp_data);
}

void compute_bulb_compositions_NNNN(e_params_t e_params,
                                    p_params_t_NNNN p_params,
                                    p_params_t p_params_ref,
                                    t_params_t t_params,
                                    int ng,
                                    int n,
                                    b_data_t_NNNN & bulb_data,
                                    b_data_t & bulb_data_ref) {
    // Organize data
    c_data_t_NNNN comp_data;
    comp_data.ng = ng;
    comp_data.n = n;

    comp_data.p_params = p_params;
    comp_data.e_params = e_params;
    comp_data.t_params = t_params;

    comp_data.bulb_data = bulb_data;
    comp_data.bulb_data_old = bulb_data;
    comp_data.bulb_data_inter = bulb_data;
    
    comp_data.bulb_data.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data.mol_fracs_bulb2.x = new double[n];
    comp_data.bulb_data_old.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data_old.mol_fracs_bulb2.x = new double[n];
    comp_data.bulb_data_inter.mol_fracs_bulb1.x = new double[n];
    comp_data.bulb_data_inter.mol_fracs_bulb2.x = new double[n];
    
    for(int c = 0; c < n; ++c) {
        comp_data.bulb_data.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
        comp_data.bulb_data_old.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_old.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
        comp_data.bulb_data_inter.mol_fracs_bulb1.x[c] = bulb_data.mol_fracs_bulb1.x[c];
        comp_data.bulb_data_inter.mol_fracs_bulb2.x[c] = bulb_data.mol_fracs_bulb2.x[c];
    }
    
    // Allocate data for tube composition
    comp_data.tube_fracs = new node_t_NNNN[ng];
    comp_data.tube_fracs_inter = new node_t_NNNN[ng];
    comp_data.tube_fracs_old = new node_t_NNNN[ng];
    
    // Allocate data for flux coefficient matrices
    comp_data.coeff_mat_boundary = new c_mat_b_t[ng + 1];
    
    // Allocate space for and initialize coefficient matrix
    for(int flux_node = 0; flux_node < ng + 1; ++flux_node) {
        comp_data.coeff_mat_boundary[flux_node].A = new double * [n];
        comp_data.coeff_mat_boundary[flux_node].A_inv = new double * [n];
        for(int c = 0; c < n; ++c) {
            comp_data.coeff_mat_boundary[flux_node].A[c] = new double[n];
            comp_data.coeff_mat_boundary[flux_node].A_inv[c] = new double[n];
            for(int c_l = 0; c_l < n; ++c_l) {
                comp_data.coeff_mat_boundary[flux_node].A[c][c_l] = 0.0;
                comp_data.coeff_mat_boundary[flux_node].A_inv[c][c_l] = 0.0;
            }
        }
    }
    
    
    // Initialize tube composition
    for(int node = 0; node < ng; ++node) {
        comp_data.tube_fracs[node].x = new double[n];
        comp_data.tube_fracs_inter[node].x = new double[n];
        comp_data.tube_fracs_old[node].x = new double[n];
        for(int c = 0; c < n; ++c) {
            comp_data.tube_fracs[node].x[c] = 1.0 / n;
            comp_data.tube_fracs_inter[node].x[c] = 1.0 / n;
            comp_data.tube_fracs_old[node].x[c] = 1.0 / n;
        }
    }
    
    // Organize reference data
    c_data_t comp_data_ref;
    comp_data_ref.ng = ng;

    comp_data_ref.p_params = p_params_ref;
    comp_data_ref.e_params = e_params;
    comp_data_ref.t_params = t_params;

    comp_data_ref.bulb_data = bulb_data_ref;
    comp_data_ref.bulb_data_old = bulb_data_ref;
    comp_data_ref.bulb_data_inter = bulb_data_ref;
    
    // Allocate data for tube composition
    comp_data_ref.tube_fracs = new node_t[ng];
    comp_data_ref.tube_fracs_inter = new node_t[ng];
    comp_data_ref.tube_fracs_old = new node_t[ng];

    // Initialize tube composition
    for(int node = 0; node < ng; ++node) {
        comp_data_ref.tube_fracs[node].x1 = 1.0 / 3;
        comp_data_ref.tube_fracs[node].x2 = 1.0 / 3;
        comp_data_ref.tube_fracs[node].x3 = 1.0 / 3;

        comp_data_ref.tube_fracs_inter[node].x1 = 1.0 / 3;
        comp_data_ref.tube_fracs_inter[node].x2 = 1.0 / 3;
        comp_data_ref.tube_fracs_inter[node].x3 = 1.0 / 3;

        comp_data_ref.tube_fracs_old[node].x1 = 1.0 / 3;
        comp_data_ref.tube_fracs_old[node].x2 = 1.0 / 3;
        comp_data_ref.tube_fracs_old[node].x3 = 1.0 / 3;
    }
    
    
    // Allocate space for and initialize coefficient matrix
    comp_data_ref.coeff_mat_boundary = new c_mat_b_t[ng + 1];
    for(int flux_node = 0; flux_node < ng + 1; ++flux_node) {
        comp_data_ref.coeff_mat_boundary[flux_node].A = new double * [n];
        comp_data_ref.coeff_mat_boundary[flux_node].A_inv = new double * [n];
        for(int c = 0; c < n; ++c) {
            comp_data_ref.coeff_mat_boundary[flux_node].A[c] = new double[n];
            comp_data_ref.coeff_mat_boundary[flux_node].A_inv[c] = new double[n];
            for(int c_l = 0; c_l < n; ++c_l) {
                comp_data_ref.coeff_mat_boundary[flux_node].A[c][c_l] = 0.0;
                comp_data_ref.coeff_mat_boundary[flux_node].A_inv[c][c_l] = 0.0;
            }
        }
    }
    
    // Compute composition
    int max_out_it = MAX_OUT;
    int max_in_it = MAX_IN;
    
    double t = t_params.to;
    double dt = (t_params.tf - t_params.to) / t_params.nt;
    
    // Loop to time t = tf
    while(t < t_params.tf) {
        
        // Update bulb composition of previous time step
        for(int c = 0; c < n; ++c) {
            comp_data.bulb_data_old.mol_fracs_bulb1.x[c] = comp_data.bulb_data.mol_fracs_bulb1.x[c];
            comp_data.bulb_data_old.mol_fracs_bulb2.x[c] = comp_data.bulb_data.mol_fracs_bulb2.x[c];
        }
        
        // Update tube composition of previous time step
        for(int node = 0; node < ng; ++node) {
            for(int c = 0; c < n; ++c) {
                comp_data.tube_fracs_old[node].x[c] = comp_data.tube_fracs[node].x[c];
            }
        }
        
        // Update bulb composition of previous time step
        comp_data_ref.bulb_data_old = comp_data_ref.bulb_data;
        
        // Update tube composition of previous time step
        for(int node = 0; node < ng; ++node)
            comp_data_ref.tube_fracs_old[node] = comp_data_ref.tube_fracs[node];
        
        // Outer Gauss-Seidel iterations
        int out_it = 0;
        while(out_it < max_out_it) {
            
            // Update intermediate bulb composition
            for(int c = 0; c < n; ++c) {
                comp_data.bulb_data_inter.mol_fracs_bulb1.x[c] = comp_data.bulb_data.mol_fracs_bulb1.x[c];
                comp_data.bulb_data_inter.mol_fracs_bulb2.x[c] = comp_data.bulb_data.mol_fracs_bulb2.x[c];
            }
            
            // Update intermediate tube composition
            for(int node = 0; node < ng; ++node) {
                for(int c = 0; c < n; ++c) {
                    comp_data.tube_fracs_inter[node].x[c] = comp_data.tube_fracs[node].x[c];
                }
            }
            
            // Update intermediate bulb composition
            comp_data_ref.bulb_data_inter = comp_data_ref.bulb_data;
            
            // Update intermediate tube composition
            for(int node = 0; node < ng; ++node)
                comp_data_ref.tube_fracs_inter[node] = comp_data_ref.tube_fracs[node];
        
            reset_coefficients(comp_data);
            
            compute_coefficients(comp_data);
                 
            // Inner Gauss-Seidel iterations
            int in_it = 0;
            while(in_it < max_in_it) {
       
                // Update estimates bulb and tube composition
                update_composition_estimate_NNNN(comp_data, comp_data_ref);
                
//                // Print coefficient matrices
//                std::cout << "ANALYZE HERE" << std::endl;
//                for(int node_l = 0; node_l < ng + 1; node_l++) {
//                    std::cout << "A_inv node " << node_l << std::endl;
//                    print_mat(comp_data.coeff_mat_boundary[node_l].A_inv, n - 1);
//                    std::cout << "A_inv_ref node " << node_l << std::endl;
//                    print_mat(comp_data_ref.coeff_mat_boundary[node_l].A_inv, n - 1);
//                }
                
                in_it++;
            }
            
//            // Print current tube composition
//            std::cout << "tube fracs comps" << std::endl;
//            for(int node = 0; node < ng; ++node) {
//                for(int c = 0; c < n; ++c) {
//                    printf("%f, ", comp_data.tube_fracs[node].x[c]);
//                }
//                printf("\n");
//            }
//
//            // Print current ref tube composition
//            std::cout << "ref 3 comps" << std::endl;
//            for(int node = 0; node < ng; ++node) {
//                printf("%f, ", comp_data_ref.tube_fracs[node].x1);
//                printf("%f, ", comp_data_ref.tube_fracs[node].x2);
//                printf("%f, ", comp_data_ref.tube_fracs[node].x3);
//                printf("\n");
//            }
            
            out_it++;
        }
        
        t = t + dt;
    }
    
    // Print current tube composition
    std::cout << "tube fracs comps" << std::endl;
    for(int node = 0; node < ng; ++node) {
        for(int c = 0; c < n; ++c) {
            printf("%f, ", comp_data.tube_fracs[node].x[c]);
        }
        printf("\n");
    }

    // Print current ref tube composition
    std::cout << "ref 3 comps" << std::endl;
    for(int node = 0; node < ng; ++node) {
        printf("%f, ", comp_data_ref.tube_fracs[node].x1);
        printf("%f, ", comp_data_ref.tube_fracs[node].x2);
        printf("%f, ", comp_data_ref.tube_fracs[node].x3);
        printf("\n");
    }
    
    // Set bulb data
    bulb_data = comp_data.bulb_data;
    
    // Set bulb data
    bulb_data_ref = comp_data_ref.bulb_data;
    
    for(int node = 0; node < ng; ++node) {
        delete [] comp_data.tube_fracs[node].x;
        delete [] comp_data.tube_fracs_inter[node].x;
        delete [] comp_data.tube_fracs_old[node].x;
    }
    
    delete [] comp_data.tube_fracs;
    delete [] comp_data.tube_fracs_inter;
    delete [] comp_data.tube_fracs_old;
}

void testing() {
    
    int n = 5;
    
    double ** mat = new double * [n];
    double ** mat_inv = new double * [n];
    double ** prod = new double * [n];
    
    for(int i = 0; i < n; ++i) {
        mat[i] = new double[n];
        mat_inv[i] = new double[n];
        prod[i] = new double[n];
    }
    
    srand((unsigned) time(NULL));
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            mat[i][j] = (double) rand() / RAND_MAX * 100;
        }
    }
    
    compute_mat_inv(mat, n, mat_inv);
    
    mat_prod(mat, mat_inv, n, prod);
    
    print_mat(mat, n);
    
    print_mat(prod, n);
    
}
