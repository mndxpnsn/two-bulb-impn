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


void init_diffusivities(p_params_t & p_params, int n) {
    p_params.D = new double * [n];
    for(int j = 0; j < n; ++j) {
        p_params.D[j] = new double[n];
        for(int c = 0; c < n; ++c) {
            p_params.D[j][c] = 0.0;
        }
    }
}

void set_frac_comp3(c_data_t & comp_data) {
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

void free_mat2D(double ** mat, int n) {

    for(int i = 0; i < n; ++i) {
        delete [] mat[i];
    }
    
    delete [] mat;
}

void compute_coeff_matrix(int flux_node,
                          c_data_t & comp_data,
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

void compute_coefficients(c_data_t & comp_data) {
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
}

void reset_coefficients(c_data_t & comp_data) {
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

void bulb1(c_data_t & comp_data_N) {

    int flux_node = 0;
    int n = comp_data_N.n - 1;
    double V = comp_data_N.e_params.V;
    double A = comp_data_N.e_params.A;
    double dz = comp_data_N.e_params.dz;
    double dt = comp_data_N.t_params.dt;
    
    double ** A_inv = comp_data_N.coeff_mat_boundary[flux_node].A_inv;
    
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
}

void bulb2(c_data_t & comp_data_N) {
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
}

void tube0(c_data_t & comp_data_N) {
    int flux_node_w = 0;
    int flux_node_e = 1;
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
}

void tube_mid(c_data_t & comp_data_N) {
    
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
    }
}

void tube_n(c_data_t & comp_data_N) {
    int ng = comp_data_N.ng;
    
    int flux_node_w = ng - 1;
    int flux_node_e = ng;
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
}

void comp_data_ref(c_data_t & comp_data) {

    bulb1(comp_data);
    
    tube0(comp_data);
    
    tube_mid(comp_data);
    
    tube_n(comp_data);
    
    bulb2(comp_data);
    
    set_frac_comp3(comp_data);
}

void compute_bulb_compositions(e_params_t e_params,
                               p_params_t p_params,
                               t_params_t t_params,
                               int ng,
                               int n,
                               b_data_t & bulb_data) {
    // Organize data
    c_data_t comp_data;
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
    comp_data.tube_fracs = new node_t[ng];
    comp_data.tube_fracs_inter = new node_t[ng];
    comp_data.tube_fracs_old = new node_t[ng];
    
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
            
            reset_coefficients(comp_data);
            
            compute_coefficients(comp_data);
                 
            // Inner Gauss-Seidel iterations
            int in_it = 0;
            while(in_it < max_in_it) {
       
                // Update estimates bulb and tube composition
                comp_data_ref(comp_data);
                
                in_it++;
            }
            
            out_it++;
        }
        
        t = t + dt;
    }
    
    // Set bulb data
    bulb_data = comp_data.bulb_data;
    
    for(int node = 0; node < ng; ++node) {
        delete [] comp_data.tube_fracs[node].x;
        delete [] comp_data.tube_fracs_inter[node].x;
        delete [] comp_data.tube_fracs_old[node].x;
    }
    
    delete [] comp_data.tube_fracs;
    delete [] comp_data.tube_fracs_inter;
    delete [] comp_data.tube_fracs_old;
}

