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

void get_order(double ** mat, int n, double * order_arr) {
    for(int row = 0; row < n; ++row) {
        int order = 0;
        for(int c = 0; c < n; ++c) {
            if(mat[row][c] == 0) {
                order++;
            }
        }
        order_arr[row] = order;
    }
}

void merge(oa_elem_t A[], int p, int q, int r) {
    int size_r, size_l;
    int i, j;
    size_l = q - p + 1;
    size_r = r - q;
    oa_elem_t L[size_l+1];
    oa_elem_t R[size_r+1];
    i = 0;
    j = 0;
    for(int n = p; n < q + 1; ++n) {
        L[i] = A[n];
        ++i;
    }
    L[size_l].val = MAX_INT;
    for(int n = q + 1; n < r + 1; ++n) {
        R[j] = A[n];
        ++j;
    }
    R[size_r].val = MAX_INT;
    i = 0;
    j = 0;
    for(int n = p; n < r + 1; ++n) {
        if(L[i].val < R[j].val) {
            A[n] = L[i];
            ++i;
        }
        else {
            A[n] = R[j];
            ++j;
        }
    }
}

void merge_sort(oa_elem_t A[], int p, int r) {
    int q;
    if(p < r) {
        q = (p + r)/2;
        merge_sort(A, p, q);
        merge_sort(A, q + 1, r);
        merge(A, p, q, r);
    }
}

void mergesort(oa_elem_t A[], int size) {
    merge_sort(A, 0, size - 1);
}



void make_ordered_mat(double ** mat, int n, double * order_arr, double ** ordered_mat) {

    oa_elem_t * order_array = new oa_elem_t[n];

    for(int row = 0; row < n; ++row) {
        order_array[row].old_row = row;
        order_array[row].val = order_arr[row];
    }

    mergesort(order_array, n);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            int old_row = order_array[row].old_row;
            ordered_mat[row][c] = mat[old_row][c];
        }
    }
}

int count_leading_zeros(double ** mat, int n, int row) {

    int count_lz = 0;

    for(int c = 0; c < n; ++c) {
        if(mat[row][c] == 0) {
            count_lz++;
        }
        else {
            break;
        }
    }

    return count_lz;
}

void mat_mult_sq(double ** A, double ** A_inv, int n, double ** mat_res) {

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            double sum_loc = 0;

            for(int k = 0; k < n; ++k) {
                sum_loc = sum_loc + A[i][k] * A_inv[k][j];
            }

            mat_res[i][j] = sum_loc;
        }
    }
}


void gauss_jordan(double ** mat, int n, double ** mat_inv) {

    double ** mat_check = mat2D(n);
    double ** mat_check_ordered = mat2D(n);
    double ** mat_ref = mat2D(n);
    double ** mat_ordered = mat2D(n);
    double ** mat_inv_ordered = mat2D(n);
    double * order_arr = new double[n];


    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            if(c == row) {
                mat_inv[row][c] = 1.0;
            }
            else {
                mat_inv[row][c] = 0.0;
            }
        }
    }

    get_order(mat, n, order_arr);

    make_ordered_mat(mat, n, order_arr, mat_ordered);

    make_ordered_mat(mat, n, order_arr, mat_check_ordered);

    make_ordered_mat(mat_inv, n, order_arr, mat_inv_ordered);

    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            mat_check[row][c] = mat_check_ordered[row][c];
            mat_ref[row][c] = mat_ordered[row][c];
            mat_inv[row][c] = mat_inv_ordered[row][c];
        }
    }

    // Convert to row echelon form
    for(int c = 0; c < n; ++c) {

        // Sort if under threshold
        if(fabs(mat_ref[c][c]) < 1e-8) {
            get_order(mat_ref, n, order_arr);

            make_ordered_mat(mat_ref, n, order_arr, mat_ordered);

            make_ordered_mat(mat_inv, n, order_arr, mat_inv_ordered);

            make_ordered_mat(mat_check, n, order_arr, mat_check_ordered);

            for(int row = 0; row < n; ++row) {
                for(int c = 0; c < n; ++c) {
                    mat_ref[row][c] = mat_ordered[row][c];
                    mat_inv[row][c] = mat_inv_ordered[row][c];
                    mat_check[row][c] = mat_check_ordered[row][c];
                }
            }
        }

        // Normalize matrix row
        for(int col = c + 1; col < n; ++col) {
            mat_ref[c][col] = mat_ref[c][col] / (mat_ref[c][c] + 1e-13);
        }

        // Update row matrix inverse
        for(int col = 0; col < n; ++col) {
            mat_inv[c][col] = mat_inv[c][col] / (mat_ref[c][c] + 1e-13);
        }

        mat_ref[c][c] = 1.0;

        // Delete elements in rows below
        for(int row = c + 1; row < n; ++row) {
            if(mat_ref[row][c] != 0) {
                for(int col = c + 1; col < n; ++col) {
                    mat_ref[row][col] = -1.0 * mat_ref[row][c] * mat_ref[c][col] + mat_ref[row][col];
                }
                for(int col = 0; col < n; ++col) {
                    mat_inv[row][col] = -1.0 * mat_ref[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat_ref[row][c] = 0;
            }
        }

    }

    // Backtracing
    for(int c = n - 1; c > 0; --c) {
        for(int row = c - 1; row > -1; row--) {
            if(mat_ref[row][c] != 0) {
                for(int col = 0; col < n; col++) {
                    mat_inv[row][col] = -1.0 * mat_ref[row][c] * mat_inv[c][col] + mat_inv[row][col];
                }
                mat_ref[row][c] = 0;
            }
        }
    }

}

void compute_mat_inv2(double ** mat, int n, double ** mat_inv) {

    gauss_jordan(mat, n, mat_inv);
}

void init_diffusivities(p_params_t & p_params, int n) {
    p_params.D = mat2D(n);
    
    for(int j = 0; j < n; ++j) {
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

double determinant(double ** a, int n) {
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

           det += pow(-1.0, j1 + 2.0) * a[0][j1] * determinant(m, n - 1);

           for (i = 0; i < n - 1; i++)
               delete [] m[i];

           delete [] m;
      }
   }

   return det;

}

void co_factor(double ** a, int n, double ** b)
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
            det = determinant(c, n - 1);

            // Fill in the elements of the cofactor
            b[i][j] = pow(-1.0, i + j + 2.0) * det;
        }
    }

    for (i = 0; i < n - 1; i++)
        delete [] c[i];

    delete [] c;
}


void transpose(double ** a, int n) {
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
    double det = determinant(mat, n);


    // If matrix is singular, exit
    if(!det) {
        printf("\nTransformation matrix is singular\n");
        exit(4);
    }


    // Calculating co-factor M_adj of M
    co_factor(mat, n, mat_adj);


    // Transposing M_adj
    transpose(mat_adj, n);


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
    
    compute_mat_inv2(comp_data.coeff_mat_boundary[flux_node].A,
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
        comp_data.coeff_mat_boundary[flux_node].A = mat2D(n);
        comp_data.coeff_mat_boundary[flux_node].A_inv = mat2D(n);
        for(int c = 0; c < n; ++c) {
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
    for(int c = 0; c < n; ++c) {
        bulb_data.mol_fracs_bulb1.x[c] = comp_data.bulb_data.mol_fracs_bulb1.x[c];
        bulb_data.mol_fracs_bulb2.x[c] = comp_data.bulb_data.mol_fracs_bulb2.x[c];
    }
    
    // Free allocated space
    delete [] comp_data.bulb_data.mol_fracs_bulb1.x;
    delete [] comp_data.bulb_data.mol_fracs_bulb2.x;
    delete [] comp_data.bulb_data_old.mol_fracs_bulb1.x;
    delete [] comp_data.bulb_data_old.mol_fracs_bulb2.x;
    delete [] comp_data.bulb_data_inter.mol_fracs_bulb1.x;
    delete [] comp_data.bulb_data_inter.mol_fracs_bulb2.x;
    
    for(int flux_node = 0; flux_node < ng + 1; ++flux_node) {
        free_mat2D(comp_data.coeff_mat_boundary[flux_node].A, n);
        free_mat2D(comp_data.coeff_mat_boundary[flux_node].A_inv, n);
    }
    
    delete [] comp_data.coeff_mat_boundary;
    
    for(int node = 0; node < ng; ++node) {
        delete [] comp_data.tube_fracs[node].x;
        delete [] comp_data.tube_fracs_inter[node].x;
        delete [] comp_data.tube_fracs_old[node].x;
    }
    
    delete [] comp_data.tube_fracs;
    delete [] comp_data.tube_fracs_inter;
    delete [] comp_data.tube_fracs_old;
}

