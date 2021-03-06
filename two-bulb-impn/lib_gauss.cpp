//
//  lib_gauss.cpp
//  two-bulb-impn
//
//  Created by Derek Harrison on 29/04/2022.
//

#include <math.h>

#include "lib.hpp"
#include "lib_gauss.hpp"
#include "user_types.h"

void get_order(double ** mat, int n, double * order_arr) {
    for(int row = 0; row < n; ++row) {
        int order = 0;
        while(fabs(mat[row][order]) <= SMALL_NUM && order < n) {
            mat[row][order] = 0.0;
            order++;
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
            if(fabs(ordered_mat[row][c]) <= SMALL_NUM) {
                ordered_mat[row][c] = 0.0;
            }
        }
    }

    delete [] order_array;
}

int count_leading_zeros(double ** mat, int n, int row) {

    int count = 0;

    while(fabs(mat[row][count]) <= SMALL_NUM && count < n) {
        count++;
    }

    return count;
}

void singularity_check(double ** mat_ref, int n, bool & mat_is_singular) {
    // Check for singularity
    for(int row = 0; row < n; ++row) {
        bool all_zeros_r = true;
        for(int col = 0; col < n; ++col) {
            if(fabs(mat_ref[row][col]) > SMALL_NUM)
                all_zeros_r = false;
        }
        if(all_zeros_r && !mat_is_singular) {
            printf("Matrix is singular\n");
            mat_is_singular = true;
        }
    }

    for(int col = 0; col < n; ++col) {
        bool all_zeros_c = true;
        for(int row = 0; row < n; ++row) {
            if(fabs(mat_ref[row][col]) > SMALL_NUM)
                all_zeros_c = false;
        }
        if(all_zeros_c && !mat_is_singular) {
            printf("Matrix is singular\n");
            mat_is_singular = true;
        }
    }
}

void gauss_jordan(double ** mat, int n, double ** mat_inv) {

    double ** mat_ref = mat2D(n);
    double ** mat_ordered = mat2D(n);
    double ** mat_inv_ordered = mat2D(n);
    double * order_arr = new double[n];

    // Initialize matrix inverse
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

    // Get row orders of input matrix
    get_order(mat, n, order_arr);

    // Get sorted matrices
    make_ordered_mat(mat, n, order_arr, mat_ordered);

    make_ordered_mat(mat_inv, n, order_arr, mat_inv_ordered);

    // Set matrices to ordered matrices
    for(int row = 0; row < n; ++row) {
        for(int c = 0; c < n; ++c) {
            mat_ref[row][c] = mat_ordered[row][c];
            mat_inv[row][c] = mat_inv_ordered[row][c];
        }
    }

    // Initialize singularity flag
    bool mat_is_singular = false;

    // Check if input matrix is singular
    singularity_check(mat_ref, n, mat_is_singular);

    // Convert to row echelon form
    for(int c = 0; c < n; ++c) {

        // Sort if under threshold
        if(fabs(mat_ref[c][c]) <= SMALL_NUM) {
            // Get row orders of coefficient matrix
            get_order(mat_ref, n, order_arr);

            // Get sorted matrices
            make_ordered_mat(mat_ref, n, order_arr, mat_ordered);

            make_ordered_mat(mat_inv, n, order_arr, mat_inv_ordered);

            // Set matrices to ordered matrices
            for(int row = 0; row < n; ++row) {
                for(int c = 0; c < n; ++c) {
                    mat_ref[row][c] = mat_ordered[row][c];
                    mat_inv[row][c] = mat_inv_ordered[row][c];
                }
            }
        }

        // Normalize matrix row
        for(int col = c + 1; col < n; ++col) {
            mat_ref[c][col] = mat_ref[c][col] / (mat_ref[c][c] + SMALL_NUM);
        }

        // Update row matrix inverse
        for(int col = 0; col < n; ++col) {
            mat_inv[c][col] = mat_inv[c][col] / (mat_ref[c][c] + SMALL_NUM);
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

            int num_lead_zeros = count_leading_zeros(mat_ref, n, row);

            if(num_lead_zeros == n && !mat_is_singular) {
                printf("Matrix is singular\n");
                mat_is_singular = true;
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

    // Check if matrix is singular
    singularity_check(mat_ref, n, mat_is_singular);

    // Free allocated space
    free_mat2D(mat_ref, n);
    free_mat2D(mat_ordered, n);
    free_mat2D(mat_inv_ordered, n);
    delete [] order_arr;
}
