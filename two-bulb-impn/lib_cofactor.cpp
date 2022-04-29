//
//  lib_cofactor.cpp
//  twin-bulb-impn
//
//  Created by Derek Harrison on 29/04/2022.
//

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "lib_cofactor.hpp"

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
