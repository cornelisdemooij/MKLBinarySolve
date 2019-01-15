#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mkl_pardiso.h"
#include "mkl_types.h"

using namespace std;

vector<double> solveWithPardiso(vector<int> ARows, vector<int> ACols, vector<double> AVals, vector<double> bVals) {
    // Convert the input to Pardiso-compatible CSC format:
    int n0 = (int)bVals.size();       // Height & width of the matrix.
    int nnz = (int)AVals.size();     // Number of nonzero values.

    int *ia0 = new int[n0 + 1];
    int *ja0 = new int[nnz];
    double *a0 = new double[nnz];

    // Construct ia:
    ia0[0] = 0; // 1; // Uncomment for 1-based indexing.
    ia0[1] = 0; // 1; // Uncomment for 1-based indexing.
    int iia = 1;
    for (int row : ARows) {
        while (row > iia - 1) {
            iia++;
            ia0[iia] = ia0[iia - 1];
        }
        ia0[iia]++;
    }

    // Copy ACols to ja:
    for (int i = 0; i < nnz; i++) {
        ja0[i] = ACols[i];  // +1; // Uncomment for 1-based indexing.
        a0[i] = AVals[i];
    }

    


    /* Matrix data in CSC format. */
    MKL_INT n = 4;
    const MKL_INT * ia = ia0;// { 1, 4, 7, 9, 12, 14};
    for (int i = 0; i < n0 + 1; i++) {
        cout << to_string(ia[i]) << endl;
    }

    const MKL_INT * ja = ja0;
    for (int i = 0; i < nnz; i++) {
        cout << to_string(ja[i]) << endl;
    }

    const double * a = a0;
    for (int i = 0; i < nnz; i++) {
        cout << to_string(a[i]) << endl;
    }

    MKL_INT mtype = 11;       /* Real unsymmetric matrix */
    /* RHS and solution vectors. */
    double b[4], x[4], bs[4], res, res0;
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i, j;
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    char *uplo;
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ ) {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    iparm[26] = 1;          // 1 = matrix checker on.
    iparm[34] = 1;          // Use 0-based indexing.
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ ) {
        pt[i] = 0;
    }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ) {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ) {
        printf ("\nERROR during numerical factorization: %d", error);
        exit (2);
    }
    printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Solution phase. */
/* -------------------------------------------------------------------- */
    phase = 33;
    /* Set right hand side to one. */
    for ( i = 0; i < n; i++ ) {
        b[i] = 1;
    }

// Transpose solve is used for systems in CSC format
    iparm[11] = 2;

    printf ("\n\nSolving the system in CSC format...\n");
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 ) {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }

    printf ("\nThe solution of the system is: ");
    for ( j = 0; j < n; j++ ) {
        printf ("\n x [%d] = % f", j, x[j]);
    }
    printf ("\n");

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);




    // Convert x to vector and return it:
    vector<double> result(4);// (n);
    for (int i = 0; i < 4; i++) {
        result[i] = 0.0;// x[i];
    }
    return result;
}