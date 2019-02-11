#pragma once
#include <vector>
#include <iostream>
#include <string>

#include "mkl_pardiso.h"

using namespace std;

vector<double> solveWithPardiso(vector<int> ARows, vector<int> ACols, vector<double> AVals, vector<double> bVals) {
    // Convert the input to Pardiso-compatible CSC format:
    int n = (int)bVals.size();       // Height & width of the matrix.
    int nnz = (int)AVals.size();     // Number of nonzero values.

    int *ia = new int[n + 1];
    int *ja = new int[nnz];
    double *a = new double[nnz];

    // Construct ia:
    ia[0] = 0;
    ia[1] = 0;
    int iia = 1;
    for (int row : ARows) {
        while (row > iia - 1) {
            iia++;
            ia[iia] = ia[iia - 1];
        }
        ia[iia]++;
    }

    // Copy ACols to ja:
    for (int i = 0; i < nnz; i++) {
        ja[i] = ACols[i];
        a[i] = AVals[i];
    }

    double *b = new double[n];  // RHS vector.
    double *x = new double[n];  // Solution vector.

    // PARDISO input parameters:
    MKL_INT mtype = 11;         // Real unsymmetric matrix.
    double bs[4], res, res0;
    MKL_INT nrhs = 1;           // Number of right hand sides.
    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures:
    void *pt[64];
    
    // Set up Pardiso control parameters:
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT i, j;
    double ddum;        // Double dummy.
    MKL_INT idum;       // Integer dummy.
    char *uplo;
    for ( i = 0; i < 64; i++ ) {
        iparm[i] = 0;
    }
    iparm[0] = 1;       // No solver default.
    iparm[1] = 2;       // Fill-in reordering from METIS.
    iparm[2] = 4;       // Numbers of processors, value of OMP_NUM_THREADS.
    iparm[3] = 0;       // No iterative-direct algorithm.
    iparm[4] = 0;       // No user fill-in reducing permutation.
    iparm[5] = 0;       // Write solution into x.
    iparm[6] = 0;       // Not in use.
    iparm[7] = 4;// 2;       // Max numbers of iterative refinement steps.
    iparm[8] = 0;       // Not in use.
    iparm[9] = 13;      // Perturb the pivot elements with 1E-13.
    iparm[10] = 1;      // Use nonsymmetric permutation and scaling MPS.
    iparm[11] = 0;      // Conjugate/transpose solve.
    iparm[12] = 2;// 1;      // Maximum weighted matching algorithm is switched-on (default for non-symmetric).
    iparm[13] = 0;      // Output: Number of perturbed pivots.
    iparm[14] = 0;      // Not in use.
    iparm[15] = 0;      // Not in use.
    iparm[16] = 0;      // Not in use.
    iparm[17] = -1;     // Output: Number of nonzeros in the factor LU.
    iparm[18] = -1;     // Output: Mflops for LU factorization.
    iparm[19] = 0;      // Output: Numbers of CG Iterations.
    iparm[26] = 1;      // 1 = matrix checker on.
    iparm[34] = 1;      // Use 0-based indexing.
    iparm[59] = 2;      // Use out of core (OOC) version.

    maxfct = 1;         // Maximum number of numerical factorizations.
    mnum = 1;           // Which factorization to use.
    msglvl = 0;// 1;         // Print statistical information in file.
    error = 0;          // Initialize error flag.
    
    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver:
    for ( i = 0; i < 64; i++ ) {
        pt[i] = 0;
    }
    
    // Reordering and Symbolic Factorization. This step also allocates
    // all memory that is necessary for the factorization:
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, 
             iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ) {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
    //printf ("\nReordering completed ... ");
    //printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    //printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
    
    // Numerical factorization:
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, 
             iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ) {
        printf ("\nERROR during numerical factorization: %d", error);
        exit (2);
    }
    //printf ("\nFactorization completed ... ");
    
    // Solution phase:
    phase = 33;
    for ( i = 0; i < n; i++ ) {
        b[i] = bVals[i];
    }

    // Transpose solve is used for systems in CSC format
    iparm[11] = 2;      // Transpose solve.
    //printf ("\n\nSolving the system in CSC format...\n");
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, 
             iparm, &msglvl, b, x, &error);
    if (error != 0) {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }

    //printf ("\nThe solution of the system is: ");
    for (j = 0; j < n; j++) {
        //printf ("\n x [%d] = % f", j, x[j]);
    }
    //printf ("\n");

    // Termination and release of memory:
    phase = -1;         // Release internal memory.
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);

    // Convert x to vector and return it:
    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = x[i];
    }
    return result;
}