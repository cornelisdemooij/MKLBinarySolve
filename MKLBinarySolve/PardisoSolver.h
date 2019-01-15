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

    // PARDISO input parameters:
    //int nnz = ia[n];      // Number of nonzeros.
    const int mtype = 11;   // Real unsymmetric matrix.
    double *b = new double[n];  // RHS vector.
    double *x = new double[n];  // Solution vector.
    int nrhs = 1;           // Number of right hand sides.
    void *pt[64];           // Internal solver memory pointer pt. 32-bit: int pt[64]; 64-bit: long int pt[64]. Both: void *pt[64].
    for (int i = 0; i < 64; i++) {
        pt[i] = 0;
    }
    // PARDISO control parameters:
    MKL_INT iparm[64];
    //double dparm[64];
    int solver;
    MKL_INT maxfct, mnum, phase, error, msglvl;
    //int num_procs;          // Number of processors.
    //char *var;              // Auxiliary variable.
    //size_t len = 0;
    double ddum;            // Double dummy variable.
    int idum;               // Integer dummy variable.

    // Set up PARDISO control parameters:
    error = 0;
    solver = 0;     // Use sparse direct solver.
    pardisoinit(pt, &mtype, &solver);// , iparm, dparm, &error);

    //_dupenv_s(&var, &len, "OMP_NUM_THREADS");   // var = getenv("OMP_NUM_THREADS");
    //if (var != NULL) {
    //    sscanf_s(var, "%d", &num_procs);        // sscanf(var, "%d", &num_procs);
    //}
    //else {
    //    printf("Set environment OMP_NUM_THREADS to 1");
    //    exit(1);
    //}

    for (int i = 0; i < 64; i++) {
    iparm[i] = 0;
    }
    //iparm[2] = 1;   // num_procs;
    iparm[0] = 1;           // No solver default
    iparm[1] = 2;           // Fill-in reordering from METIS
    // Numbers of processors, value of OMP_NUM_THREADS
    iparm[2] = 4;
    iparm[3] = 0;           // No iterative-direct algorithm
    iparm[4] = 0;           // No user fill-in reducing permutation
    iparm[5] = 0;           // Write solution into x
    iparm[6] = 0;           // Not in use
    iparm[7] = 2;           // Max numbers of iterative refinement steps
    iparm[8] = 0;           // Not in use
    iparm[9] = 13;          // Perturb the pivot elements with 1E-13
    iparm[10] = 1;          // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 0;          // Conjugate transposed/transpose solve
    iparm[12] = 1;          // Maximum weighted matching algorithm is switched-on (default for non-symmetric)
    iparm[13] = 0;          // Output: Number of perturbed pivots
    iparm[14] = 0;          // Not in use
    iparm[15] = 0;          // Not in use
    iparm[16] = 0;          // Not in use
    iparm[17] = -1;         // Output: Number of nonzeros in the factor LU
    iparm[18] = -1;         // Output: Mflops for LU factorization
    iparm[19] = 0;          // Output: Numbers of CG Iterations
    //iparm[60 - 1] = 1;      // out of core version
    //iparm[26] = 0;
    //iparm[28] = 0;
    iparm[34] = 1;
    iparm[59] = 2;  // Use out of core (OOC) version.

    maxfct = 1;     // Maximum number of numerical factorizations.
    mnum = 1;       // Which factorization to use.

    msglvl = 1;     // Print statistical information.
    error = 0;      // Initialize error flag.

    // Convert matrix from 0-based C-notation to Fortran 1-based notation:
    cout << "a: " << endl;
    for (int i = 0; i < nnz; i++) {
        cout << to_string(a[i]) << " ";
    }
    cout << endl << "ia: " << endl;
    for (int i = 0; i < n + 1; i++) {
        //ia[i] += 1;
        cout << to_string(ia[i]) << " ";
    }
    cout << endl << "ja: " << endl;
    for (int i = 0; i < nnz; i++) {
        //ja[i] += 1;
        cout << to_string(ja[i]) << " ";
    }
    cout << endl;

    // Copy contents of b vector to RHS b array:
    for (int i = 0; i < n; i++) {
        b[i] = bVals[i];
    }

    // Reordering and Symbolic Factorization. This step also allocates all memory for the factorization:
    phase = 11;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
    &n, a, ia, ja, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
    printf("\nERROR during symbolic factorization: %d", error);
    exit(1);
    }
    printf("\nReordering completed ...\n");
    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

    // Numerical factorization:
    phase = 22;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
    &n, a, ia, ja, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
    exit(2);
    }
    printf("\nFactorization completed ...\n");

    // Back substitution and iterative refinement:
    phase = 33;
    iparm[7] = 1;       // Max numbers of iterative refinement steps.
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
    &n, a, ia, ja, &idum, &nrhs,
    iparm, &msglvl, b, x, &error);
    if (error != 0) {
    printf("\nERROR during solution: %d", error);
    exit(3);
    }
    printf("\nSolve completed ...\n");
    printf("\nThe solution of the system is: ");
    for (int i = 0; i < n; i++) {
    printf("\n x [%d] = %f", i, x[i]);
    }
    printf("\n");

    // Back substitution with transposed matrix A^t x = b
    phase = 33;
    iparm[7] = 1;       // Max numbers of iterative refinement steps.
    iparm[11] = 1;      // Solving with transpose matrix.
    // Copy contents of b vector to RHS b array:
    for (int i = 0; i < n; i++) {
    b[i] = bVals[i];
    }
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
    &n, a, ia, ja, &idum, &nrhs,
    iparm, &msglvl, b, x, &error);
    if (error != 0) {
    printf("\nERROR during solution: %d", error);
    exit(3);
    }
    printf("\nSolve completed ...\n");
    printf("\nThe solution of the system is: ");
    for (int i = 0; i < n; i++) {
    printf("\n x [%d] = %f", i, x[i]);
    }
    printf("\n");

    // Convert matrix back to 0-based C-notation:
    for (int i = 0; i < n + 1; i++) {
    //ia[i] -= 1;
    }
    for (int i = 0; i < nnz; i++) {
    //ja[i] -= 1;
    }

    // Termination and release of memory:
    phase = -1;     // Release internal memory.
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
    &n, &ddum, ia, ja, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);

    // Convert x to vector and return it:
    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = x[i];
    }
    return result;
}