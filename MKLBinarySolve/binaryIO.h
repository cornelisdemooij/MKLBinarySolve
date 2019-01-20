#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class COOMatrix {
public:
    vector<int> rows;
    vector<int> cols;
    vector<double> vals;
};

COOMatrix getCOOMatrix(string filename) {
    // Read matrix from file. Matrix must be in COO format, sorted in row-first order.
    ifstream file;
    file.open(filename, ios::in | ios::binary);
    
    int nnz;
    file.read((char*)&(nnz), sizeof(int));

    COOMatrix matrix = COOMatrix();
    matrix.rows = vector<int>(nnz);
    matrix.cols = vector<int>(nnz);
    matrix.vals = vector<double>(nnz);
    
    int row;
    int col;
    double val;
    for (int i = 0; i < nnz; i++) {
        file.read((char*)&(row), sizeof(int));
        file.read((char*)&(col), sizeof(int));
        file.read((char*)&(val), sizeof(double));

        matrix.rows[i] = row;
        matrix.cols[i] = col;
        matrix.vals[i] = val;
    }

    file.close();
    return matrix;
}
vector<double> getVector(string filename) {
    // Read vector from file.
    ifstream file;
    file.open(filename, ios::in | ios::binary);

    int n;
    file.read((char*)&(n), sizeof(int));

    vector<double> v = vector<double>(n);

    double val;
    for (int i = 0; i < n; i++) {
        file.read((char*)&(val), sizeof(double));
        v[i] = val;
    }

    file.close();
    return v;
}
