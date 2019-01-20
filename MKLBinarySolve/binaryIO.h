#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class COOMatrix {
public:
    COOMatrix(int nnz) {
        rows = vector<int>(nnz);
        cols = vector<int>(nnz);
        vals = vector<double>(nnz);
    }

    vector<int> rows;
    vector<int> cols;
    vector<double> vals;
};

COOMatrix readCOOMatrix(string filename) {
    // Matrix must be in COO format, sorted in row-first order.
    ifstream file;
    file.open(filename, ios::in | ios::binary);
    
    int nnz;
    file.read((char*)&(nnz), sizeof(int));
    COOMatrix matrix = COOMatrix(nnz);
    
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
vector<double> readVector(string filename) {
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

void writeCOOMatrix(string filename, COOMatrix matrix) {
    ofstream file;
    file.open(filename, ios::out | ios::binary);

    int nonzeros = (int)matrix.vals.size();
    file.write((char*)&nonzeros, sizeof(int));

    for (int i = 0; i < nonzeros; i++) {
        file.write((char*)&(matrix.rows[i]), sizeof(int));
        file.write((char*)&(matrix.cols[i]), sizeof(int));
        file.write((char*)&(matrix.vals[i]), sizeof(double));
    }

    file.flush();
    file.close();
}
void writeVector(string filename, vector<double> v) {
    ofstream file;
    file.open(filename, ios::out | ios::binary);

    int n = (int)v.size();
    file.write((char*)&(n), sizeof(int));

    for (int row = 0; row < n; row++) {
        double val = v[row];
        file.write((char*)&(val), sizeof(double));
    }

    file.flush();
    file.close();
}