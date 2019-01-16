#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

vector<double> getVector(string filename, int length) {
    // Read vector from result file:
    ifstream file;
    file.open(filename, ios::in | ios::binary);
    vector<double> v = vector<double>(0);
    int i = 0;
    double val;
    for (int i = 0; i < length; i++) {
        file.read((char*)&(val), sizeof(double));
        v.push_back(val);
    }
    file.close();
    return v;
}

