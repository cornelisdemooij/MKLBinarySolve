#include "stdafx.h"
#include "binaryIO.h"
#include "PardisoSolver.h"

using namespace std;

int main(int argc, char **argv) {
    string AFilename = argv[1];
    string bFilename = argv[2];
    string xFilename = argv[3];

    COOMatrix AMatrix = readCOOMatrix(AFilename);
    vector<double> bVector = readVector(bFilename);

    vector<int> ARows = AMatrix.rows;
    vector<int> ACols = AMatrix.cols;
    vector<double> AVals = AMatrix.vals;
    vector<double> bVals = bVector;

    //vector<double> xVector = readVector(xFilename);

    /*vector<int> ARows = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 };
    vector<int> ACols = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
    vector<double> AVals = { 0.3, 4.5, -7.3, 8.1, 0.2, -1.8, 9.7, -2.7, 6.6, -0.3, 10.9, 8.7, -1.1, 6.5, -4.1, 8.9 };
    vector<double> bVals = { 1, 0.1, 0.01, 0.001 };*/
    // Solution should be: x = { -3.937, -2.975, 0.746, 1.952 };

    vector<double> x = solveWithPardiso(ARows, ACols, AVals, bVals);

    writeVector(xFilename, x);

    //cout << "Done. Press any key to exit." << endl;
    //cin.get();	// Press any key to continue.
    exit(0);
}

/*
function [] = binaryReadSolveWrite(KFilename, fFilename, outFilename)
    KFile = fopen(KFilename);
    h = fread(KFile,1,'int');
    w = fread(KFile,1,'int');
    rows = [];
    cols = [];
    vals = [];
    count = 0;
    while(~feof(KFile))
        rowMaybe = fread(KFile,1,'int');
        if (isempty(rowMaybe))
            break;
        end
        count = count + 1;
        rows(count) = rowMaybe+1;
        cols(count) = fread(KFile,1,'int')+1;
        vals(count) = fread(KFile,1,'double');
        %disp([num2str(rows(count)) ' ' num2str(cols(count)) ' ' num2str(vals(count))]);
    end
    Ks = sparse(rows,cols,vals,h,w);
    %size(Ks)
    fclose(KFile);
    
    fFile = fopen(fFilename);
    vals = [];
    count = 0;
    while(~feof(fFile))
        val = fread(fFile,1,'double');
        if (isempty(val))
            break;
        end
        count = count + 1;
        vals(count) = val;
    end
    f = vals';
    %size(f)
    fclose(fFile);
    
    x = Ks\f;
    %disp(x(1:24))
    
    xFile = fopen(outFilename,'w');
    fwrite(xFile,x,'double');
    fclose(xFile);
end
*/