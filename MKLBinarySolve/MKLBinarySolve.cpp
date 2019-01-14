#include "MKLBinarySolve.h";

int main() {


    cout << "Done. Press any key to exit." << endl;
    cin.get();	// Press any key to continue.
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