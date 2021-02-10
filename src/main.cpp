#include <iostream>

#include "fp.h"

using namespace std;

int main(int argc, char** argv)
{
    if( argc != 3 ) {
        cout << "Usage:\n"
             << "\t./fp InputFile OutputFile\n";
        return 0;
    }

    FP fp;
    fp.readHeatMap( string(argv[1]) );

    return 1;
}