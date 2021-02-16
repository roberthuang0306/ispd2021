#include "fp.h"

FP::FP()
{
    _Dimension = 0;
    _HeatMap = 0;
    for(int i = 0; i < D; ++i) _HeatMapSize[i] = 0;
    for(int i = 0; i < 2; ++i) _PESize[i] = 0;
    for(int i = 0; i < 2; ++i) _CostCoefficients[i] = 0.0;
    _S = 0;
    _Tiles = 0;
    _n_Tiles = 0;
}

FP::~FP()
{
    if( _HeatMap ) {
        for(int i = 0; i < _HeatMapSize[0] + 1; ++i) {
            for(int j = 0; j < _HeatMapSize[1] + 1; ++j) {
                delete [] _HeatMap[i][j];
            }
            delete [] _HeatMap[i];
        }
        delete [] _HeatMap;
    }
    if( _Tiles ) {
        delete [] _Tiles;
    }
}

int FP::readHeatMap(string filename)
{
    ifstream file(filename);
    if(!file.is_open()){
        cout << "Unable to open file \"" << filename << "\"\n";
        return 0;
    }
    string line;
    int n_line = 0;
    size_t start, end;
    int x = 0, y = 0, z = 0;
    while(getline(file, line)) {
        if(line[0] == '#') continue;

        if(n_line == 0) {
            // Dimension
            _Dimension = stoi(line);
        }
        else if(n_line == 1) {
            // HeatMap size
            start = 0;
            for(int i = 0; i < _Dimension; ++i) {
                end = line.find(' ', start);
                _HeatMapSize[i] = stoi( line.substr(start, end - start) );
                start = end + 1;
            }
            _HeatMap = new double**[_HeatMapSize[0] + 1];
            for(int i = 0; i < _HeatMapSize[0] + 1; ++i) {
                _HeatMap[i] = new double*[_HeatMapSize[1] + 1];
                for(int j = 0; j < _HeatMapSize[1] + 1; ++j) {
                    _HeatMap[i][j] = new double[_HeatMapSize[2] + 1];
                }
            }
        }
        else if(n_line == 2) {
            // PE size
            start = 0;
            for(int i = 0; i < 2; ++i) {
                end = line.find(' ', start);
                _PESize[i] = stoi( line.substr(start, end - start) );
                start = end + 1;
            }
        }
        else if(n_line == 3) {
            // Cost parameters
            start = 0;
            for(int i = 0; i < 2; ++i) {
                end = line.find(' ', start);
                _CostCoefficients[i] = stod( line.substr(start, end - start) );
                start = end + 1;
            }
        }
        else {
            // HeatMap
            _HeatMap[x][y][z] = stod(line);
            if( x == _HeatMapSize[0] ) {
                if( y == _HeatMapSize[1] ) {
                    ++z;
                    y = 0;
                }
                else {
                    ++y;
                }
                x = 0;
            }
            else {
                ++x;
            }
        }
        ++n_line;
    }
    return 1;
}

void FP::printHeatMap()
{
    for(int k = 0; k < _HeatMapSize[2] + 1; ++k) {
        for(int j = 0; j < _HeatMapSize[1] + 1; ++j) {
            for(int i = 0; i < _HeatMapSize[0] + 1; ++i) {
                cout << i << " " << j << " " << k << "\t" << _HeatMap[i][j][k] << endl;
            }
        }
    }
}