#include "fp.h"

int Tile::_Dimension = 0;

FP::FP()
{
    _Dimension = 0;
    _HeatMap = 0;
    _HeatMapSize = 0;
    for(int i = 0; i < 2; ++i) _PESize[i] = 0;
    for(int i = 0; i < 2; ++i) _CostCoefficients[i] = 0.0;
    _S = 0;
    _ComputeRatio = 0.0;
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
    if( _HeatMapSize ) {
        delete [] _HeatMapSize;
    }
}


// Public Funcs

int FP::readHeatMap(const string& filename)
{
    ifstream file(filename);
    if(!file.is_open()){
        cout << "Unable to open file \"" << filename << "\"\n";
        return 0;
    }
    string line;
    int n_line = 0;
    size_t start;
    int x = 0, y = 0, z = 0;
    while(getline(file, line)) {
        if(line[0] == '#') continue;

        if(n_line == 0) {
            // Dimension
            _Dimension = stoi(line);
            assert(_Dimension == 2 || _Dimension == 3);
        }
        else if(n_line == 1) {
            // HeatMap size
            _HeatMapSize = new int[_Dimension];
            for(int i = 0; i < _Dimension; ++i) {
                _HeatMapSize[i] = stoi( getLineNextToken(line, start, !i) );
            }
            _HeatMap = new double**[_HeatMapSize[0] + 1];
            for(int i = 0; i < _HeatMapSize[0] + 1; ++i) {
                _HeatMap[i] = new double*[_HeatMapSize[1] + 1];
                int tmp = (_Dimension == 2 ? 1 : _HeatMapSize[2] + 1);
                for(int j = 0; j < _HeatMapSize[1] + 1; ++j) {
                    _HeatMap[i][j] = new double[tmp];
                }
            }
        }
        else if(n_line == 2) {
            // PE size
            for(int i = 0; i < 2; ++i) {
                _PESize[i] = stoi( getLineNextToken(line, start, !i) );
            }
        }
        else if(n_line == 3) {
            // Cost parameters
            for(int i = 0; i < 2; ++i) {
                _CostCoefficients[i] = stod( getLineNextToken(line, start, !i) );
            }
        }
        else {
            // HeatMap
            _HeatMap[x][y][z] = stod(line);
            if(x == _HeatMapSize[0]) {
                if(y == _HeatMapSize[1]) {
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

void FP::Calculate_S()
{
    // TODO
    _S = 0;
}

void FP::Init_Partition()
{
    assert(_S%10 == 0);
    Tile::setDimension(_Dimension);
    _Tiles.clear();
    int id;
    int* pos = new int[_Dimension];
    for(int z = 0; z < (_Dimension == 2 ? 1 : _S * _HeatMapSize[2]); ++z) {
        if (_Dimension == 3) pos[2] = z;
        for(int y = 0; y < _S * _HeatMapSize[1]; ++y) {
            pos[1] = y;
            for(int x = 0; x < _S * _HeatMapSize[0]; ++x) {
                pos[0] = x;
                _Tiles.emplace_back(id++, pos);
            }
        }
    }
    delete [] pos;

    _Connectivities.clear();
    int n_sampling_points = _S * _HeatMapSize[0] * _S * _HeatMapSize[1];
    if(_Dimension == 3) n_sampling_points *= (_S * _HeatMapSize[2]);
    for(int i = 0; i < n_sampling_points; ++i) _Connectivities.emplace_back();
    bool left, right, front, back, down, up;
    id = 0;
    for(int z = 0; z < (_Dimension == 2 ? 1 : _S * _HeatMapSize[2]); ++z) {
        down = (z != 0);
        up = (z != _S * _HeatMapSize[2] - 1);
        for(int y = 0; y < _S * _HeatMapSize[1]; ++y) {
            front = (y != 0);
            back = (y != _S * _HeatMapSize[1] - 1);
            for(int x = 0; x < _S * _HeatMapSize[0]; ++x) {
                left = (x != 0);
                right = (x != _S * _HeatMapSize[0] - 1);
                if(left) _Connectivities[id].push_back(id - 1);
                if(right) _Connectivities[id].push_back(id + 1);
                if(front) _Connectivities[id].push_back(id - _S * _HeatMapSize[0]);
                if(back) _Connectivities[id].push_back(id + _S * _HeatMapSize[0]);
                if(_Dimension == 3 && down) _Connectivities[id].push_back(id - _S * _HeatMapSize[0] * _S * _HeatMapSize[1]);
                if(_Dimension == 3 && up) _Connectivities[id].push_back(id + _S * _HeatMapSize[0] * _S * _HeatMapSize[1]);
                ++id;
            }
        }
    }
}

void FP::setComputeRatio(const double& r)
{
    _ComputeRatio = r;
}

void FP::Partition()
{

}

inline double FP::HeatMap_Resolution(int x, int y, int z)
{
    return _HeatMap[x][y][z];
}


// Private Funcs

void FP::printHeatMap()
{
    for(int k = 0; k < (_Dimension == 2 ? 1 : _HeatMapSize[2] + 1); ++k) {
        for(int j = 0; j < _HeatMapSize[1] + 1; ++j) {
            for(int i = 0; i < _HeatMapSize[0] + 1; ++i) {
                cout << i << " " << j;
                if(_Dimension == 3) cout << " " << k;
                cout << "    " << _HeatMap[i][j][k] << endl;
            }
        }
    }
}

double FP::Integral_Resolution3()
{
    // TODO
    return 0.0;
}

double FP::Integral_Resolution(int* origin, int* width)
{
    // TODO
    assert(_S);
    return 0.0;
}


// Utils

const string& FP::getLineNextToken(const string& line, size_t& start, bool begin)
{
    size_t end = line.find(' ', begin ? 0 : start);
    static string s = line.substr(start, end - start);
    start = end + 1;
    return s;
}
