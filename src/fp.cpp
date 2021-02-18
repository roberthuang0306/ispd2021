#include "fp.h"


// Tile
// ************************************************************

int Tile::_Dimension = 0;

Tile::Tile(int Id, int* HeatMapPos)
{
    _Id = Id;
    _R = 0;
    _HeatMapPos = new int [_Dimension];
    for(int i = 0; i < _Dimension; ++i) _HeatMapPos[i] = HeatMapPos[i];
}

Tile::Tile(const Tile& t)
{
    _Id = t._Id;
    _R = t._R;
    _HeatMapPos = new int [_Dimension];
    for(int i = 0; i < _Dimension; ++i) _HeatMapPos[i] = t._HeatMapPos[i];
    for(int i = 0; i < 2; ++i) _PEPos[i] = t._PEPos[i];
}

Tile::~Tile()
{
    delete [] _HeatMapPos;
}

void Tile::setDimension(int Dimension)
{
    assert(Dimension == 2 || Dimension == 3);
    _Dimension = Dimension;
}

void Tile::print()
{
    cout << _Id << "\t" << _R << "\t";
    for(int i = 0; i < _Dimension; ++i) {
        cout << _HeatMapPos[i] << " ";
    }
    cout << endl;
}

// ************************************************************





// Adjacents
// ************************************************************

int Adjacents::_Dimension = 0;

Adjacents::Adjacents()
{
    _Adjacents = new int*[FaceNum(_Dimension)];
    for(int i = 0; i < FaceNum(_Dimension); ++i) {
        _Adjacents[i] = new int[MaxNeighborNum(_Dimension)];
    }
    _Size = new int[FaceNum(_Dimension)]();
}

Adjacents::~Adjacents()
{
    for(int i = 0; i < FaceNum(_Dimension); ++i) {
        delete [] _Adjacents[i];
    }
    delete _Adjacents;
    delete [] _Size;
}

int Adjacents::entry(int face, int i)
{
    assert(face < FaceNum(_Dimension));
    assert(i < size(face));
    return _Adjacents[face][i];
}

int Adjacents::size(int face)
{
    assert(face < FaceNum(_Dimension));
    return _Size[face];
}

void Adjacents::push(int id, int face)
{
    assert(face < FaceNum(_Dimension));
    assert(_Size[face] < MaxNeighborNum(_Dimension));
    _Adjacents[face][_Size[face]++] = id;
}

void Adjacents::setDimension(int Dimension)
{
    assert(Dimension == 2 || Dimension == 3);
    _Dimension = Dimension;
}

// ************************************************************





// FP
// ************************************************************

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
    Tile::setDimension(_Dimension);
    int n_sampling_points = SamplingSpaceWidth(0) * SamplingSpaceWidth(1);
    if(_Dimension == 3) n_sampling_points *= SamplingSpaceWidth(2);

    // initialize tiles according to _S
    _Tiles.clear();
    _Tiles.reserve(n_sampling_points);
    int id = 0;
    int* pos = new int[_Dimension];
    for(int z = 0; z < (_Dimension == 2 ? 1 : SamplingSpaceWidth(2)); ++z) {
        if (_Dimension == 3) pos[2] = _S * z;
        for(int y = 0; y < SamplingSpaceWidth(1); ++y) {
            pos[1] = _S * y;
            for(int x = 0; x < SamplingSpaceWidth(0); ++x) {
                pos[0] = _S * x;
                _Tiles.emplace_back(id++, pos);
            }
        }
    }
    delete [] pos;

    // initialize connectivity
    Adjacents::setDimension(_Dimension);
    _Connectivities.clear();
    _Connectivities.resize(n_sampling_points);
    bool left, right, front, back, down, up;
    id = 0;
    for(int z = 0; z < (_Dimension == 2 ? 1 : SamplingSpaceWidth(2)); ++z) {
        if(_Dimension == 3) {
            down = (z != 0);
            up = (z != SamplingSpaceWidth(2) - 1);
        }
        for(int y = 0; y < SamplingSpaceWidth(1); ++y) {
            front = (y != 0);
            back = (y != SamplingSpaceWidth(1) - 1);
            for(int x = 0; x < SamplingSpaceWidth(0); ++x) {
                left = (x != 0);
                right = (x != SamplingSpaceWidth(0) - 1);
                if(left) _Connectivities[id].push(id - 1, 0);
                if(right) _Connectivities[id].push(id + 1, 1);
                if(front) _Connectivities[id].push(id - SamplingSpaceWidth(0), 2);
                if(back) _Connectivities[id].push(id + SamplingSpaceWidth(0), 3);
                if(_Dimension == 3 && down) _Connectivities[id].push(id - SamplingSpaceWidth(0) * SamplingSpaceWidth(1), 4);
                if(_Dimension == 3 && up) _Connectivities[id].push(id + SamplingSpaceWidth(0) * SamplingSpaceWidth(1), 5);
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
    int R = 1;
    while(true)
    {

        ++R;
    }
    
}

double FP::HeatMap_Resolution(int x, int y, int z)
{
    return _HeatMap[x][y][z];
}


// Private Funcs

// For calculating heatmap resolution

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

// For Partition

int FP::SamplingSpaceWidth(int axis)
{
    assert(_S);
    assert(axis < _Dimension);
    int tmp = _S * _HeatMapSize[axis];
    if(tmp % 10) return tmp/10 + 1;
    else return tmp/10;
}

bool FP::isEdge(int tid)
{
    for(int face = 0; face < 2*_Dimension; ++face) {
        if(_Connectivities[tid].size(face) == MaxNeighborNum(_Dimension)) {
            return true;
        }
    }
    return false;
}

// Debug

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

void FP::printTiles()
{
    cout << "Tiles:" << endl;
    for(size_t i = 0; i < _Tiles.size(); ++i)
        _Tiles[i].print();
    cout << endl;
}

void FP::printConnectivities()
{
    cout << "Connectivities:" << endl;
    for(size_t id = 0; id < _Connectivities.size(); ++id) {
        for(int face = 0; face < 2*_Dimension; ++face) {
            for(int i = 0; i < _Connectivities[id].size(face); ++i) {
                cout << id << "\t" << face << "\t"  << _Connectivities[id].entry(face, i)<< endl;
            }
        }
    }
    cout << endl;
}


// Utils

const string& FP::getLineNextToken(const string& line, size_t& start, bool begin)
{
    size_t end = line.find(' ', begin ? 0 : start);
    static string s = line.substr(start, end - start);
    start = end + 1;
    return s;
}

// ************************************************************
