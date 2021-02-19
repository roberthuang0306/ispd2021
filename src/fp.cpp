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
    // inv_sampling_step
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

double FP::HeatMap_Resolution(double x, double y, double z)
{
    int xc = ceil(x), xf = floor(x), yc = ceil(y), yf = floor(y), zc = ceil(z), zf = floor(z);
    double xl = x-(double)xf, xr = (double)xc-x,
           yl = y-(double)yf, yr = (double)yc-y,
           zl = z-(double)zf, zr = (double)zc-z;
    if( x == (int)x){xl = 1; xr = 0; xc = x; xf = x;}
    if( y == (int)y){yl = 1; yr = 0; yc = y; yf = y;}
    if( z == (int)z){zl = 1; zr = 0; zc = z; zf = z;}
    double interpolate_value = HeatMap_Resolution(xf, yf, zf)*xr*yr*zr +
                               HeatMap_Resolution(xf, yf, zc)*xr*yr*zl +
                               HeatMap_Resolution(xf, yc, zf)*xr*yl*zr +
                               HeatMap_Resolution(xc, yf, zf)*xl*yr*zr +
                               HeatMap_Resolution(xf, yc, zc)*xr*yl*zl +
                               HeatMap_Resolution(xc, yf, zc)*xl*yr*zl +
                               HeatMap_Resolution(xc, yc, zf)*xl*yl*zr +
                               HeatMap_Resolution(xc, yc, zc)*xl*yl*zl;

    return interpolate_value;
}

// Private Funcs

// For calculating heatmap resolution

double FP::Integral_Resolution3()
{
    // TODO
    double integral = 0.0;
    if( _Dimension == 2){
        for( int i = 0; i < _HeatMapSize[0]-1; ++i){
            for ( int j = 0; j < _HeatMapSize[1]-1; ++j){
                // ax+by+cxy+d = r
                // compute r^2
                double r_0 = HeatMap_Resolution(i, j, 0), 
                       r_1 = HeatMap_Resolution(i+1, j, 0), 
                       r_2 = HeatMap_Resolution(i, j+1, 0), 
                       r_3 = HeatMap_Resolution(i+1, j+1, 0);
                double a = r_1-r_0,
                       b = r_2-r_0,
                       c = r_3-r_1-r_2+r_0,
                       d = r_0;

                vector<double> r_Constants = {a, b, c, d};
                vector<int> r_xPowers = {1, 0, 1, 0};
                vector<int> r_yPowers = {0, 1, 1, 0};
                vector<double> r2_Constants;
                vector<double> r2_xPowers;
                vector<double> r2_yPowers;

                for( auto& C1: r_Constants){
                    for( auto& C2: r_Constants){
                        r2_Constants.push_back(C1*C2);
                    }
                }
                for( auto& P1: r_xPowers){
                    for( auto& P2: r_xPowers){
                        r2_xPowers.push_back(P1+P2);
                    }
                }
                for( auto& P1: r_yPowers){
                    for( auto& P2: r_yPowers){
                        r2_yPowers.push_back(P1+P2);
                    }
                }
                assert(r2_Constants.size() == 16);
                assert(r2_xPowers.size() == 16);
                assert(r2_yPowers.size() == 16);

                for( int cnt = 0; cnt < 16; ++cnt){
                    integral += r2_Constants[cnt]* 1/((double)r2_xPowers[cnt]+1)* 1/((double)r2_yPowers[cnt]+1);
                }
            }
        }
    }
    else if( _Dimension == 3){
        for( int i = 0; i < _HeatMapSize[0]-1; ++i){
            for ( int j = 0; j < _HeatMapSize[1]-1; ++j){
                for ( int k = 0; k < _HeatMapSize[2]-1; ++k){
                    // ax+by+cz+dxy+exz+fyz+gxyz+h = r
                    // compute r^3
                    double r_0 = HeatMap_Resolution(i, j, k),
                           r_1 = HeatMap_Resolution(i+1, j, k), 
                           r_2 = HeatMap_Resolution(i, j+1, k), 
                           r_3 = HeatMap_Resolution(i, j, k+1), 
                           r_4 = HeatMap_Resolution(i+1, j+1, k), 
                           r_5 = HeatMap_Resolution(i+1, j, k+1), 
                           r_6 = HeatMap_Resolution(i, j+1, k+1), 
                           r_7 = HeatMap_Resolution(i+1, j+1, k+1);
                    double a = r_1-r_0,
                           b = r_2-r_0,
                           c = r_3-r_0,
                           d = r_4-r_1-r_2+r_0,
                           e = r_5-r_1-r_3+r_0,
                           f = r_6-r_2-r_3+r_0,
                           g = r_7-r_6-r_5-r_4+r_1+r_2+r_3-r_0,
                           h = r_0;


                    vector<double> r_Constants = {a, b, c, d, e, f, g, h};
                    vector<int> r_xPowers = {1, 0, 0, 1, 1, 0, 1, 0};
                    vector<int> r_yPowers = {0, 1, 0, 1, 0, 1, 1, 0};
                    vector<int> r_zPowers = {0, 0, 1, 0, 1, 1, 1, 0};
                    vector<double> r3_Constants;
                    vector<double> r3_xPowers;
                    vector<double> r3_yPowers;
                    vector<double> r3_zPowers;

                    for( auto& C1: r_Constants){
                        for( auto& C2: r_Constants){
                            for( auto& C3: r_Constants){
                                r3_Constants.push_back(C1*C2*C3);
                            }
                        }
                    }
                    for( auto& P1: r_xPowers){
                        for( auto& P2: r_xPowers){
                            for( auto& P3: r_xPowers){
                                r3_xPowers.push_back(P1+P2+P3);
                            }
                        }
                    }
                    for( auto& P1: r_yPowers){
                        for( auto& P2: r_yPowers){
                            for( auto& P3: r_yPowers){
                                r3_yPowers.push_back(P1+P2+P3);
                            }
                        }
                    }
                    for( auto& P1: r_zPowers){
                        for( auto& P2: r_zPowers){
                            for( auto& P3: r_zPowers){
                                r3_zPowers.push_back(P1+P2+P3);
                            }
                        }
                    }
                    assert(r3_Constants.size() == 512);
                    assert(r3_xPowers.size() == 512);
                    assert(r3_yPowers.size() == 512);
                    assert(r3_zPowers.size() == 512);

                    for( int cnt = 0; cnt < 512; ++cnt){
                        integral += r3_Constants[cnt]* 1/((double)r3_xPowers[cnt]+1)* 1/((double)r3_yPowers[cnt]+1)* 1/((double)r3_zPowers[cnt]+1);;
                    }
                }
            }
        }

    }
    return integral;
}

double FP::Integral_Resolution(const int* const origin, int width)
{
    // TODO
    assert(_S);
    double integral = 0.0;
    double *h_origin = new double[_Dimension];
    for( int i = 0; i < _Dimension; ++i) h_origin[i] = ((double)origin[i])/_S;
    double h_width = ((double)width)/_S;

    if( _Dimension == 2){
        for( double i = h_origin[0]; i < h_origin[0]+h_width; ){
            // compute next_i
            double next_i = 0.0;
            if( i == (int)i )   next_i = i+1;
            else    next_i = ceil(i);
            if( next_i > h_origin[0]+h_width)   
                next_i = h_origin[0]+h_width;
            if( next_i > h_origin[0]+h_width) break;

            for( double j = h_origin[1]; j < h_origin[1]+h_width; ){
                // compute next_j
                double next_j = 0.0;
                if( j == (int)j )   next_j = j+1;
                else    next_j = ceil(j);
                if( next_j > h_origin[1]+h_width)
                    next_j = h_origin[1]+h_width;
                if( next_j > h_origin[1]+h_width) break;
                // compute interpolate
                integral += ((HeatMap_Resolution(i, j)+
                              HeatMap_Resolution(i, next_j)+
                              HeatMap_Resolution(next_i, j)+
                              HeatMap_Resolution(next_i, next_j))/4)*(next_i-i)*(next_j-j);

                j = next_j;
            }
            i = next_i;
        }
    }
    else if( _Dimension == 3){
        for( double i = h_origin[0]; i < h_origin[0]+h_width; ){
            // compute next_i
            double next_i = 0.0;
            if( i == (int)i )   next_i = i+1;
            else    next_i = ceil(i);
            if( next_i > h_origin[0]+h_width)   
                next_i = h_origin[0]+h_width;

            for( double j = h_origin[1]; j < h_origin[1]+h_width; ){
                // compute next_j
                double next_j = 0.0;
                if( j == (int)j )   next_j = j+1;
                else    next_j = ceil(j);
                if( next_j > h_origin[1]+h_width)
                    next_j = h_origin[1]+h_width;
                    for( double k = h_origin[2]; k < h_origin[2]+h_width; ){
                        // compute next_k
                        double next_k = 0.0;
                        if( k == (int)k )   next_k = k+1;
                        else    next_k = ceil(k);
                        if( next_k > h_origin[2]+h_width)
                            next_k = h_origin[2]+h_width;
                        // compute interpolate
                        integral += ((HeatMap_Resolution(i, j, k)+
                                      HeatMap_Resolution(i, j, next_k)+
                                      HeatMap_Resolution(i, next_j, k)+
                                      HeatMap_Resolution(i, next_j, next_k)+
                                      HeatMap_Resolution(next_i, j, k)+
                                      HeatMap_Resolution(next_i, j, next_k)+
                                      HeatMap_Resolution(next_i, next_j, k)+
                                      HeatMap_Resolution(next_i, next_j, next_k))/8)*(next_i-i)*(next_j-j)*(next_k-k);

                        k = next_k;
                    }
                j = next_j;
            }
            i = next_i;
        }
    }

    delete[] h_origin;
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
