#ifndef FP_H
#define FP_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

using namespace std;

class Tile;
class FP;

// Max Dimension
typedef vector< vector<int> > adjacent_list;


class Tile
{
public:
    Tile(int Id, int* HeatMapPos) {
        _Id = Id;
        _R = 0;
        _HeatMapPos = new int [_Dimension];
        for(int i = 0; i < _Dimension; ++i) _HeatMapPos[i] = HeatMapPos[i];
    }

    ~Tile() {
        delete [] _HeatMapPos;
    }

    static void setDimension(int Dimension) {
        _Dimension = Dimension;
    }

private:
    int _Id;
    int _R;                 // 2^(-R) = Resolution

    int* _HeatMapPos;       // x10 -> sampling space
    int _PEPos[2];

    static int _Dimension;
};

class FP
{
public:
    FP();
    ~FP();
    int             readHeatMap(const string& filename);            // Read an input HeatMap
    void            Calculate_S();                                  // Calculate inverse sampling step
    void            Init_Partition();                               // Initialize Tiles and Connectivities
    void            setComputeRatio(const double& r);
    void            Partition();                                    // Find a legal partition of HeatMap
    void            Place();                                        // Place all the tiles and adpaters on PE array

    inline double   HeatMap_Resolution(int x, int y, int z=0);

    // Debug
    void            printHeatMap();
    
private:
    // Input
    int             _Dimension;                                     // Dimension of the HeatMap
    double***       _HeatMap;
    int*            _HeatMapSize;
    int             _PESize[2];

    // Cost parameters
    double          _CostCoefficients[2];                           // Alpha and Beta
    
    // Solution
    int             _S;                                             // x10 -> inverse sampling step
    vector<Tile>    _Tiles;

    adjacent_list   _Connectivities;

    // todo
    // Adapters


    // Parameters
    double          _ComputeRatio;


    // Funcs
    double          Integral_Resolution3();                         // Integrate Resolution^3 over the entire heatmap
    double          Integral_Resolution(int* origin, int* width);   // Integrate Resolution over some area in sampling space


    // Utils
    const string&   getLineNextToken(const string& line, size_t& start, bool begin=false);
};

#endif