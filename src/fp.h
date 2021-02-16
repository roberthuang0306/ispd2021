#ifndef FP_H
#define FP_H

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class Tile;
class FP;

// Max Dimension
#define D 3


class Tile
{
public:

private:
    int _Id;
    int _R; // 2^(-r) = resolution
    int _Size[D];

    int _PEPos[2];
    int _HeatMapPos[D];
};


class FP
{
public:
    FP();
    ~FP();
    int readHeatMap(string filename);

    // Debug
    void printHeatMap();
private:
    // Input
    int         _Dimension;
    double***   _HeatMap;
    int         _HeatMapSize[D];
    int         _PESize[2];

    // Cost parameters
    double      _CostCoefficients[2];
    
    // Solution
    int         _S; // inverse sampling step
    Tile*       _Tiles;
    int         _n_Tiles;

    // todo
    // Connections
    // Adapters


    // Parameters
    double      _ComputeRatio;
};

#endif