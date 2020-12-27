#ifndef CELLHASH_H
#define CELLHASH_H

#include <unordered_map>
#include <vector>
#include "gridcell.h"
#include <iostream>
using namespace std;


// manages hashing gridcell objects
// cells are stored in a chaied hash
// insert, remove, and retrieve gridcells from the hash

class CellHash
{
public:
    CellHash();
    bool isGridCellInHash(int i, int j, int k);
    void insertGridCell(GridCell *cell);
    void removeGridCell(GridCell *cell);
    GridCell* getGridCell(int i, int j, int k);
    GridCell* findGridCell(int i, int j, int k, bool *isGridCellFound);
    void getGridCells(std::vector<GridCell*> *cells);

private:
    inline long computeHash(int i, int j, int k);

    long maxNumHashValues;
    std::unordered_map<long, std::vector<GridCell*>> cellMap;
};

#endif // CELLHASH_H
