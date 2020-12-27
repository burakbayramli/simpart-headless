#ifndef SPATIALGRID_H
#define SPATIALGRID_H

#include <vector>
#include <unordered_map>
#include <cmath>
#include <time.h>
#include "glm/glm.hpp"
#include "cellhash.h"
#include "utils.h"
#include "gridcell.h"

#include <iostream>
using namespace std;

// Dynamic spatial grid
// Grid cells only exist where a gridpoint exists
// After point insertion, update position/remove using returned reference id

class SpatialGrid
{
public:
    SpatialGrid();
    SpatialGrid(double cell_size);
    int insertPoint(glm::vec3 point);
    void movePoint(int id, glm::vec3 position);
    void removePoint(int id);
    std::vector<glm::vec3> getObjectsInRadiusOfPoint(int ref, double radius);
    std::vector<int> getIDsInRadiusOfPoint(int ref, double radius);
    void update();
    void draw();

private:
    int generateUniqueGridPointID();
    void initializeFreeCells();
    void insertGridPointIntoGrid(GridPoint *p);
    void positionToIJK(glm::vec3 p, int *i, int *j, int *k);
    glm::vec3 IJKToPosition(int i, int j, int k);
    GridCell* getNewGridCell(int i, int j, int k);
    void updateGridPointCellOffset(GridPoint *gp, int i, int j, int k);
    std::vector<int> fastIDNeighbourSearch(int ref, double r, GridPoint *gp);
    void removeGridPointsMarkedForRemoval();

    double size;
    int currentGridPointID;
    std::vector<GridPoint*> points;
    std::unordered_map<int,GridPoint*> gridPointsByID;
    std::vector<GridCell*> freeCells;
    int numInitialFreeCells;
    CellHash cellHashTable;
    bool isCellRemoved = false;
};

#endif // SPATIALGRID_H













