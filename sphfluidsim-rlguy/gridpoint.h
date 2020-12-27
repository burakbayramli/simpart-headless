#ifndef GRIDPOINT_H
#define GRIDPOINT_H

#endif // GRIDPOINT_H

#include <vector>
#include "glm/glm.hpp"
#include <iostream>
using namespace std;


// definition of a point in the spatial grid

struct GridPoint {
    glm::vec3 position;
    int id;
    double tx, ty, tz;
    int i, j, k;
    bool isInGridCell = false;
    bool isMarkedForRemoval = false;
};
