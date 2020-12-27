#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H

#endif // SPHPARTICLE_H

#include <vector>
#include "glm/glm.hpp"

#include <iostream>
using namespace std;

struct SPHParticle {
    glm::vec3 position;
    glm::vec3 prevPosition;   // only used for obstacles
    glm::vec3 velocity;
    glm::vec3 velocityAtHalfTimeStep;  // for leap frog integration
    glm::vec3 acceleration;
    double soundSpeed;
    double mass;
    double density;
    double densityVelocity;
    double pressure;
    std::vector<SPHParticle*> neighbours;
    int gridID;  // used for spatial grid lookup
    bool isHalfTimeStepVelocityInitialized;
    bool isObstacle;
    bool isMarkedForRemoval = false;
    bool isVisible = true;
    double zdistance = 0.0;

    // graphics
    glm::vec3 color;
    double colorDensity;
    double colorVelocity = 0.0;
    bool isStuckInBoundary = false;
    double boundaryAlphaValue = 1.0;
    double alpha = 1.0;
};


