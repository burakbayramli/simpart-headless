#ifndef SPHOBSTACLE_H
#define SPHOBSTACLE_H

#endif // SPHOBSTACLE_H

#include <vector>
#include "glm/glm.hpp"

// obstacle is a collection of fluid particles.

struct SPHObstacle {
    glm::vec3 position;
    std::vector<SPHParticle*> particles;
    bool isVisible = true;
    int id;
};
