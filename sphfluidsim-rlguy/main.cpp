#include <GL/glu.h>
#include <vector>
#include <tuple>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <cmath>
#include "glm/glm.hpp"
#include "utils.h"
#include "spatialgrid.h"
#include "sphfluidsimulation.h"

float simulationFPS = 30.0;
SPHFluidSimulation fluidSim;

int main(int argc, char *argv[])
{
    int simulationFPS = 20;

    double radius = 0.2f;
    fluidSim = SPHFluidSimulation(radius);

    double minx = 2.50f;
    double maxx = 2.50f;
    double miny = -2.25f;
    double maxy = 5.0f;
    double minz = -2.50f;
    double maxz = 2.50f;
    double rx = 0.90;
    double ry = 0.2;
    double rz = 1.0;


    fluidSim.setBounds(minx, maxx, miny, maxy, minz, maxz);
    std::vector<glm::vec3> points;
    double n = 1000.f;
    for (int i=0; i<n; i++) {
        float x = minx + rx*((float)rand()/RAND_MAX) * (maxx - minx);
        float y = miny + ry*((float)rand()/RAND_MAX) * (maxy - miny);
        float z = minz + rz*((float)rand()/RAND_MAX) * (maxz - minz);
        glm::vec3 p = glm::vec3(x, y, z);
        points.push_back(p);
    }
    fluidSim.addFluidParticles(points);

    double damp = 2.0f;
    fluidSim.setDampingConstant(damp);


    // create obstacles

    float len = 4.0;
    float height = 4.5;

    float r = 0.6*radius;
    int layers = 1;
    bool isStaggered = true;


    glm::vec3 p = glm::vec3((maxx-minx)/2, 1.0, (maxz-minz)/2);

    std::vector<glm::vec3> o1, o2;
    o1 = utils::createPointPanel(len, height, r, layers,
                                 glm::vec3(1.0,0.0,0.0), glm::vec3(0.0, 0.0, 1.0),
                                 isStaggered);
    o2 = utils::createPointPanel(len, height, r, layers,
                                 glm::vec3(1.0,0.0,0.0), glm::vec3(0.0, 1.0, 0.0),
                                 isStaggered);

    int id1 = fluidSim.addObstacleParticles(o1);
    int id2 = fluidSim.addObstacleParticles(o2);

    //fluidSim.compute();

    float dt = 1.0 / simulationFPS;
    for (int i=0;i<50;i++){	
        float speed = 2.5;
        fluidSim.rotateObstacle(0, Quaternion(glm::vec3(1.0, 0.0, 0.0), speed*dt));
        fluidSim.rotateObstacle(1, Quaternion(glm::vec3(1.0, 0.0, 0.0), speed*dt));
        fluidSim.update(dt);
        fluidSim.compute();
	std::cout << "---------------------" << std::endl;
    }

    
}
