#ifndef SPHFLUIDSIMULATION_H
#define SPHFLUIDSIMULATION_H

#include <vector>
#include <time.h>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "glm/glm.hpp"
#include "spatialgrid.h"
#include "sphparticle.h"
#include "sphobstacle.h"
#include "quaternion.h"
#include <iostream>
#include <fstream> 
using namespace std;

class SPHFluidSimulation
{
public:
    SPHFluidSimulation();
    SPHFluidSimulation(double smoothingRadius);
    void update(float dt);
    void compute();

    // user functions
    void addFluidParticles(std::vector<glm::vec3> points);
    void addFluidParticle(glm::vec3 pos, glm::vec3 velocity);
    int addObstacleParticles(std::vector<glm::vec3> points);
    void removeObstacle(int id);
    void setBounds(double xmin, double xmax,
                   double ymin, double ymax,
                   double zmin, double zmax);
    void setDampingConstant(double c);
    std::vector<SPHParticle*> getFluidParticles();
    std::vector<SPHParticle*> getObstacleParticles();
    std::vector<SPHParticle*> getAllParticles();
    float getParticleSize();    // smoothing length
    float getInitialDensity();
    void setObstaclePosition(int id, glm::vec3 pos);
    void translateObstacle(int id, glm::vec3 trans);
    void rotateObstacle(int id, Quaternion q);

    ofstream foutx;
    ofstream fouty;
    ofstream foutz;
    
private:

    // init
    void initSimulationConstants();
    void initKernelConstants();

    SPHParticle* createSPHParticle(glm::vec3 pos, glm::vec3 velocity);
    SPHParticle* createSPHObstacleParticle(glm::vec3 pos);
    SPHParticle* addObstacleParticle(glm::vec3 pos);
    int getUniqueObstacleID();
    int currentObstacleID = 0;

    // simulation
    inline double evaluateSpeedOfSound(SPHParticle *sp);
    inline double evaluateSpeedOfSoundSquared(SPHParticle *sp);
    void removeSPHParticlesMarkedForRemoval();

    // if activated, lines fluid simulation boundaries with obstacle particles
    void initializeBoundaryParticles();
    void updateFluidConstants();               // updated values pulled from
                                               // Lua script
    void updateObstacleVelocity(double dt);    // calc velocity based on new
                                               // positions from translation/rotation
    // notify grid of new positions for fluid/obstacle particles
    void updateGrid();

    // retrieve nearest neighbours for all particles
    void updateNearestNeighbours();
    void updateFluidDensityAndPressure();

    // acceleration on a particle due to being close to boundary
    glm::vec3 calculateBoundaryAcceleration(SPHParticle *sp);
    void updateFluidAcceleration();       //acceleration due to density/pressure,
                                          // viscostiy

    // using leap frog integration
    void updateFluidPosition(double dt);

    // make sure particle has not crossed fluid simulation bounds
    void enforceFluidParticlePositionBounds(SPHParticle *p);

    // don't enforce if bounds has changed since last frame or particles may stick
    // to the boundaries
    bool isEnforcingFluidParticlePositionBoundsThisTimeStep = false;

    // make sure time step is small engoubh so no particle moves more than one
    // smoothing length in a frame
    double calculateTimeStep();


    glm::vec3 calculateFluidParticleColor(SPHParticle *sp);
    bool isFluidParticleStuckToBoundary(SPHParticle *sp);
    std::vector<std::array<double, 3>> fluidGradient;

    // for smoothing color changes
    double maxColorVelocity = 1.0;
    double maxColorAcceleration = 1.0;
    double minColorDensity = 0.0;
    double maxColorDensity = 100.0;
    double colorArrivalRadius = 0.5;
    double stuckToBoundaryRadius = 0.01;
    double stuckToBoundaryAlphaVelocity = 1.0;

    // simulation constants. Controlled by lua script
    double h;                              // smoothing radius
    double hsq;                            // radius squared
    glm::vec3 gravityForce;
    double courantSafetyFactor = 1.0;
    double minTimeStep = 1.0/240.0;
    double gravityMagnitude;
    double initialDensity;
    double pressureCoefficient;
    double particleMass;
    bool isMotionDampingEnabled = false;
    bool displayConsoleOutput = false;
    double motionDampingCoefficient;
    double boundaryDampingCoefficient;
    double ratioOfSpecificHeats;
    double viscosityCoefficient;
    double maximumVelocity;
    double maximumAcceleration;


    // kernel constants
    double poly6Coefficient;
    double spikeyGradCoefficient;
    double viscosityLaplacianCoefficient;

    // boundary constraints
    double boundaryForceRadius = 0.1;
    double minBoundaryForce = 0.0;
    double maxBoundaryForce = 0.0;
    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;
    double zmin = 0.0;
    double zmax = 1.0;
    int boundaryObstacleID;
    bool isBoundaryObstacleInitialized = false;
    bool isHiddenBoundaryParticlesEnabled = true;
    bool isBoundaryParticlesEnabled = false;

    SpatialGrid grid;
    std::vector<SPHParticle*> fluidParticles;
    std::vector<SPHParticle*> obstacleParticles;
    std::vector<SPHParticle*> allParticles;
    std::vector<SPHObstacle*> obstacles;
    std::unordered_map<int,SPHParticle*> particlesByGridID;
    std::unordered_map<int,SPHObstacle*> obstaclesByID;
    bool isSPHParticleRemoved = false;
    glm::vec3 cameraPosition;

    // timing metrics
    double neighbourSearchTime = 0.0;
    double simulationTime = 0.0;
    double graphicsUpdateTime = 0.0;
    double graphicsDrawTime = 0.0;

};

#endif // SPHFLUIDSIMULATION_H
