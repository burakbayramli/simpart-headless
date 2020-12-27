// How to compile, and other info
// https://burakbayramli.github.io/dersblog/sk/2020/08/sph.html
// g++ simsph2.cpp -std=c++1z  -g -O2 -o /tmp/a.exe; /tmp/a.exe
// 
#include <map>
#include <iostream>
#include <fstream> 
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

const static Vector3d G(0.f, 0.f, 12000*-9.8f); 
const static float REST_DENS = 100.f; // rest density
const static float GAS_CONST = 200.f; // const for equation of state
const static float H = 16.f; 
const static float DIST = 10.f;
const static float HSQ = H*H; 
const static float MASS = 100.f;
const static float VISC = 200.f; 
const static float DT = 0.01f; 
static int LOOP = 40;
static float MAX_COORD = 500;
ofstream foutx;
ofstream fouty;
ofstream foutz;

// puruzsuzlestirici cekirdek ve turevleri
const static float POLY6 = 315.f/(65.f*M_PI*pow(H, 9.f));
const static float SPIKY_GRAD = -45.f/(M_PI*pow(H, 6.f));
const static float VISC_LAP = 45.f/(M_PI*pow(H, 6.f));

const static float EPS = H;
const static float BOUND_DAMPING = -0.5f;

struct Particle {
    Particle(float _x, float _y, float _z) : x(_x, _y, _z),
                                   v(0.f, 0.f, 0.f),
                                   f(0.f, 0.f, 0.f),
                                   rho(0.f),
				   p(0.f) { }
    Vector3d x, v, f;
    float rho, p;
    
};


static vector<Particle> particles;

void InitSPH(void)
{

    foutx.open ("/tmp/simsph-x.csv");
    fouty.open ("/tmp/simsph-y.csv");
    foutz.open ("/tmp/simsph-z.csv");
    
    int balls = 0;
    for(float x = 0.f; x < 50.f; x += 5.f) { 
        for(float y = 450.f; y < 500.f; y += 5.f) {     
            for(float z = 0.f; z < 50.f; z += 5.f) {
                Particle p(x,y,z);
                particles.push_back(p);
                balls ++;
            }
        }
    }
    
    std::cout << "balls:" << balls << std::endl;
}


void Integrate(void)
{
    for(auto &p : particles)
    {
        // ileri Euler entegrasyonu
        if (p.rho > 0.0f) p.v += DT*p.f/p.rho;
        p.x += DT*p.v;

        // sinir sartlarini kontrol et
        if(p.x(0)-EPS < 0.0f)
        {
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = 0.0f;
        }
        if(p.x(0)+EPS > 500.0f) 
        {
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = 500.0f-EPS;
        }
        
        if(p.x(1)-EPS < 0.0f)
        {
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = 0.0f;
        }
        if(p.x(1)+EPS > 500.0f)
        {
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = 500.0f-EPS;
        }

        if(p.x(2)-EPS < 0.0f)
        {
            p.v(2) *= BOUND_DAMPING;
            p.x(2) = 0.0f;
        }
        if(p.x(2)+EPS > 500.0f)
        {
            p.v(2) *= BOUND_DAMPING;
            p.x(2) = 500.0f-EPS;
        }
        
    }
}

void ComputeDensityPressure(void)
{
    for(auto &pi : particles)
    {
        pi.rho = 0.f;
        for(auto &pj : particles)
        {
            Vector3d rij = pj.x - pi.x;
            float r2 = rij.squaredNorm();

            if(r2 < HSQ)
            {
                pi.rho += MASS*POLY6*pow(HSQ-r2, 3.f);
            }
        }
        pi.p = GAS_CONST*(pi.rho - REST_DENS);
    }
}

void ComputeForces(void)
{
    for(auto &pi : particles)
    {
        Vector3d fpress(0.f, 0.f, 0.f);
        Vector3d fvisc(0.f, 0.f, 0.f);
        for(auto &pj : particles)
        {
            if(&pi == &pj)
                continue;

            Vector3d rij = pj.x - pi.x;
            float r = rij.norm();

            if(r < DIST)
            {
                fpress += -rij.normalized()*MASS*(pi.p + pj.p)/(2.f * pj.rho)
                    * SPIKY_GRAD*pow(H-r,2.f);
                fvisc += VISC*MASS*(pj.v - pi.v)/pj.rho * VISC_LAP*(H-r);
            }
        }
        Vector3d fgrav = G * pi.rho;
        pi.f = fpress + fvisc + fgrav;
    }
}

void Update(void)
{
    std::cout << "---" << std::endl;
    ComputeDensityPressure();
    ComputeForces();
    Integrate();
    
    for(auto &pi : particles)
    {
        foutx << pi.x[0] << ";";
        fouty << pi.x[1] << ";";
        foutz << pi.x[2] << ";";
    }
    foutx << '\n';
    fouty << '\n';
    foutz << '\n';
}


int main(int argc, char** argv)
{
    InitSPH();

    for (int i=0;i<LOOP;i++) {
        Update();       
    }
    
    return 0;
}
