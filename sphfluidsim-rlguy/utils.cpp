#include "utils.h"

namespace utils {
    
    // finds the set of points on a half circle centred at position pos
    // with radius r around axis vector. Circle half is cut off by normal plane.
    // points are stared in *points
    // this method is so convoluted
    void getHalfCirclePoints(glm::vec3 pos, float r,
                                   glm::vec3 axis, glm::vec3 plane,
                                   std::vector<glm::vec3> *points) {

        glm::vec3 v1 = glm::normalize(glm::vec3(-axis.y, axis.x, 0.0));
        glm::vec3 v2 = glm::normalize(glm::cross(axis, v1));

        std::vector<glm::vec3> rawPoints;

        int n = 100;
        for(int i = 0; i < n; i++)
        {
            float theta = 2.0f * 3.1415926 * float(i) / n;
            glm::vec3 p = glm::vec3( r*cos(theta)*v1.x + r*sin(theta)*v2.x,
                                     r*cos(theta)*v1.y + r*sin(theta)*v2.y,
                                     r*cos(theta)*v1.z + r*sin(theta)*v2.z) + pos;
            glm::vec3 op = p - pos;
            if (glm::dot(op, plane) > 0) {
                rawPoints.push_back(p);
            }
        }

        int idx = 0;
        float d = 2*r - 0.1 * r;
        float dsq = d*d;
        if (glm::distance(rawPoints[0], rawPoints[rawPoints.size()-1]) > d) {
            for (int i=0; i<(int)rawPoints.size(); i++) {
                points->push_back(rawPoints[i]);
            }
            return;
        }

        for (int i=0; i<(int)rawPoints.size()-1; i++) {
            glm::vec3 v = rawPoints[i] - rawPoints[i+1];
            if (glm::dot(v, v) > dsq) {
              idx = i;
              break;
            }
        }

        for (int i=idx+1; i<(int)rawPoints.size(); i++) {
            points->push_back(rawPoints[i]);
        }
        for (int i=0; i<=idx; i++) {
            points->push_back(rawPoints[i]);
        }

    }



    // closest distance from point p to line with origin o and direction dir
    float pointToLineDistance(glm::vec3 p, glm::vec3 o, glm::vec3 dir) {
        return glm::length(((o - p) - glm::dot((o - p), dir) * dir));
    }

    // finds closest point on line with origin o and direction dir from point p
    glm::vec3 closestPointOnLineFromPoint(glm::vec3 p, glm::vec3 o,
                                                glm::vec3 dir) {
        return p + ((o - p) - glm::dot((o - p), dir) * dir);
    }

    // linear interpolation with t in [0,1]
    float lerp(float x1, float x2, float t) {
        return x1 + t*(x2 - x1);
    }

    float smoothstep(float t) {
        return t*t*(3 - 2*t);
    }

    // From csc 486 Summer 2014 lecture slides
    /**
     *  For ease in/out motion. Formula from interpolation lecture slides
     *  @param	t  time to calculate position for in range [0,1]
     *  @param	t1 time at which motion stops accelerating in range [0,1] 0 < t1 < t2
     *  @param	t2 time at which motion starts decelerating in range [0,1] t1 < t2 < 1
     *  @return	position of point at time t
     */
    float easeInOut(float t, float t1, float t2) {
        float v0;
        float d;

        v0 = 2/(1+t2-t1); /* constant velocity attained */
        if (t<t1) {
            d = v0*t*t/(2*t1);
        } else {
            d = v0*t1/2;
            if (t<t2) {
                d += (t-t1)*v0;
            } else {
                d += (t2-t1)*v0;
                d += (t-t*t/2-t2+t2*t2/2)*v0/(1-t2);
            }
        }
        return (d);
    }

    // find point on plane intesected by a line
    // lOrigin - point on line
    // lDir - direction a line
    // pOrigin - point on plane
    // pNormal - normal of plane
    glm::vec3 linePlaneIntersection(glm::vec3 lOrigin, glm::vec3 lDir,
                                          glm::vec3 pOrigin, glm::vec3 pNormal) {
        float dot = glm::dot(lDir, pNormal);
        if (dot == 0.0) {
            return glm::vec3(0.0, 0.0, 0.0);
        }
        float d = glm::dot((pOrigin - lOrigin), pNormal) / dot;

        return lOrigin + d * lDir;
    }

    std::vector<glm::vec3> createPointPanel(float width, float height,
                                            float spacing, int numLayers,
                                            glm::vec3 w, glm::vec3 h,
                                            bool isStaggered) {

        // adjust spacing to fit width/height constraints
        int spaces = floor(width/spacing);
        float xpad = width/spaces;
        int nx = spaces+1;

        spaces = floor(height/spacing);
        float ypad = height/spaces;
        int ny = spaces + 1;

        float zpad = spacing;
        int nz = numLayers;

        w = glm::normalize(w);
        h = glm::normalize(h);
        glm::vec3 normal = glm::cross(w, h);

        std::vector<glm::vec3> points;
        glm::vec3 start = (float)(-0.5*width)*w -
                          (float)(0.5*height)*h -
                          (float)(0.5*(nz-1)*zpad)*normal;
        for (int k=0; k<nz; k++) {
            for (int j=0; j<ny; j++) {
                for (int i=0; i<nx; i++) {
                    glm::vec3 p = start + i*xpad*w + j*ypad*h + k*zpad*normal;
                    points.push_back(p);

                    if (isStaggered && i != nx-1 && j != ny-1 && (k != nz-1 || nz == 1)) {
                        p = p + (float)(0.5*xpad)*w +
                                (float)(0.5*ypad)*h +
                                (float)(0.5*zpad)*normal;
                        points.push_back(p);
                    }
                }
            }
        }

        return points;
    }

    std::vector<glm::vec3> rotatePoints(std::vector<glm::vec3> points, Quaternion q) {
        glm::mat4 rot = q.getRotationMatrix();
        std::vector<glm::vec3> newPoints;
        for (uint i=0; i<points.size(); i++) {
            glm::vec4 t = glm::vec4(points[i].x, points[i].y, points[i].z, 1.0);
            t = rot*t;
            newPoints.push_back(glm::vec3(t.x, t.y, t.z));
        }

        return newPoints;
    }


    std::vector<glm::vec3> translatePoints(std::vector<glm::vec3> points, glm::vec3 trans) {
        std::vector<glm::vec3> newPoints;
        for (uint i=0; i<points.size(); i++) {
            newPoints.push_back(points[i] + trans);
        }

        return newPoints;
    }

    std::vector<glm::vec3> mergePoints(std::vector<glm::vec3> points1,
                                       std::vector<glm::vec3> points2) {
        for (uint i=0; i<points2.size(); i++) {
            points1.push_back(points2[i]);
        }

        return points1;
    }

}
















