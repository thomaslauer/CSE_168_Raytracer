#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include <random>

#include "Scene.h"

class Integrator {

protected:

    Scene* _scene;

public:

    void setScene(Scene* scene)
    {
        _scene = scene;
    }

    virtual glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction) = 0;

};


class RayTracerIntegrator : public Integrator {

private:

    glm::vec3 computeShading(
        glm::vec3 incidentDirection,
        glm::vec3 toLight,
        glm::vec3 normal,
        glm::vec3 lightBrightness,
        const material_t& material);

    glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction, int depth);

public:

    virtual glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction);

};

class AnalyticDirectIntegrator : public Integrator {
private:
    float calculateSubtendedAngle(glm::vec3 r, glm::vec3 p1, glm::vec3 p2);
    glm::vec3 calculateSubtendedAngleNormal(glm::vec3 r, glm::vec3 p1, glm::vec3 p2);
public:
    virtual glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction);
};

class MonteCarloDirectIntegrator : public Integrator {
private:
    glm::vec3 brdf(
        material_t mat,
        glm::vec3 w_in,
        glm::vec3 w_out,
        glm::vec3 surfaceNormal);
    float occlusion(glm::vec3 origin, glm::vec3 target);
    float geometry(
        glm::vec3 origin,
        glm::vec3 direction,
        glm::vec3 surfaceNormal,
        glm::vec3 lightDirection,
        glm::vec3 lightNormal);

    static thread_local std::default_random_engine rng;
    std::uniform_real_distribution<float> gen;

public:
    MonteCarloDirectIntegrator();
    virtual glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction);
};