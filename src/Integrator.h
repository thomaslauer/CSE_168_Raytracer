#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include <random>

#include "Scene.h"
#include "BRDF.h"

class Integrator {

protected:

    Scene* _scene;

public:

    virtual void setScene(Scene* scene)
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


class PathTracerIntegrator : public Integrator {
private:
    // random number generation
    static thread_local std::default_random_engine rng;
    std::uniform_real_distribution<float> gen;

    int numRaysPerBounce = 1;

    BRDF* _phongBRDF;
    BRDF* _ggxBRDF;

    glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction, int numBounces);

    // separate functions for direct and indirect lighting components

    // direct lighting from light sources
    glm::vec3 directLighting(
        glm::vec3   position,
        glm::vec3   normal,
        material_t  material,
        glm::vec3   origin);

    // indirect light calculated as the integral over the unit hemisphere
    glm::vec3 indirectLighting(
        glm::vec3   position,
        glm::vec3   normal,
        material_t  material,
        glm::vec3   origin,
        int         numBounces);


    glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, float& pdfNormalization);

    // calculate BRDF for a material
    glm::vec3 brdf(
        material_t mat,
        glm::vec3 w_in,
        glm::vec3 w_out,
        glm::vec3 surfaceNormal);

    float geometry(
        glm::vec3 direction,
        glm::vec3 surfaceNormal,
        glm::vec3 lightDirection,
        glm::vec3 lightNormal);

    float occlusion(glm::vec3 origin, glm::vec3 target);

    glm::vec3 ggxBRDF(
        material_t mat,
        glm::vec3 w_in,
        glm::vec3 w_out,
        glm::vec3 surfaceNormal);

    float ggxMicrofacetDistribution(
        material_t mat,
        float halfAngle);

    float ggxMicrofacetSelfShadowing(
        material_t mat,
        glm::vec3 normal,
        glm::vec3 view);

    glm::vec3 ggxFresnel(
        material_t mat,
        glm::vec3 w_in,
        glm::vec3 halfVector);

public:
    PathTracerIntegrator();
    void setScene(Scene* scene);
    virtual glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction);
};