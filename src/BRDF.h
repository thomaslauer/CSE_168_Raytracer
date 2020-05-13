#pragma once

#include <glm/glm.hpp>

#include <random>

#include "Scene.h"

class BRDF
{
protected:

    Scene* _scene;

    glm::vec3 sphereCoordsToVector(float theta, float phi, glm::vec3 samplingSpaceCenter);

    float averageVector(glm::vec3 vec);


    // random devices
    static thread_local std::default_random_engine rng;
    std::uniform_real_distribution<float> gen;

public:

    BRDF();

    void setScene(Scene* scene)
    {
        _scene = scene;
    }


    /*************************************
     * COMMANDS REQUIRED FOR BRDF OUTPUT *
     *************************************/

    virtual glm::vec3 brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) = 0;
    virtual glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material) = 0;
    virtual float pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) = 0;
};

class PhongBRDF : public BRDF {
public:
    glm::vec3 brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);
    glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material);
    float pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);
};

class GGXBRDF : public BRDF {
public:
    glm::vec3 brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);
    glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material);
    float pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);

    float microfacetDistribution(float halfAngle, material_t material);
    float microfacetSelfShadowing(glm::vec3 normal, glm::vec3 view, material_t material);
    glm::vec3 fresnel(glm::vec3 w_in, glm::vec3 halfVector, material_t material);
};