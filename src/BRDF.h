#pragma once

#include <glm/glm.hpp>

#include <random>
#include <set>

#include "Scene.h"

class BRDF
{
protected:
    Scene *_scene;

    // random devices
    static thread_local std::default_random_engine rng;
    std::uniform_real_distribution<float> gen;

public:
    BRDF();

    void setScene(Scene *scene)
    {
        _scene = scene;
    }

    /*************************************
     * COMMANDS REQUIRED FOR BRDF OUTPUT *
     *************************************/

    virtual glm::vec3 brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) = 0;
    virtual glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, float &pdfNormalization) = 0;
    virtual float pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) = 0;
};

class PhongBRDF : public BRDF
{
public:
    glm::vec3 brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);
    glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, float &pdfNormalization);
    float pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);
};

class GGXBRDF : public BRDF
{
public:
    glm::vec3 brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);
    glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, float &pdfNormalization);
    float pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material);
};

class VolumetricBSDF
{
protected:
    Scene *_scene;

    // random devices
    static thread_local std::default_random_engine rng;
    std::uniform_real_distribution<float> gen;

public:
    VolumetricBSDF();

    glm::vec3 importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, std::set<std::string> &volumes, float &pdfNormalization);

    void setScene(Scene *scene)
    {
        _scene = scene;
    }
};