#include "Integrator.h"

#include <iostream>
#include <random>

#include "Constants.h"

thread_local std::default_random_engine MonteCarloDirectIntegrator::rng = std::default_random_engine((std::random_device())());

MonteCarloDirectIntegrator::MonteCarloDirectIntegrator()
{
    gen = std::uniform_real_distribution<float>(0.0f, 1.0f);
}

glm::vec3 MonteCarloDirectIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction)
{
    glm::vec3 outputColor = glm::vec3(0);

    glm::vec3 hitPosition;
    glm::vec3 hitNormal;
    material_t hitMaterial;
    bool hit = _scene->castRay(origin, direction, &hitPosition, &hitNormal, &hitMaterial);

    if (hit)
    {
        if (hitMaterial.light)
            return hitMaterial.emission;

        int stratifyGridSize = _scene->stratifyGridSize;

        for (auto light : _scene->quadLights)
        {

            for (int n = 0; n < _scene->lightSamples; n++)
            {
                // TODO: select random location on light

                float offsetB = gen(rng);
                float offsetC = gen(rng);

                if (_scene->lightStratify)
                {
                    int cellB = n % stratifyGridSize;
                    int cellC = n / stratifyGridSize;

                    offsetB = (cellB + offsetB) / (float)stratifyGridSize;
                    offsetC = (cellC + offsetC) / (float)stratifyGridSize;
                }

                glm::vec3 lightPosition = light.a + offsetB * light.ab + offsetC * light.ac;
                glm::vec3 lightNormal = glm::normalize(glm::cross(light.ac, light.ab));
                float lightArea = glm::length(glm::cross(light.ab, light.ac));

                glm::vec3 w_in = glm::normalize(lightPosition - hitPosition);
                glm::vec3 w_out = glm::normalize(origin - hitPosition);

                glm::vec3 F = brdf(hitMaterial, w_in, w_out, hitNormal);
                float V = occlusion(hitPosition, lightPosition);
                float G = geometry(hitPosition, hitNormal, lightPosition, lightNormal);

                outputColor += lightArea * light.intensity * F * V * G / (float)_scene->lightSamples;
                //outputColor += brdf(hitMaterial) * V;
            }
        }
    }

    return outputColor;
}

float MonteCarloDirectIntegrator::geometry(
    glm::vec3 surfacePoint,
    glm::vec3 surfaceNormal,
    glm::vec3 lightPoint,
    glm::vec3 lightNormal)
{
    //float cosThetai = glm::dot(surfacePoint - origin, surfaceNormal) / glm::length(surfacePoint-origin);
    //float cosThetal = -glm::dot(lightPoint - surfacePoint, lightNormal) / glm::length(lightPoint-surfacePoint);

    glm::vec3 lightVector = glm::normalize(lightPoint - surfacePoint);

    float cosThetai = glm::dot(lightVector, surfaceNormal);
    float cosThetal = glm::dot(lightVector, -lightNormal);

    cosThetai = (cosThetai < 0) ? 0 : cosThetai;
    cosThetal = (cosThetal < 0) ? 0 : cosThetal;

    //return cosThetai * cosThetal;
    return cosThetai * cosThetal / glm::pow(glm::length(lightPoint - surfacePoint), 2);
}

glm::vec3 MonteCarloDirectIntegrator::brdf(
    material_t mat,
    glm::vec3 w_in,
    glm::vec3 w_out,
    glm::vec3 surfaceNormal)
{
    glm::vec3 reflection = 2 * glm::dot(surfaceNormal, w_out) * surfaceNormal - w_out;

    glm::vec3 diffuse = mat.diffuse * INV_PI;
    glm::vec3 specular = mat.specular * (mat.shininess + 2) / (2 * PI) * glm::pow(glm::dot(reflection, w_in), mat.shininess);

    return diffuse + specular;
}

float MonteCarloDirectIntegrator::occlusion(glm::vec3 origin, glm::vec3 target)
{
    glm::vec3 direction = target - origin;
    bool occluded = _scene->castOcclusionRay(origin, glm::normalize(direction), glm::length(direction));
    if (occluded)
        return 0;
    return 1;
}