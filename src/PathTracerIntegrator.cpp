#include "Integrator.h"

#include <iostream>

#include <glm/glm.hpp>

#include "Constants.h"

thread_local std::default_random_engine PathTracerIntegrator::rng = std::default_random_engine((std::random_device())());

PathTracerIntegrator::PathTracerIntegrator()
{
    gen = std::uniform_real_distribution<float>(0.0f, 1.0f);
}

glm::vec3 PathTracerIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction)
{
    return traceRay(origin, direction, 0);
}

glm::vec3 PathTracerIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction, int numBounces)
{
    glm::vec3 outputColor = glm::vec3(0);

    glm::vec3 hitPosition;
    glm::vec3 hitNormal;

    material_t hitMaterial;

    bool hit = _scene->castRay(origin, direction, &hitPosition, &hitNormal, &hitMaterial);

    if (hit)
    {
        if (_scene->nextEventEstimation) {
            // end conditions for NEE
            if (hitMaterial.light) {
                if (numBounces == 0)
                    outputColor += hitMaterial.emission;
                else
                    return glm::vec3(0);
            }
            if (!_scene->russianRoulette && numBounces >= _scene->maxDepth) return glm::vec3(0);
        } else {
            if (hitMaterial.light || (!_scene->russianRoulette && numBounces > _scene->maxDepth)) return hitMaterial.emission;
        }

        if (_scene->nextEventEstimation) {
            outputColor += directLighting(
                hitPosition,
                hitNormal,
                hitMaterial,
                origin);
        }

        outputColor += indirectLighting(
            hitPosition,
            hitNormal,
            hitMaterial,
            origin,
            numBounces + 1);

        outputColor += hitMaterial.emission;
    }

    return outputColor;
}

glm::vec3 PathTracerIntegrator::directLighting(
    glm::vec3 position,
    glm::vec3 normal,
    material_t material,
    glm::vec3 origin)
{
    glm::vec3 outputColor = glm::vec3(0);
    int stratifyGridSize = _scene->stratifyGridSize;
    for (auto light : _scene->quadLights)
    {
        for (int n = 0; n < _scene->lightSamples; n++)
        {
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

            glm::vec3 w_in = glm::normalize(lightPosition - position);
            glm::vec3 w_out = glm::normalize(origin - position);

            glm::vec3 F = brdf(material, w_in, w_out, normal);
            float V = occlusion(position, lightPosition);
            float G = geometry(position, normal, lightPosition, lightNormal);

            outputColor += lightArea * light.intensity * F * V * G / (float)_scene->lightSamples;
        }
    }
    return outputColor;
}

/**
 * Computes the indirect lighting by recursively tracing relflected rays.
 *
 * @param position Where light is being sampled.
 * @param normal The surface normal vector.
 * @param material The object material.
 * @param origin The origin or eye position of the ray which hit the object.
 * @param numBounces The depth of the recursion, used as the end condition when Russian Roulette is disabled.
 */
glm::vec3 PathTracerIntegrator::indirectLighting(
    glm::vec3 position,
    glm::vec3 normal,
    material_t material,
    glm::vec3 origin,
    int numBounces)
{
    glm::vec3 outputColor = glm::vec3(0);
    glm::vec3 w_out = glm::normalize(position - origin);
    for (int i = 0; i < numRaysPerBounce; i++) {
        glm::vec3 w_in = sampleHemisphere(normal);
        glm::vec3 f = brdf(material, w_in, w_out, normal);

        glm::vec3 T = TWO_PI * f * glm::dot(w_in, normal);

        if (_scene->russianRoulette) {
            float p = 1 - glm::min(glm::max(T.x, glm::max(T.y, T.z)), 1.0f);

            if (p > gen(rng)) {
                // kill ray
                continue;
            } else {
                // boost ray
                float boost = 1.0f / (1.0f - p);
                outputColor += boost * T * traceRay(position, w_in, numBounces);
            }
        } else {
            outputColor += T * traceRay(position, w_in, numBounces);
        }

    }
    return outputColor / ((float) numRaysPerBounce);
}

/**
 * The modified Blinn-Phong BRDF function.
 *
 * @param mat           The surface material.
 * @param w_in          The incoming direction from the light.
 * @param w_out         The outgoing direction to the observer.
 * @param surfaceNormal The normal vector of the lit surface.
 */
glm::vec3 PathTracerIntegrator::brdf(
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

float PathTracerIntegrator::geometry(
    glm::vec3 surfacePoint,
    glm::vec3 surfaceNormal,
    glm::vec3 lightPoint,
    glm::vec3 lightNormal)
{
    glm::vec3 lightVector = glm::normalize(lightPoint - surfacePoint);

    float cosThetai = glm::dot(lightVector, surfaceNormal);
    float cosThetal = glm::dot(lightVector, -lightNormal);

    cosThetai = (cosThetai < 0) ? 0 : cosThetai;
    cosThetal = (cosThetal < 0) ? 0 : cosThetal;

    //return cosThetai * cosThetal;
    return cosThetai * cosThetal / glm::pow(glm::length(lightPoint - surfacePoint), 2);
}

float PathTracerIntegrator::occlusion(glm::vec3 origin, glm::vec3 target)
{
    glm::vec3 direction = target - origin;
    bool occluded = _scene->castOcclusionRay(origin, glm::normalize(direction), glm::length(direction));
    if (occluded)
        return 0;
    return 1;
}

glm::vec3 PathTracerIntegrator::sampleHemisphere(glm::vec3 normal)
{
    float epsilon1 = gen(rng);
    float epsilon2 = gen(rng);

    float theta = glm::acos(epsilon1);
    float phi = TWO_PI * epsilon2;

    // a sample over the unit hemisphere
    glm::vec3 s = glm::vec3(glm::cos(phi) * glm::sin(theta), glm::sin(phi) * glm::sin(theta), glm::cos(theta));

    // calculate the new coordinate frame
    glm::vec3 w = glm::normalize(normal);

    // create arbitrariy a vector. it's a secret tool we will need later
    glm::vec3 a = glm::vec3(0, 1, 0);

    // if a isn't a good choice, pick a new one
    if (glm::length(w - a) < 0.001 || glm::length(w + a) < 0.001)
        a = glm::vec3(0, 0, 1);

    glm::vec3 u = glm::normalize(glm::cross(a, w));
    glm::vec3 v = glm::cross(w, u);

    glm::vec3 result = s.x * u + s.y * v + s.z * w;

    return result;
}