#include "Integrator.h"

#include <iostream>
#include <set>

#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>

#include "BRDF.h"
#include "Constants.h"
#include "MathUtils.h"

#define EPSILON 0.0001f

thread_local std::default_random_engine VolumetricPathTracerIntegrator::rng = std::default_random_engine((std::random_device())());

VolumetricPathTracerIntegrator::VolumetricPathTracerIntegrator()
{
    gen = std::uniform_real_distribution<float>(0.0f, 1.0f);

    _phongBRDF = new PhongBRDF();
    _ggxBRDF = new GGXBRDF();
    _volumetricBSDF = new VolumetricBSDF();
}

glm::vec3 VolumetricPathTracerIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction)
{
    std::set<std::string> volumes;
    volumes.insert(_scene->defaultVolume);
    return traceRay(origin, direction, volumes, 0);
}

glm::vec3 VolumetricPathTracerIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction, std::set<std::string> &volumes, int numBounces)
{
    glm::vec3 hitPosition;
    glm::vec3 hitNormal;

    material_t hitMaterial;

    bool hit = _scene->castRay(origin, direction, &hitPosition, &hitNormal, &hitMaterial);

    // find the highest priority volume for the current ray
    volume_t volume = highestPriorityVolume(_scene, volumes);

    float t = -glm::log(gen(rng)) * volume.meanScatterDistance;

    float backfaceDistance = glm::length(hitPosition - origin);

    if (hit && (volume.meanScatterDistance == -1 || t > backfaceDistance))
    {
        // calculate surface lighting
        if (hitMaterial.light || numBounces > _scene->maxDepth)
        {
            if (glm::dot(hitNormal, direction) < 0)
                return attenuate(hitMaterial.emission * absdot(hitNormal, direction), backfaceDistance, volume);
            else
                return glm::vec3(0);
        }
        else
            return attenuate(indirectLighting(
                                 hitPosition,
                                 hitNormal,
                                 hitMaterial,
                                 origin,
                                 volumes,
                                 numBounces + 1),
                             backfaceDistance, volume);
    }

    if (numBounces > _scene->maxDepth)
        return glm::vec3(0);

    // move hitPosition to the sampled location
    hitPosition = origin + t * direction;

    glm::vec3 newDirection;
    float pdfNormalization;

    if (volume.scatterDirectionality != 0)
    {
        newDirection = hgScatter(direction, volume, pdfNormalization);
        //std::cout << pdfNormalization << std::endl;
    }
    else
    {
        float theta = glm::acos(2 * gen(rng) - 1);
        float phi = TWO_PI * gen(rng);

        newDirection = sphereCoordsToVector(theta, phi, direction);
        pdfNormalization = 1;
    }

    glm::vec3 scatteredLight = traceRay(hitPosition, newDirection, volumes, numBounces + 1);

    return attenuate(scatteredLight, t, volume) * pdfNormalization;
}

glm::vec3 VolumetricPathTracerIntegrator::indirectLighting(
    glm::vec3 position,
    glm::vec3 normal,
    material_t material,
    glm::vec3 origin,
    std::set<std::string> &volumes,
    int numBounces)
{
    glm::vec3 outputColor = glm::vec3(0);
    glm::vec3 w_out = glm::normalize(origin - position);

    for (int i = 0; i < numRaysPerBounce; i++)
    {
        float pdfNormalization = 1;
        glm::vec3 w_in = importanceSample(normal, w_out, material, volumes, pdfNormalization);

        glm::vec3 f = brdf(normal, w_in, w_out, material);
        glm::vec3 T = f * absdot(w_in, normal) / pdfNormalization;

        //if(material.brdf == GGX_VOLUMETRIC) std::cout << glm::to_string(f) << glm::to_string(T) << std::endl;

        if (_scene->russianRoulette)
        {
            float p = 1 - glm::min(glm::max(T.x, glm::max(T.y, T.z)), 1.0f);

            if (p > gen(rng))
            {
                // kill ray
                continue;
            }
            else
            {
                // boost ray
                float boost = 1.0f / (1.0f - p);
                outputColor += boost * T * traceRay(position + w_in * EPSILON, w_in, volumes, numBounces);
            }
        }
        else
        {
            outputColor += T * traceRay(position + w_in * EPSILON, w_in, volumes, numBounces);
        }
    }
    return outputColor / ((float)numRaysPerBounce);
}

inline glm::vec3 VolumetricPathTracerIntegrator::brdf(
    glm::vec3 surfaceNormal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t mat)
{
    if (mat.brdf == GGX)
        return _ggxBRDF->brdf(surfaceNormal, w_in, w_out, mat);
    else if (mat.brdf == GGX_VOLUMETRIC)
        return glm::vec3(1);
    else
        return _phongBRDF->brdf(surfaceNormal, w_in, w_out, mat);
}

glm::vec3 VolumetricPathTracerIntegrator::importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, std::set<std::string> &volumes, float &pdfNormalization)
{
    glm::vec3 w_in;

    if (material.brdf == GGX)
        w_in = _ggxBRDF->importanceSample(normal, w_out, material, pdfNormalization);
    else if (material.brdf == GGX_VOLUMETRIC)
        w_in = _volumetricBSDF->importanceSample(normal, w_out, material, volumes, pdfNormalization);
    else
        w_in = _phongBRDF->importanceSample(normal, w_out, material, pdfNormalization);

    return w_in;
}

glm::vec3 VolumetricPathTracerIntegrator::hgScatter(
    glm::vec3 direction,
    volume_t volume,
    float &normalization)
{
    float g = volume.scatterDirectionality;
    float g2 = g * g;

    normalization = 1;
    if (g == 1)
        return direction;
    if (g == -1)
        return -direction;

    float e1 = gen(rng);
    float e2 = gen(rng);

    float s = 2 * e1 - 1;
    float mu = (1 / (2 * g)) * (1 + g2 - glm::pow((1 - g2) / (1 + g * s), 2));
    float theta = glm::acos(mu);

    float phi = 2 * PI * e2;

    glm::vec3 out = sphereCoordsToVector(theta, phi, direction);


    normalization = 1 / FOUR_PI * (1 - g2) / glm::pow(1 + g2 - 2 * g * mu, 3/2);

    return out;
}

void VolumetricPathTracerIntegrator::setScene(Scene *scene)
{
    _scene = scene;

    _phongBRDF->setScene(scene);
    _ggxBRDF->setScene(scene);
    _volumetricBSDF->setScene(scene);
}