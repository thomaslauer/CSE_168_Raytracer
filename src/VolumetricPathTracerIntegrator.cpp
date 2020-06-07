#include "Integrator.h"

#include <iostream>

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
    return traceRay(origin, direction, 0);
}

glm::vec3 VolumetricPathTracerIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction, int numBounces)
{
    glm::vec3 outputColor = glm::vec3(0);

    glm::vec3 hitPosition;
    glm::vec3 hitNormal;

    material_t hitMaterial;

    bool hit = _scene->castRay(origin, direction, &hitPosition, &hitNormal, &hitMaterial);

    if (hit)
    {
        if (hitMaterial.light || numBounces > _scene->maxDepth)
            return hitMaterial.emission;

        outputColor += indirectLighting(
            hitPosition,
            hitNormal,
            hitMaterial,
            origin,
            numBounces + 1);
    }

    return outputColor;
}

glm::vec3 VolumetricPathTracerIntegrator::indirectLighting(
    glm::vec3 position,
    glm::vec3 normal,
    material_t material,
    glm::vec3 origin,
    int numBounces)
{
    glm::vec3 outputColor = glm::vec3(0);
    glm::vec3 w_out = glm::normalize(origin - position);

    for (int i = 0; i < numRaysPerBounce; i++)
    {
        float pdfNormalization = 1;
        glm::vec3 w_in = importanceSample(normal, w_out, material, pdfNormalization);

        glm::vec3 f = brdf(normal, w_in, w_out, material);
        glm::vec3 T = f * glm::abs(glm::dot(w_in, normal)) / pdfNormalization;

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
                outputColor += boost * T * traceRay(position + w_in * EPSILON, w_in, numBounces);
            }
        }
        else
        {
            outputColor += T * traceRay(position + w_in * EPSILON, w_in, numBounces);
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
        return _volumetricBSDF->brdf(surfaceNormal, w_in, w_out, mat);
    else
        return _phongBRDF->brdf(surfaceNormal, w_in, w_out, mat);
}

/**
 * Geometry term for rendering equation
 */
float VolumetricPathTracerIntegrator::geometry(
    glm::vec3 surfacePoint,
    glm::vec3 surfaceNormal,
    glm::vec3 lightPoint,
    glm::vec3 lightNormal)
{
    glm::vec3 lightVector = glm::normalize(lightPoint - surfacePoint);

    float cosThetai = glm::abs(glm::dot(lightVector, surfaceNormal));
    float cosThetal = glm::abs(glm::dot(lightVector, -lightNormal));

    cosThetai = (cosThetai < 0) ? 0 : cosThetai;
    cosThetal = (cosThetal < 0) ? 0 : cosThetal;

    //return cosThetai * cosThetal;
    return cosThetai * cosThetal / glm::pow(glm::length(lightPoint - surfacePoint), 2);
}

glm::vec3 VolumetricPathTracerIntegrator::importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, float &pdfNormalization)
{
    glm::vec3 w_in;

    float epsilon1 = gen(rng);
    float epsilon2 = gen(rng);

    float theta = 0;
    float phi = 0;

    //float k_s = averageVector(material.specular);
    //float k_d = averageVector(material.diffuse);

    glm::vec3 samplingSpaceCenter = normal;
    //glm::vec3 reflection = (2 * glm::dot(normal, w_out) * normal - w_out);

    if (_scene->importanceSampling == COSINE_SAMPLING)
    {

        theta = glm::acos(glm::sqrt(epsilon1));
        phi = TWO_PI * epsilon2;

        w_in = sphereCoordsToVector(theta, phi, samplingSpaceCenter);
        pdfNormalization = glm::abs(glm::dot(normal, w_in)) / PI;
    }
    else if (_scene->importanceSampling == BRDF_SAMPLING)
    {
        if (material.brdf == GGX)
            w_in = _ggxBRDF->importanceSample(normal, w_out, material, pdfNormalization);
        else if (material.brdf == GGX_VOLUMETRIC)
            w_in = _volumetricBSDF->importanceSample(normal, w_out, material, pdfNormalization);
        else
            w_in = _phongBRDF->importanceSample(normal, w_out, material, pdfNormalization);
    }
    else
    {
        // hemisphere sampling
        theta = glm::acos(epsilon1);
        phi = TWO_PI * epsilon2;

        w_in = sphereCoordsToVector(theta, phi, samplingSpaceCenter);

        pdfNormalization = 1.0f / TWO_PI;
    }

    return w_in;
}

void VolumetricPathTracerIntegrator::setScene(Scene *scene)
{
    _scene = scene;

    _phongBRDF->setScene(scene);
    _ggxBRDF->setScene(scene);
    _volumetricBSDF->setScene(scene);
}