#include "Integrator.h"

#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>

#include "BRDF.h"
#include "Constants.h"
#include "MathUtils.h"

thread_local std::default_random_engine PathTracerIntegrator::rng = std::default_random_engine((std::random_device())());

PathTracerIntegrator::PathTracerIntegrator()
{
    gen = std::uniform_real_distribution<float>(0.0f, 1.0f);

    _phongBRDF = new PhongBRDF();
    _ggxBRDF = new GGXBRDF();
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
        if (_scene->nextEventEstimation || _scene->MIS)
        {
            // end conditions for NEE
            if (hitMaterial.light)
            {
                if (numBounces == 0 && glm::dot(hitNormal, direction) <= 0) {
                    outputColor += hitMaterial.emission;
                }
                else {
                    return glm::vec3(0);
                }
            }
            if (!_scene->russianRoulette && numBounces >= _scene->maxDepth)
                return glm::vec3(0);
        }
        else
        {
            if (hitMaterial.light || (!_scene->russianRoulette && numBounces > _scene->maxDepth))
                return hitMaterial.emission;
        }

        if (_scene->MIS) {
            float neeWeighting;
            glm::vec3 neeColor = nextEventEstimation(
                hitPosition,
                hitNormal,
                hitMaterial,
                origin,
                neeWeighting);

            float brdfWeighting;
            glm::vec3 brdfColor = ggxDirect(
                hitPosition,
                hitNormal,
                hitMaterial,
                origin,
                brdfWeighting);
            
            //brdfWeighting = 1 - brdfWeighting;

            //std::cout << brdfWeighting << " " << neeWeighting << std::endl;

            //outputColor += neeWeighting * neeColor;
            outputColor += neeColor;
            outputColor += brdfWeighting * brdfColor;
            //outputColor += brdfWeighting * glm::vec3(1);

        } else if (_scene->nextEventEstimation) {
            float neePDF;
            outputColor += nextEventEstimation(
                hitPosition,
                hitNormal,
                hitMaterial,
                origin,
                neePDF);
        }

        outputColor += indirectLighting(
            hitPosition,
            hitNormal,
            hitMaterial,
            origin,
            numBounces + 1);

        /*
        outputColor += hitMaterial.emission;
        */
    }

    return outputColor;
}

glm::vec3 PathTracerIntegrator::nextEventEstimation(
    glm::vec3 position,
    glm::vec3 normal,
    material_t material,
    glm::vec3 origin,
    float& pdfNormalization)
{
    glm::vec3 outputColor = glm::vec3(0);
    int stratifyGridSize = _scene->stratifyGridSize;
    pdfNormalization = 0;
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

            glm::vec3 F = brdf(normal, w_in, w_out, material);
            float V = occlusion(position, lightPosition);
            float G = geometry(position, normal, lightPosition, lightNormal);

            float currentpdf = neePDF(position, w_in);

            pdfNormalization = brdfMisWeighting(position, normal, w_in, w_out, material, false);
            outputColor += lightArea * light.intensity * F * V * G * pdfNormalization / ((float)_scene->lightSamples);
        }
    }
    pdfNormalization /= (float)_scene->quadLights.size();
    return outputColor;
}

float PathTracerIntegrator::neePDF(glm::vec3 position, glm::vec3 w_in) {
    float p = 0;
    for (auto light : _scene->quadLights) {
        // dummy variable I don't actually use
        glm::vec2 baryPos;

        const float EPSILON1 = 0.01f;
        const float EPSILON2 = 1.01f;

        // sample the first triangle
        float distA = 0;
        bool hitA = glm::intersectRayTriangle(
            position,
            w_in, 
            light.a - EPSILON1 * (light.ab + light.ac),
            light.a + EPSILON2 * light.ab,
            light.a + EPSILON2 * light.ac,
            baryPos,
            distA);

        // sample the second triangle
        float distB = 0;
        bool hitB = glm::intersectRayTriangle(
            position,
            w_in, 
            light.a + EPSILON2 * (light.ab + light.ac),
            light.a + EPSILON2 * light.ab,
            light.a + EPSILON2 * light.ac,
            baryPos,
            distB);
        
        float actualDist;
        if (hitA && hitB) {
            actualDist = glm::max(distA, distB);
        } else if (!hitA && hitB) {
            actualDist = distB;
        } else if (hitA && !hitB) {
            actualDist = distA;
        } else {
            // missed both, shouldn't contribute at all
            continue;
        }

        float lightArea = glm::length(glm::cross(light.ab, light.ac));
        glm::vec3 lightNormal = glm::normalize(glm::cross(light.ac, light.ab));

        p += glm::pow(actualDist, 2) / (lightArea * glm::abs(glm::dot(w_in, lightNormal)));
    }
    p = p / (float)_scene->quadLights.size();

    return p;
}

float PathTracerIntegrator::brdfMisWeighting(glm::vec3 position, glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material, bool use_brdf) {
    float neePDFValue = neePDF(position, w_in);
    float brdfPDFValue = pdf(normal, w_in, w_out, material);

    //std::cout << "NEE PDF is " << neePDFValue << " brdf pdf is " << brdfPDFValue << std::endl;

    float neePDF2 = glm::pow(neePDFValue, 2);
    float brdfPDF2 = glm::pow(brdfPDFValue, 2);

    //std::cout << "pdf: " << neePDFValue << " " << brdfPDFValue << std::endl;

    float denom = neePDF2 + brdfPDF2;

    if (use_brdf) {
        return brdfPDF2 / denom;
    } else {
        return neePDF2 / denom;
    }
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
    glm::vec3 w_out = glm::normalize(origin - position);
    
    for (int i = 0; i < numRaysPerBounce; i++)
    {
        float pdfNormalization;
        glm::vec3 w_in = importanceSample(normal, w_out, material, pdfNormalization);

        glm::vec3 f = brdf(normal, w_in, w_out, material);
        glm::vec3 T = f * glm::dot(w_in, normal) / pdfNormalization;

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
                outputColor += boost * T * traceRay(position, w_in, numBounces);
            }
        }
        else
        {
            outputColor += T * traceRay(position, w_in, numBounces);
        }
    }
    return outputColor / ((float)numRaysPerBounce);
}

glm::vec3 PathTracerIntegrator::ggxDirect(
    glm::vec3   position,
    glm::vec3   normal,
    material_t  material,
    glm::vec3   origin,
    float&      pdfNormalization) 
{
    glm::vec3 outputColor = glm::vec3(0);
    glm::vec3 w_out = glm::normalize(origin - position);

    float dummy;
    glm::vec3 w_in = importanceSample(normal, w_out, material, dummy);
    
    // hit parameters
    glm::vec3 hitPosition;
    glm::vec3 hitNormal;

    material_t hitMaterial;

    bool hit = _scene->castRay(position, w_in, &hitPosition, &hitNormal, &hitMaterial);

    if (hit) {
        glm::vec3 f = brdf(normal, w_in, w_out, material);
        glm::vec3 T = f * glm::dot(w_in, normal);
        float G = geometry(position, normal, hitPosition, hitNormal);
        pdfNormalization = brdfMisWeighting(position, normal, w_in, w_out, material, true);
        outputColor = T * hitMaterial.emission / pdf(normal, w_in, w_out, material);
    }

    return outputColor;
}

inline glm::vec3 PathTracerIntegrator::brdf(
    glm::vec3 surfaceNormal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t mat)
{
    if (mat.brdf == GGX)
        return _ggxBRDF->brdf(surfaceNormal, w_in, w_out, mat);
    else
        return _phongBRDF->brdf(surfaceNormal, w_in, w_out, mat);
}

inline float PathTracerIntegrator::pdf(
    glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) 
{
    if (material.brdf == GGX)
        return _ggxBRDF->pdf(normal, w_in, w_out, material);
    else
        return _phongBRDF->pdf(normal, w_in, w_out, material);
}

/**
 * Geometry term for rendering equation
 */
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

/**
 * Occlusion term for rendering equation
 */
float PathTracerIntegrator::occlusion(glm::vec3 origin, glm::vec3 target)
{
    glm::vec3 direction = target - origin;
    bool occluded = _scene->castOcclusionRay(origin, glm::normalize(direction), glm::length(direction));
    if (occluded)
        return 0;
    return 1;
}


glm::vec3 PathTracerIntegrator::importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, float &pdfNormalization) {
    glm::vec3 w_in;

    float epsilon1 = gen(rng);
    float epsilon2 = gen(rng);

    float theta = 0;
    float phi = 0;

    //float k_s = averageVector(material.specular);
    //float k_d = averageVector(material.diffuse);

    glm::vec3 samplingSpaceCenter = normal;
    //glm::vec3 reflection = (2 * glm::dot(normal, w_out) * normal - w_out);

    if (_scene->importanceSampling == COSINE_SAMPLING) {

        theta = glm::acos(glm::sqrt(epsilon1));
        phi = TWO_PI * epsilon2;

        w_in = sphereCoordsToVector(theta, phi, samplingSpaceCenter);
        pdfNormalization = glm::dot(normal, w_in) / PI;

    } else if (_scene->importanceSampling == BRDF_SAMPLING) {
        if (material.brdf == GGX) {
            w_in = _ggxBRDF->importanceSample(normal, w_out, material, pdfNormalization);
        } else {
            w_in = _phongBRDF->importanceSample(normal, w_out, material, pdfNormalization);
        }
    } else {
        // hemisphere sampling
        theta = glm::acos(epsilon1);
        phi = TWO_PI * epsilon2;

        w_in = sphereCoordsToVector(theta, phi, samplingSpaceCenter);

        pdfNormalization = 1.0f / TWO_PI;
    }

    return w_in;
}

void PathTracerIntegrator::setScene(Scene* scene) {
    _scene = scene;

    _phongBRDF->setScene(scene);
    _ggxBRDF->setScene(scene);
}