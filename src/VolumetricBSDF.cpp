#include "BRDF.h"

#include <iostream>
#include <glm/gtx/string_cast.hpp>

#include "MathUtils.h"
#include "Constants.h"

glm::vec3 VolumetricBSDF::brdf(
    glm::vec3 normal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t material)
{
    float ior_in;
    float ior_out;

    bool frontface;

    if (glm::dot(w_out, normal) > 0)
    {
        // front face hit
        frontface = true;
        ior_out = 1.0f;
        ior_in = material.ior;
    }
    else
    {
        // back face hit
        frontface = false;
        ior_out = material.ior;
        ior_in = 1.0f;
    }

    return glm::vec3(1);
    // begin reflection calculation

    if (frontface)
    {
        if (glm::dot(w_in, normal) > 0 && glm::dot(w_out, normal) > 0)
        {
            // reflection ray
            glm::vec3 h_r = glm::normalize(glm::sign(glm::dot(w_out, normal)) * (w_in + w_out));
            float halfAngle_r = glm::acos(glm::min(1.0f, glm::dot(h_r, normal)));

            float D = microfacetDistribution(halfAngle_r, material);
            float G = microfacetSelfShadowing(normal, w_in, material) * microfacetSelfShadowing(normal, w_out, material);
            float F = fresnelIOR(w_out, normal, material.ior, 1.0f);

            float normalization = 4 * glm::dot(w_in, normal) * glm::dot(w_out, normal);

            glm::vec3 retval = glm::vec3(1) * D * G * glm::pow(F, 2.0f) / normalization;
            return 1.0f / retval;
        }
        else
        {
            // transmission ray
            // TODO: Add correct brsf from paper
            float F = fresnelIOR(w_out, normal, material.ior, 1.0f);
            return glm::vec3(1);
        }
    }
    else
    {
        // backface hit
        // TODO: also add brsf here, switching for backfaces
        return glm::vec3(1);
    }
}

glm::vec3 VolumetricBSDF::importanceSample(
    glm::vec3 normal,
    glm::vec3 w_out,
    material_t material,
    float &pdfNormalization)
{
    glm::vec3 w_in_refraction;
    glm::vec3 w_in_reflection;
    float f;

    float epsilon1 = gen(rng);
    float epsilon2 = gen(rng);
    float epsilon3 = gen(rng);
    float theta = glm::atan(material.roughness * glm::sqrt(epsilon1), glm::sqrt(1 - epsilon1));
    float phi = TWO_PI * epsilon2;

    glm::vec3 halfVector = sphereCoordsToVector(theta, phi, normal);

    f = fresnelIOR(w_out, halfVector, material.ior, 1.0f);

    if (glm::dot(w_out, halfVector) > 0)
    {
        // front face hit
        w_in_refraction = calculateRefraction(halfVector, -w_out, 1.0f, material.ior);
        w_in_reflection = glm::reflect(-w_out, halfVector);
    }
    else
    {
        // back face hit
        w_in_refraction = calculateRefraction(-halfVector, -w_out, material.ior, 1.0f);
        pdfNormalization = absdot(w_in_refraction, normal);
        return w_in_refraction;
    }

    glm::vec3 retval;

    if (f > epsilon3)
    {
        retval = w_in_reflection;
    }
    else
    {
        retval = w_in_refraction;
    }

    pdfNormalization = absdot(w_in_reflection, normal) * f + absdot(w_in_refraction, normal) * (1 - f);
    return retval;
}

float VolumetricBSDF::pdf(
    glm::vec3 normal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t material)
{
    glm::vec3 halfVector = glm::normalize(w_in + w_out);
    float halfAngle = glm::acos(glm::min(1.0f, glm::dot(halfVector, normal)));
    float F = fresnelIOR(w_out, normal, material.ior, 1.0f);

    float ggxTerm = glm::pow(F, 2.0f) * microfacetDistribution(halfAngle, material) * glm::dot(normal, halfVector) / (4.0f * glm::dot(halfVector, w_in));

    return ggxTerm;
}