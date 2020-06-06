#include "BRDF.h"

#include <iostream>
#include <glm/gtx/string_cast.hpp>

#include "MathUtils.h"
#include "Constants.h"

glm::vec3 reflectionBRDF(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material)
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
        normal = -normal;
        return glm::vec3(0);
    }

    glm::vec3 halfVector = glm::normalize((w_in + w_out));
    float halfAngle = glm::acos(glm::min(1.0f, glm::dot(halfVector, normal)));

    float D = microfacetDistribution(halfAngle, material);
    float G = microfacetSelfShadowing(halfVector, w_in, material) * microfacetSelfShadowing(halfVector, w_out, material);
    float F = fresnelIOR(w_out, halfVector, ior_in, ior_out);

    float normalization = 4 * absdot(w_in, normal) * absdot(w_out, normal);

    glm::vec3 retval = glm::vec3(1) * F * D * G / normalization;

    return glm::vec3(0);

    return retval;
}

glm::vec3 transmissionBTDF(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material)
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
        normal = -normal;
    }

    glm::vec3 halfVector = glm::normalize(ior_in * w_in + ior_out * w_out);
    float halfAngle = glm::acos(glm::min(1.0f, glm::dot(halfVector, normal)));

    float D = microfacetDistribution(halfAngle, material);
    float G = glm::abs(microfacetSelfShadowing(normal, w_in, material)) * glm::abs(microfacetSelfShadowing(normal, w_out, material));
    float F = fresnelIOR(w_out, normal, ior_in, ior_out);

    float dots = absdot(w_out, halfVector) * absdot(w_in, halfVector) / (absdot(w_out, normal) * absdot(w_in, halfVector));

    float rest = glm::pow(ior_in, 2.0f) * (1 - F) * G * D / glm::pow(ior_out * glm::dot(w_out, halfVector) + ior_in * glm::dot(w_in, halfVector), 2.0f);

    //std::cout << D << "  " << G << "  " << F << std::endl;
    //std::cout << dots << " " << rest << std::endl;

    return glm::vec3(1) * dots * rest;
}

glm::vec3 VolumetricBSDF::brdf(
    glm::vec3 normal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t material)
{
    return reflectionBRDF(normal, w_in, w_out, material) + transmissionBTDF(normal, w_in, w_out, material);
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
        normal = -normal;
    }

    glm::vec3 halfVector = sphereCoordsToVector(theta, phi, normal);

    f = fresnelIOR(w_out, halfVector, ior_in, ior_out);

    /*
    if (glm::dot(w_out, halfVector) < 0)
    {
        // back face hit
        w_in_refraction = calculateRefraction(-halfVector, w_out, material.ior, 1.0f);
        pdfNormalization = absdot(w_in_refraction, normal);
        return w_in_refraction;
    }
    */

    // front face hit
    w_in_refraction = calculateRefraction(halfVector, w_out, ior_in, ior_out);
    w_in_reflection = glm::reflect(-w_out, halfVector);

    glm::vec3 retval;

    f = 0;

    if (f > epsilon3)
    {
        retval = w_in_reflection;
    }
    else
    {
        retval = w_in_refraction;
    }

    //pdfNormalization = absdot(w_in_reflection, normal) * f + absdot(w_in_refraction, normal) * (1 - f);

    float halfAngle = glm::acos(glm::min(1.0f, glm::dot(halfVector, normal)));
    float ggxTerm = f * microfacetDistribution(halfAngle, material) * glm::dot(normal, halfVector) / (4.0f * glm::dot(halfVector, w_in_reflection));
    float transmissionTerm = 0;

    //pdfNormalization = ggxTerm + transmissionTerm;

    pdfNormalization = 1.0f / (absdot(w_out, halfVector) / (absdot(w_out, normal) * absdot(normal, halfVector)));

    return retval;
}

float VolumetricBSDF::pdf(
    glm::vec3 normal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t material)
{
    bool frontface;
    float ior_out;
    float ior_in;

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
        normal = -normal;
    }

    std::cout << "asdf" << std::endl;
    glm::vec3 halfVector = glm::normalize(w_in + w_out);
    float halfAngle = glm::acos(glm::min(1.0f, glm::dot(halfVector, normal)));

    float f = fresnelIOR(w_out, halfVector, ior_in, ior_out);
    float ggxTerm = f * microfacetDistribution(halfAngle, material) * glm::dot(normal, halfVector) / (4.0f * glm::dot(halfVector, w_in));

    float transmissionTerm = 0.2;

    return ggxTerm + transmissionTerm;
}