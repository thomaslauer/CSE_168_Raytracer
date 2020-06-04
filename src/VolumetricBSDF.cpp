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
        //normal = -normal;
    }

    // begin reflection calculation

    glm::vec3 h_r = glm::normalize(glm::sign(glm::dot(w_out, normal)) * (w_in + w_out));
    float halfAngle_r = glm::acos(glm::min(1.0f, glm::dot(h_r, normal)));

    float f_r_numerator = fresnelIOR(w_out, h_r, ior_in, ior_out) * microfacetDistribution(halfAngle_r, material) * microfacetSelfShadowing(normal, w_in, material) * microfacetSelfShadowing(normal, w_out, material);
    float f_r_denominator = 4 * absdot(w_in, normal) * absdot(w_out, normal);

    float f_r = f_r_numerator / f_r_denominator;

    // begin transmission calculation

    glm::vec3 h_t = glm::normalize((ior_in * w_in + ior_out * w_out));
    //std::cout << glm::dot(h_t, normal) << std::endl;
    float halfAngle_t = glm::acos(glm::min(1.0f, glm::dot(h_t, -normal)));

    float dotprods = absdot(w_out, h_t) * absdot(w_in, h_t) / (absdot(w_out, normal) * absdot(w_in, normal));

    float F_t = fresnelIOR(w_out, h_t, ior_in, ior_out);
    float G_t = microfacetSelfShadowing(-h_t, w_out, material) * microfacetSelfShadowing(h_t, w_in, material);
    float D_t = microfacetDistribution(halfAngle_t, material);

    float f_t = glm::pow(ior_in, 2) * (1 - F_t) * G_t * D_t / glm::pow(ior_out * glm::dot(w_out, h_t) + ior_in * glm::dot(w_in, h_t), 2);

    //return glm::vec3(1) * absdot(normal, w_out);
    return glm::vec3(1);
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

    if (glm::dot(w_out, halfVector) > 0)
    {
        // front face hit
        w_in_refraction = calculateRefraction(halfVector, -w_out, 1.0f, material.ior);
        w_in_reflection = glm::reflect(halfVector, -w_out);
        f = fresnelIOR(w_out, halfVector, material.ior, 1.0f);
    }
    else
    {
        // back face hit
        w_in_refraction = calculateRefraction(-halfVector, -w_out, material.ior, 1.0f);
        w_in_reflection = glm::reflect(-halfVector, -w_out);
        f = fresnelIOR(w_out, halfVector, 1.0f, material.ior);
    }

    glm::vec3 retval;
    if (f < epsilon3)
    {
        retval = w_in_reflection;
    }
    else
    {
        retval = w_in_refraction;
    }
    //retval = w_in_refraction;
    retval = w_in_reflection;

    // calculate weighting from paper

    float weight = absdot(w_out, halfVector) / (absdot(w_out, normal) * absdot(halfVector, normal));

    //pdfNormalization = absdot(w_in_refraction, normal);
    pdfNormalization = absdot(w_in_reflection, normal);

    return retval;
}

float VolumetricBSDF::pdf(
    glm::vec3 normal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t material)
{
    return 1;
}