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
    glm::vec3 halfVector = glm::normalize(w_in + w_out);

    return fresnel(w_in, halfVector, material);
}

glm::vec3 VolumetricBSDF::importanceSample(
    glm::vec3 normal,
    glm::vec3 w_out,
    material_t material,
    float &pdfNormalization)
{
    glm::vec3 newDirection;

    float epsilon1 = gen(rng);
    float epsilon2 = gen(rng);
    float theta = glm::atan(material.roughness * glm::sqrt(epsilon1), glm::sqrt(1 - epsilon1));
    float phi = TWO_PI * epsilon2;

    glm::vec3 halfVector = sphereCoordsToVector(theta, phi, normal);

    if (glm::dot(w_out, halfVector) > 0)
    {
        // front face hit
        newDirection = calculateRefraction(halfVector, -w_out, 1.0f, material.ior);
    }
    else
    {
        // back face hit
        newDirection = calculateRefraction(-halfVector, -w_out, material.ior, 1.0f);
    }

    pdfNormalization = glm::min(glm::abs(glm::dot(normal, w_out)), 1.0f);
    return newDirection;
}

float VolumetricBSDF::pdf(
    glm::vec3 normal,
    glm::vec3 w_in,
    glm::vec3 w_out,
    material_t material)
{
    return 1;
}