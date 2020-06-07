#include "BRDF.h"

#include <iostream>
#include <glm/gtx/string_cast.hpp>

#include "MathUtils.h"
#include "Constants.h"

thread_local std::default_random_engine VolumetricBSDF::rng = std::default_random_engine((std::random_device())());
VolumetricBSDF::VolumetricBSDF()
{
    gen = std::uniform_real_distribution<float>(0.0f, 1.0f);
}

glm::vec3 VolumetricBSDF::importanceSample(
    glm::vec3 normal,
    glm::vec3 w_out,
    material_t material,
    std::set<std::string> &volumes,
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

    volume_t hitVolume = _scene->volumeMap[material.volumeID];

    glm::vec3 halfVector = sphereCoordsToVector(theta, phi, normal);

    f = fresnelIOR(w_out, halfVector, hitVolume.ior, 1.0f);

    if (glm::dot(w_out, halfVector) < 0)
    {
        // back face hit
        volumes.erase(hitVolume.id);
        w_in_refraction = calculateRefraction(-halfVector, w_out, hitVolume.ior, 1.0f);
        pdfNormalization = absdot(w_in_refraction, normal);
        return w_in_refraction;
    }

    // front face hit
    w_in_refraction = calculateRefraction(halfVector, w_out, 1.0f, hitVolume.ior);
    w_in_reflection = glm::reflect(-w_out, halfVector);

    glm::vec3 retval;

    if (f > epsilon3)
    {
        retval = w_in_reflection;
    }
    else
    {
        retval = w_in_refraction;
        volumes.insert(material.volumeID);
    }

    pdfNormalization = absdot(w_in_reflection, normal) * f + absdot(w_in_refraction, normal) * (1 - f);
    return retval;
}