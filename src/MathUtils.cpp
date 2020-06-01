#include "MathUtils.h"

#include <glm/glm.hpp>
#include <iostream>

#include "Constants.h"
#include "Scene.h"

glm::vec3 sphereCoordsToVector(float theta, float phi, glm::vec3 samplingSpaceCenter) {
    // a sample over the unit hemisphere
    glm::vec3 s = glm::vec3(glm::cos(phi) * glm::sin(theta), glm::sin(phi) * glm::sin(theta), glm::cos(theta));

    // calculate the new coordinate frame
    glm::vec3 w = samplingSpaceCenter;

    // create arbitrariy a vector. it's a secret tool we will need later
    glm::vec3 a = glm::vec3(0, 1, 0);

    // if a isn't a good choice, pick a new one
    if (glm::length(w - a) < 0.001 || glm::length(w + a) < 0.001)
        a = glm::vec3(0, 0, 1);

    glm::vec3 u = glm::normalize(glm::cross(a, w));
    glm::vec3 v = glm::cross(w, u);

    return s.x * u + s.y * v + s.z * w;
}

float averageVector(glm::vec3 vec) {
    float avg = 0;
    avg += vec.x;
    avg += vec.y;
    avg += vec.z;
    return avg / 3.0f;
}

glm::vec3 calculateRefraction(glm::vec3 halfVector, glm::vec3 w, float ior_in, float ior_out) {

    glm::vec3 refraction = glm::refract(w, halfVector, ior_in / ior_out);
    if (refraction == glm::vec3(0)) {
        return glm::reflect(w, halfVector);
    }
    return refraction;

}

float microfacetDistribution(float halfAngle, material_t material) {
    float a_squared = material.roughness * material.roughness;

    float denominator = PI * glm::pow(glm::cos(halfAngle), 4.0f) * glm::pow(a_squared + glm::pow(glm::tan(halfAngle), 2.0f), 2.0f);
    return a_squared / denominator;
}

float microfacetSelfShadowing(glm::vec3 normal, glm::vec3 view, material_t material) {
    if (glm::dot(view, normal) <= 0)
        return 0;

    float thetaV = glm::acos(glm::dot(view, normal));

    return 2.0f / (1.0f + glm::sqrt(glm::max(0.0f, 1.0f + glm::pow(material.roughness, 2.0f) * glm::pow(glm::tan(thetaV), 2.0f))));
}

glm::vec3 fresnel(glm::vec3 w_in, glm::vec3 halfVector, material_t material) {
    return material.specular + (glm::vec3(1.0f) - material.specular) * glm::pow(1.0f - glm::dot(w_in, halfVector), 5.0f);
}