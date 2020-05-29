#include "MathUtils.h"

#include <glm/glm.hpp>

#include <iostream>

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

    /*
    // TODO: Check direction on w
    float cosTheta1 = glm::dot(halfVector, w);
    if (cosTheta1 < 0) cosTheta1 = -cosTheta1;

    float sinTheta2 = (ior_in / ior_out) * glm::sqrt(1 - glm::pow(cosTheta1, 2));

    float radicand = 1 - glm::pow(ior_in / ior_out, 2) * (1 - glm::pow(cosTheta1, 2));

    if (radicand < 0) {
        // total internal reflection, make sure this isn't bugged
        return glm::reflect(w, halfVector);
    }

    float cosTheta2 = glm::sqrt(radicand);
    return glm::vec3(0);
    */
}