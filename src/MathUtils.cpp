#include "MathUtils.h"

#include <glm/glm.hpp>

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