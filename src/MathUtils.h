#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <glm/glm.hpp>

glm::vec3 sphereCoordsToVector(float theta, float phi, glm::vec3 samplingSpaceCenter);

float averageVector(glm::vec3 vec);

glm::vec3 calculateRefraction(glm::vec3 halfVector, glm::vec3 w, float ior_in, float ior_out);

#endif