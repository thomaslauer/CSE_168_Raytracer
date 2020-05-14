#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <glm/glm.hpp>

glm::vec3 sphereCoordsToVector(float theta, float phi, glm::vec3 samplingSpaceCenter);

float averageVector(glm::vec3 vec);

#endif