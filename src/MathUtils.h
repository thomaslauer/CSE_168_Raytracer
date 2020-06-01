#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <glm/glm.hpp>

#include "Scene.h"

glm::vec3 sphereCoordsToVector(float theta, float phi, glm::vec3 samplingSpaceCenter);

float averageVector(glm::vec3 vec);

glm::vec3 calculateRefraction(glm::vec3 halfVector, glm::vec3 w, float ior_in, float ior_out);

float microfacetDistribution(float halfAngle, material_t material);
float microfacetSelfShadowing(glm::vec3 normal, glm::vec3 view, material_t material);
glm::vec3 fresnel(glm::vec3 w_in, glm::vec3 halfVector, material_t material);

#endif