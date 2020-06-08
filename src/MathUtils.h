#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <glm/glm.hpp>
#include <string>
#include <set>

#include "Scene.h"

glm::vec3 sphereCoordsToVector(float theta, float phi, glm::vec3 samplingSpaceCenter);

float averageVector(glm::vec3 vec);

float microfacetDistribution(float halfAngle, material_t material);
float microfacetSelfShadowing(glm::vec3 normal, glm::vec3 view, material_t material);
glm::vec3 fresnel(glm::vec3 w_in, glm::vec3 halfVector, material_t material);

glm::vec3 calculateRefraction(glm::vec3 halfVector, glm::vec3 w, float ior_in, float ior_out);
float fresnelIOR(glm::vec3 w_out, glm::vec3 normal, float ior_in, float ior_out);

float absdot(glm::vec3 a, glm::vec3 b);

volume_t findVolumeByID(Scene *scene, std::string id);
volume_t highestPriorityVolume(Scene *scene, std::set<std::string> volumes);

glm::vec3 attenuate(glm::vec3 T, float dist, volume_t volume);

#endif