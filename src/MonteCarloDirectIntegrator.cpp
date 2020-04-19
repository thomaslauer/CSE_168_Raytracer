#include "Integrator.h"

#include <iostream>

#include "Constants.h"

MonteCarloDirectIntegrator::MonteCarloDirectIntegrator() {
    rng = boost::random::mt19937();
    gen = boost::random::uniform_real_distribution<float>(0.0f, 1.0f);
}

glm::vec3 MonteCarloDirectIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction) {
    glm::vec3 outputColor = glm::vec3(0);

    glm::vec3 hitPosition;
    glm::vec3 hitNormal;
    material_t hitMaterial;
    bool hit = _scene->castRay(origin, direction, &hitPosition, &hitNormal, &hitMaterial);

    if (hit) {
        if (hitMaterial.light) return hitMaterial.emission;

        for (auto light : _scene->quadLights) {
            // TODO: select random location on light
            float offsetB = gen(rng);
            float offsetC = gen(rng);
            glm::vec3 lightPosition = light.a + offsetB * light.ab + offsetC * light.ac;
            float lightArea = glm::length(glm::cross(light.ab, light.ac));

            float V = occlusion(hitPosition, lightPosition);

            glm::vec3 lightNormal = glm::normalize(glm::cross(light.ac, light.ab));

            float G = geometry(origin, hitPosition, hitNormal, lightPosition, lightNormal);

            G = (G < 0) ? 0 : G;

            outputColor += brdf(hitMaterial) * lightArea * light.intensity * V * G;
            //outputColor += brdf(hitMaterial) * V;

        }
    }

    return outputColor;
}

float MonteCarloDirectIntegrator::geometry(
    glm::vec3 origin,
    glm::vec3 surfacePoint,
    glm::vec3 surfaceNormal,
    glm::vec3 lightPoint,
    glm::vec3 lightNormal)
{
    //float cosThetai = glm::dot(surfacePoint - origin, surfaceNormal) / glm::length(surfacePoint-origin);
    //float cosThetal = -glm::dot(lightPoint - surfacePoint, lightNormal) / glm::length(lightPoint-surfacePoint);

    glm::vec3 lightVector = glm::normalize(lightPoint - surfacePoint);

    float cosThetai = glm::dot(lightVector, surfaceNormal);
    float cosThetal = glm::dot(lightVector, -lightNormal);

    cosThetai = (cosThetai < 0) ? 0 : cosThetai;
    cosThetal = (cosThetal < 0) ? 0 : cosThetal;

    //return cosThetai * cosThetal;
    return cosThetai * cosThetal / glm::pow(glm::length(lightPoint - surfacePoint), 2);
}

glm::vec3 MonteCarloDirectIntegrator::brdf(material_t mat) {
    glm::vec3 diffuse = mat.diffuse * INV_PI;
    return diffuse;
}

float MonteCarloDirectIntegrator::occlusion(glm::vec3 origin, glm::vec3 target) {
    glm::vec3 direction = target - origin;
    bool occluded = _scene->castOcclusionRay(origin, glm::normalize(target - origin), glm::length(target-origin));
    if (occluded) return 0;
    return 1;
}