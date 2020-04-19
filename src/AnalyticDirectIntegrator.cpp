#include "Integrator.h"

#include <iostream>

#include "Constants.h"

float AnalyticDirectIntegrator::calculateSubtendedAngle(glm::vec3 r, glm::vec3 p1, glm::vec3 p2) {
    glm::vec3 temp1 = glm::normalize(p1 - r);
    glm::vec3 temp2 = glm::normalize(p2 - r);
    float result = glm::acos(glm::dot(temp1, temp2));
    return result;
}

glm::vec3 AnalyticDirectIntegrator::calculateSubtendedAngleNormal(glm::vec3 r, glm::vec3 p1, glm::vec3 p2) {
    glm::vec3 temp1 = p1 - r;
    glm::vec3 temp2 = p2 - r;
    glm::vec3 result = glm::normalize(glm::cross(temp1, temp2));
    return result;
}

glm::vec3 AnalyticDirectIntegrator::traceRay(glm::vec3 origin, glm::vec3 direction) {
    glm::vec3 outputColor = glm::vec3(0.0f, 0.0f, 0.0f);

    glm::vec3 hitPosition;
    glm::vec3 hitNormal;
    material_t hitMaterial;
    bool hit = _scene->castRay(origin, direction, &hitPosition, &hitNormal, &hitMaterial);

    if (hit) {
        for (auto quadLight : _scene->quadLights) {
            glm::vec3 a = quadLight.a;
            glm::vec3 b = quadLight.a + quadLight.ab;
            glm::vec3 c = quadLight.a + quadLight.ac;
            glm::vec3 d = quadLight.a + quadLight.ab + quadLight.ac;

            glm::vec3 lights[4] = {a, b, d, c};

            glm::vec3 irradiance = glm::vec3(0);
            for (int i = 0; i < 4; i++) {
                irradiance += calculateSubtendedAngle(hitPosition, lights[i], lights[(i+1)%4])
                            * calculateSubtendedAngleNormal(hitPosition, lights[i], lights[(i+1)%4]);
            }
            irradiance *= 0.5f;

            glm::vec3 f = hitMaterial.diffuse * INV_PI;

            outputColor += f * quadLight.intensity * glm::dot(hitNormal, irradiance);
        }
    }

    return outputColor;
}