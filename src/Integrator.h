#pragma once

#include <glm/glm.hpp>

#include "Scene.h"

class Integrator {

protected:

    Scene* _scene;

public:

    void setScene(Scene* scene)
    {
        _scene = scene;
    }

    virtual glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction) = 0;

};


class RayTracerIntegrator : public Integrator {

private:

    glm::vec3 computeShading(
        glm::vec3 incidentDirection,
        glm::vec3 toLight,
        glm::vec3 normal,
        glm::vec3 lightBrightness,
        const material_t& material);

    glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction, int depth);

public:

    virtual glm::vec3 traceRay(glm::vec3 origin, glm::vec3 direction);

};
