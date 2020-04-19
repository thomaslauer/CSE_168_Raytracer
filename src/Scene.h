#pragma once

#include <string>
#include <vector>
#include <limits>

#include <glm/glm.hpp>
#include <embree3/rtcore.h>

class Integrator;

struct camera_t {
    glm::vec3 origin;
    glm::vec3 imagePlaneTopLeft;
    glm::vec3 pixelRight;
    glm::vec3 pixelDown;
};

struct material_t {
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
    glm::vec3 emission;
    glm::vec3 ambient;
    bool light;
};

struct directionalLight_t {
    glm::vec3 toLight;
    glm::vec3 brightness;
};

struct pointLight_t {
    glm::vec3 point;
    glm::vec3 brightness;
    glm::vec3 attenuation;
};

struct quadLight_t {
    glm::vec3 a;
    glm::vec3 ab;
    glm::vec3 ac;
    glm::vec3 intensity;
};

class Scene {

public:

    glm::uvec2 imageSize;
    int maxDepth;
    std::string outputFileName;
    camera_t camera;
    std::vector<glm::mat3> sphereNormalTransforms;
    std::vector<material_t> sphereMaterials;
    std::vector<material_t> triMaterials;
    std::vector<directionalLight_t> directionalLights;
    std::vector<pointLight_t> pointLights;
    std::vector<quadLight_t> quadLights;

    Integrator* integrator;

    RTCScene embreeScene;

    int lightSamples;
    bool lightStratify;
    int stratifyGridSize;

    bool castRay(
        glm::vec3 origin,
        glm::vec3 direction,
        glm::vec3* hitPosition,
        glm::vec3* hitNormal,
        material_t* hitMaterial) const;

    bool castOcclusionRay(
        glm::vec3 origin,
        glm::vec3 direction,
        float maxDistance = std::numeric_limits<float>::infinity()) const;

};
