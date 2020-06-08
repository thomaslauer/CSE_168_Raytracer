#pragma once

#include <string>
#include <vector>
#include <limits>

#include <glm/glm.hpp>
#include <embree3/rtcore.h>

class Integrator;

struct camera_t
{
    glm::vec3 origin;
    glm::vec3 imagePlaneTopLeft;
    glm::vec3 pixelRight;
    glm::vec3 pixelDown;
};

typedef enum brdf_t
{
    PHONG,
    GGX,
    GGX_VOLUMETRIC
} brdf_t;

struct triangleData_t
{
    bool interpolate;
    glm::vec3 normals[3];
    glm::vec3 barycentricCoords;
};

struct volume_t
{
    std::string id;
    int priority;
    float ior;
    glm::vec3 absorbsion;
    float meanScatterDistance;
    float scatterDirectionality;
};

struct material_t
{
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
    glm::vec3 emission;
    glm::vec3 ambient;
    float roughness;
    std::string volumeID;
    brdf_t brdf;
    bool light;
    triangleData_t triangleData;
};

struct directionalLight_t
{
    glm::vec3 toLight;
    glm::vec3 brightness;
};

struct pointLight_t
{
    glm::vec3 point;
    glm::vec3 brightness;
    glm::vec3 attenuation;
};

struct quadLight_t
{
    glm::vec3 a;
    glm::vec3 ab;
    glm::vec3 ac;
    glm::vec3 intensity;
};

typedef enum samplingMethod_t
{
    HEMISPHERE_SAMPLING,
    COSINE_SAMPLING,
    BRDF_SAMPLING
} samplingMethod_t;

class Scene
{

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

    std::vector<volume_t> volumeList;
    std::string defaultVolume;

    Integrator *integrator;
    samplingMethod_t importanceSampling;
    float gamma;

    RTCScene embreeScene;

    // lighting information for Monte Carlo direct integrator
    int lightSamples;
    bool lightStratify;
    int stratifyGridSize;

    int samplesPerPixel;

    bool MIS;
    bool nextEventEstimation;
    bool russianRoulette;

    bool castRay(
        glm::vec3 origin,
        glm::vec3 direction,
        glm::vec3 *hitPosition,
        glm::vec3 *hitNormal,
        material_t *hitMaterial) const;

    bool castOcclusionRay(
        glm::vec3 origin,
        glm::vec3 direction,
        float maxDistance = std::numeric_limits<float>::infinity()) const;
};
