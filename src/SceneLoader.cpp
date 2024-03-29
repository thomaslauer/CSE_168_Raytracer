#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <limits>
#include <stdexcept>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Constants.h"
#include "Scene.h"
#include "Integrator.h"

#include "SceneLoader.h"

#include "tiny_obj_loader.h"

class SceneLoader
{

private:
    RTCDevice _embreeDevice;

    glm::uvec2 _imageSize = glm::uvec2(1280, 720);
    int _maxDepth = 5;
    std::string _outputFileName = "out.png";
    glm::vec3 _cameraOrigin = glm::vec3(-1.0f, 0.0f, 0.0f);
    glm::vec3 _cameraLookAt = glm::vec3(0.0f, 0.0f, 0.0f);
    glm::vec3 _cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
    float _cameraFieldOfView = 45.0f;
    std::vector<glm::mat4> _sphereTransforms;
    std::vector<material_t> _sphereMaterials;
    std::vector<glm::vec3> _rawVertices;
    std::vector<glm::uvec3> _indices;
    std::vector<glm::vec3> _vertices;
    std::vector<material_t> _triMaterials;
    glm::mat4 curTransform = glm::mat4(1.0f);
    std::vector<glm::mat4> _transformStack;
    std::vector<directionalLight_t> _directionalLights;
    std::vector<pointLight_t> _pointLights;
    std::vector<quadLight_t> _quadLights;
    glm::vec3 _curAttenuation = glm::vec3(1.0f, 0.0f, 0.0f);

    std::string _defaultVolume;
    std::vector<volume_t> _volumeList;

    material_t _curMaterial = {
        glm::vec3(0.0f),             // diffuse
        glm::vec3(0.0f),             // specular
        1.0f,                        // shininess
        glm::vec3(0.0f),             // emission
        glm::vec3(0.2f, 0.2f, 0.2f), // ambient
        0.0f,                        // roughness
        "air",
        PHONG, // brdf type
        false, // light
        {
            // triangle interpolation data
            false,                                      // do interpolation?
            {glm::vec3(0), glm::vec3(0), glm::vec3(0)}, // vertex normals
            glm::vec3(0)                                // barycentric coords
        }};

    std::string _integratorType = "raytracer";
    int _lightSamples = 1;
    bool _lightStratify = false;

    int _samplesPerPixel = 1;

    bool _MIS = false;
    bool _nextEventEstimation = false;
    bool _russianRoulette = false;
    float _gamma = 1;
    samplingMethod_t _importanceSampling = HEMISPHERE_SAMPLING;

    void quadLightToTriangles();

public:
    SceneLoader(RTCDevice embreeDevice);
    glm::vec3 loadVec3(const std::vector<std::string> &arguments, size_t startIndex = 0);
    glm::uvec3 loadUVec3(const std::vector<std::string> &arguments, size_t startIndex = 0);
    void executeCommand(const std::string &command, const std::vector<std::string> &arguments);
    void loadSceneData(const std::string &filePath);
    Integrator *createIntegrator();
    void loadEmbreeTriangles(RTCScene embreeScene);
    void loadEmbreeSpheres(RTCScene embreeScene);
    RTCScene createEmbreeScene();
    Scene *commitSceneData();
};

SceneLoader::SceneLoader(RTCDevice embreeDevice)
    : _embreeDevice(embreeDevice)
{
}

glm::vec3 SceneLoader::loadVec3(const std::vector<std::string> &arguments, size_t startIndex)
{
    return glm::vec3(
        std::stof(arguments[startIndex]),
        std::stof(arguments[startIndex + 1]),
        std::stof(arguments[startIndex + 2]));
}

glm::uvec3 SceneLoader::loadUVec3(const std::vector<std::string> &arguments, size_t startIndex)
{
    return glm::uvec3(
        std::stoi(arguments[startIndex]),
        std::stoi(arguments[startIndex + 1]),
        std::stoi(arguments[startIndex + 2]));
}

void SceneLoader::executeCommand(
    const std::string &command,
    const std::vector<std::string> &arguments)
{
    if (command == "size")
    {

        _imageSize = glm::uvec2(std::stoi(arguments[0]), std::stoi(arguments[1]));
    }
    else if (command == "maxdepth")
    {

        _maxDepth = std::stoi(arguments[0]);
        if (_maxDepth == -1)
            _maxDepth = std::numeric_limits<int>::max();
    }
    else if (command == "output")
    {

        _outputFileName = arguments[0];
    }
    else if (command == "camera")
    {

        _cameraOrigin = loadVec3(arguments, 0);
        _cameraLookAt = loadVec3(arguments, 3);
        _cameraUp = loadVec3(arguments, 6);
        _cameraFieldOfView = std::stof(arguments[9]);
    }
    else if (command == "sphere")
    {

        glm::vec3 center = loadVec3(arguments, 0);
        float radius = std::stof(arguments[3]);

        glm::mat4 transform = glm::mat4(1.0f);
        transform = curTransform * transform;
        transform = glm::translate(transform, center);
        transform = glm::scale(transform, glm::vec3(radius));

        _curMaterial.triangleData.interpolate = false;

        _sphereTransforms.push_back(transform);

        _sphereMaterials.push_back(_curMaterial);
    }
    else if (command == "maxverts")
    {

        // ignore since we are using std::vector
    }
    else if (command == "vertex")
    {

        _rawVertices.push_back(loadVec3(arguments));
    }
    else if (command == "tri")
    {

        glm::uvec3 rawIndices = loadUVec3(arguments);

        _indices.push_back(glm::uvec3(
            _vertices.size(),
            _vertices.size() + 1,
            _vertices.size() + 2));

        _vertices.push_back(glm::vec3(curTransform * glm::vec4(_rawVertices[rawIndices.x], 1.0f)));
        _vertices.push_back(glm::vec3(curTransform * glm::vec4(_rawVertices[rawIndices.y], 1.0f)));
        _vertices.push_back(glm::vec3(curTransform * glm::vec4(_rawVertices[rawIndices.z], 1.0f)));

        _triMaterials.push_back(_curMaterial);
    }
    else if (command == "obj")
    {
        // load OBJ file using tinyobjloader

        // mostly taken fron the tinyobjloader example
        std::string filename = arguments[0];
        tinyobj::attrib_t attrib;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;

        std::string err;

        tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename.c_str());

        // Loop over shapes
        for (size_t s = 0; s < shapes.size(); s++)
        {
            // Loop over faces(polygon)
            size_t index_offset = 0;
            for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
            {
                size_t fv = shapes[s].mesh.num_face_vertices[f];

                _indices.push_back(glm::uvec3(
                    _vertices.size(),
                    _vertices.size() + 1,
                    _vertices.size() + 2));

                // Loop over vertices in the face.
                for (size_t v = 0; v < fv; v++)
                {
                    // access to vertex
                    tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                    tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
                    tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
                    tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
                    tinyobj::real_t nx = attrib.normals[3 * idx.normal_index + 0];
                    tinyobj::real_t ny = attrib.normals[3 * idx.normal_index + 1];
                    tinyobj::real_t nz = attrib.normals[3 * idx.normal_index + 2];
                    // tinyobj::real_t tx = attrib.texcoords[2 * idx.texcoord_index + 0];
                    // tinyobj::real_t ty = attrib.texcoords[2 * idx.texcoord_index + 1];
                    // Optional: vertex colors
                    // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
                    // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
                    // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];

                    glm::vec3 rawVertex = glm::vec3(vx, vy, vz);
                    //rawVertex = (rawVertex - center) / maximumDimension;
                    _curMaterial.triangleData.normals[v] = glm::vec3(nx, ny, nz);
                    _vertices.push_back(glm::vec3(curTransform * glm::vec4(rawVertex, 1.0f)));
                }
                index_offset += fv;

                _curMaterial.triangleData.interpolate = true;
                _triMaterials.push_back(_curMaterial);

                // per-face material
                shapes[s].mesh.material_ids[f];
            }
        }
        // stop interpolating new materials
        _curMaterial.triangleData.interpolate = false;
    }
    else if (command == "translate")
    {

        glm::vec3 translation = loadVec3(arguments);
        curTransform = glm::translate(curTransform, translation);
    }
    else if (command == "rotate")
    {

        glm::vec3 axis = loadVec3(arguments, 0);
        float radians = std::stof(arguments[3]) * PI / 180.0f;
        curTransform = glm::rotate(curTransform, radians, axis);
    }
    else if (command == "scale")
    {

        glm::vec3 scale = loadVec3(arguments);
        curTransform = glm::scale(curTransform, scale);
    }
    else if (command == "pushTransform")
    {

        _transformStack.push_back(curTransform);
    }
    else if (command == "popTransform")
    {

        curTransform = _transformStack.back();
        _transformStack.pop_back();
    }
    else if (command == "directional")
    {

        directionalLight_t light;
        light.toLight = glm::normalize(loadVec3(arguments, 0));
        light.brightness = loadVec3(arguments, 3);

        _directionalLights.push_back(light);
    }
    else if (command == "point")
    {

        pointLight_t light;
        light.point = loadVec3(arguments, 0);
        light.brightness = loadVec3(arguments, 3);
        light.attenuation = _curAttenuation;

        _pointLights.push_back(light);
    }
    else if (command == "quadLight")
    {

        quadLight_t light;
        light.a = loadVec3(arguments, 0);
        light.ab = loadVec3(arguments, 3);
        light.ac = loadVec3(arguments, 6);
        light.intensity = loadVec3(arguments, 9);

        _quadLights.push_back(light);
    }
    else if (command == "attenuation")
    {
        _curAttenuation = loadVec3(arguments);
    }
    else if (command == "ambient")
    {
        _curMaterial.ambient = loadVec3(arguments);
    }
    else if (command == "diffuse")
    {
        _curMaterial.diffuse = loadVec3(arguments);
    }
    else if (command == "specular")
    {
        _curMaterial.specular = loadVec3(arguments);
    }
    else if (command == "shininess")
    {
        _curMaterial.shininess = std::stof(arguments[0]);
    }
    else if (command == "emission")
    {
        _curMaterial.emission = loadVec3(arguments);
    }
    else if (command == "roughness")
    {
        _curMaterial.roughness = std::stof(arguments[0]);
    }
    else if (command == "volumetype")
    {
        std::string id = arguments[0];

        volume_t curVolume;
        curVolume.id = id;
        curVolume.priority = std::stoi(arguments[1]);
        curVolume.ior = std::stof(arguments[2]);
        curVolume.absorption = loadVec3(arguments, 3);
        curVolume.meanScatterDistance = std::stof(arguments[6]);
        curVolume.scatterDirectionality = std::stof(arguments[7]);

        _volumeList.push_back(curVolume);
    }
    else if (command == "volume")
    {
        _curMaterial.volumeID = arguments[0];
    }
    else if (command == "defaultvolume")
    {
        _defaultVolume = arguments[0];
    }
    else if (command == "integrator")
    {
        _integratorType = arguments[0];
    }
    else if (command == "lightsamples")
    {
        _lightSamples = std::stoi(arguments[0]);
    }
    else if (command == "lightstratify")
    {
        if (arguments[0] == "on")
            _lightStratify = true;
        else
            _lightStratify = false;
    }
    else if (command == "spp")
    {
        _samplesPerPixel = std::stoi(arguments[0]);
    }
    else if (command == "nexteventestimation")
    {
        if (arguments[0] == "on")
            _nextEventEstimation = true;
        else if (arguments[0] == "mis")
            _MIS = true;
        else
            _nextEventEstimation = false;
    }
    else if (command == "russianroulette")
    {
        if (arguments[0] == "on")
            _russianRoulette = true;
        else
            _russianRoulette = false;
    }
    else if (command == "importancesampling")
    {
        if (arguments[0] == "cosine")
        {
            _importanceSampling = COSINE_SAMPLING;
        }
        else if (arguments[0] == "brdf")
        {
            _importanceSampling = BRDF_SAMPLING;
        }
        else
        {
            _importanceSampling = HEMISPHERE_SAMPLING;
        }
    }
    else if (command == "brdf")
    {
        if (arguments[0] == "ggx")
        {
            _curMaterial.brdf = GGX;
        }
        else if (arguments[0] == "volumetric")
        {
            _curMaterial.brdf = GGX_VOLUMETRIC;
        }
        else
        {
            _curMaterial.brdf = PHONG;
        }
    }
    else if (command == "gamma")
    {
        _gamma = std::stof(arguments[0]);
    }
    else
    {
        std::cerr << "Unknown command in scene file: '" << command << "'" << std::endl;
    }
}

void SceneLoader::loadSceneData(const std::string &filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: '" + filePath + "'");

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream tokenStream(line);

        std::string command;
        tokenStream >> command;

        if (command.size() == 0 || command[0] == '#')
            continue;

        std::vector<std::string> arguments;
        std::string argument;
        while (tokenStream >> argument)
        {
            arguments.push_back(argument);
        }

        executeCommand(command, arguments);
    }
}

void SceneLoader::quadLightToTriangles()
{
    for (auto light : _quadLights)
    {
        int offset = _vertices.size();
        _vertices.push_back(light.a);
        _vertices.push_back(light.a + light.ab);
        _vertices.push_back(light.a + light.ab + light.ac);
        _vertices.push_back(light.a + light.ac);

        _indices.push_back(glm::uvec3(0, 2, 1) + glm::uvec3(offset));
        _indices.push_back(glm::uvec3(0, 3, 2) + glm::uvec3(offset));

        material_t mat;
        mat.ambient = glm::vec3(0);
        mat.diffuse = glm::vec3(1);
        mat.specular = glm::vec3(0);
        mat.emission = light.intensity;
        mat.shininess = 1;
        mat.light = true;

        _triMaterials.push_back(mat);
        _triMaterials.push_back(mat);
    }
}

void SceneLoader::loadEmbreeTriangles(RTCScene embreeScene)
{
    RTCGeometry embreeTriangles = rtcNewGeometry(_embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);

    glm::vec3 *embreeVertices = reinterpret_cast<glm::vec3 *>(rtcSetNewGeometryBuffer(
        embreeTriangles,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT3,
        sizeof(glm::vec3),
        _vertices.size()));
    std::memcpy(embreeVertices, _vertices.data(), _vertices.size() * sizeof(glm::vec3));

    glm::uvec3 *embreeIndices = reinterpret_cast<glm::uvec3 *>(rtcSetNewGeometryBuffer(
        embreeTriangles,
        RTC_BUFFER_TYPE_INDEX,
        0,
        RTC_FORMAT_UINT3,
        sizeof(glm::uvec3),
        _indices.size()));
    std::memcpy(embreeIndices, _indices.data(), _indices.size() * sizeof(glm::uvec3));

    rtcCommitGeometry(embreeTriangles);
    rtcAttachGeometry(embreeScene, embreeTriangles);
    rtcReleaseGeometry(embreeTriangles);
}

void SceneLoader::loadEmbreeSpheres(RTCScene embreeScene)
{
    RTCScene embreeSphereScene = rtcNewScene(_embreeDevice);

    RTCGeometry embreeSphere = rtcNewGeometry(_embreeDevice, RTC_GEOMETRY_TYPE_SPHERE_POINT);

    glm::vec4 *embreeSpherePoint = reinterpret_cast<glm::vec4 *>(rtcSetNewGeometryBuffer(
        embreeSphere,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT4,
        sizeof(glm::vec4),
        1));
    *embreeSpherePoint = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);

    rtcCommitGeometry(embreeSphere);
    rtcAttachGeometry(embreeSphereScene, embreeSphere);
    rtcReleaseGeometry(embreeSphere);
    rtcCommitScene(embreeSphereScene);

    for (glm::mat4 transform : _sphereTransforms)
    {
        RTCGeometry embreeSphereInstance = rtcNewGeometry(_embreeDevice, RTC_GEOMETRY_TYPE_INSTANCE);
        rtcSetGeometryInstancedScene(embreeSphereInstance, embreeSphereScene);
        rtcSetGeometryTimeStepCount(embreeSphereInstance, 1);
        rtcSetGeometryTransform(
            embreeSphereInstance,
            0,
            RTC_FORMAT_FLOAT4X4_COLUMN_MAJOR,
            glm::value_ptr(transform));
        rtcCommitGeometry(embreeSphereInstance);
        rtcAttachGeometry(embreeScene, embreeSphereInstance);
        rtcReleaseGeometry(embreeSphereInstance);
    }

    rtcReleaseScene(embreeSphereScene);
}

RTCScene SceneLoader::createEmbreeScene()
{
    RTCScene embreeScene = rtcNewScene(_embreeDevice);
    loadEmbreeTriangles(embreeScene);
    loadEmbreeSpheres(embreeScene);
    rtcCommitScene(embreeScene);
    return embreeScene;
}

Integrator *SceneLoader::createIntegrator()
{
    Integrator *it;

    if (_integratorType == "raytracer")
        it = new RayTracerIntegrator();
    else if (_integratorType == "direct")
        it = new MonteCarloDirectIntegrator();
    else if (_integratorType == "analyticdirect")
        it = new AnalyticDirectIntegrator();
    else if (_integratorType == "pathtracer")
        it = new PathTracerIntegrator();
    else if (_integratorType == "volumetricpathtracer")
        it = new VolumetricPathTracerIntegrator();
    else
        it = new RayTracerIntegrator();
    return it;
}

Scene *SceneLoader::commitSceneData()
{
    float aspectRatio = static_cast<float>(_imageSize.x) / _imageSize.y;
    glm::vec3 cameraLook = glm::normalize(_cameraLookAt - _cameraOrigin);
    glm::vec3 imagePlaneRight = glm::normalize(glm::cross(cameraLook, _cameraUp));
    glm::vec3 imagePlaneUp = glm::normalize(glm::cross(imagePlaneRight, cameraLook));

    camera_t camera;
    camera.origin = _cameraOrigin;
    camera.imagePlaneTopLeft =
        _cameraOrigin + cameraLook / std::tan(PI * _cameraFieldOfView / 360.0f) + imagePlaneUp - aspectRatio * imagePlaneRight;
    camera.pixelRight = (2.0f * aspectRatio / _imageSize.x) * imagePlaneRight;
    camera.pixelDown = (-2.0f / _imageSize.y) * imagePlaneUp;

    std::vector<glm::mat3> sphereNormalTransforms;
    for (size_t i = 0; i < _sphereTransforms.size(); i++)
    {
        sphereNormalTransforms.push_back(glm::inverseTranspose(glm::mat3(_sphereTransforms[i])));
    }

    quadLightToTriangles();

    Scene *scene = new Scene();
    scene->imageSize = _imageSize;
    scene->maxDepth = _maxDepth;
    scene->outputFileName = _outputFileName;
    scene->camera = camera;
    scene->sphereNormalTransforms = std::move(sphereNormalTransforms);
    scene->sphereMaterials = std::move(_sphereMaterials);
    scene->triMaterials = std::move(_triMaterials);
    scene->directionalLights = std::move(_directionalLights);
    scene->pointLights = std::move(_pointLights);
    scene->quadLights = _quadLights;
    scene->embreeScene = createEmbreeScene();
    scene->integrator = createIntegrator();
    scene->integrator->setScene(scene);
    scene->lightSamples = _lightSamples;
    scene->lightStratify = _lightStratify;
    scene->stratifyGridSize = glm::sqrt(_lightSamples);
    scene->samplesPerPixel = _samplesPerPixel;
    scene->MIS = _MIS;
    scene->nextEventEstimation = _nextEventEstimation;
    scene->russianRoulette = _russianRoulette;
    scene->importanceSampling = _importanceSampling;
    scene->gamma = _gamma;

    // volumetric data
    scene->volumeList = _volumeList;
    scene->defaultVolume = _defaultVolume;

    return scene;
}

void loadScene(
    const std::string &filePath,
    RTCDevice embreeDevice,
    Scene **scene)
{
    SceneLoader sceneLoader(embreeDevice);
    sceneLoader.loadSceneData(filePath);
    *scene = sceneLoader.commitSceneData();
}
