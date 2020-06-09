#include <limits>

#include <glm/glm.hpp>
#include <embree3/rtcore.h>
#include <iostream>

#include "Scene.h"

bool Scene::castRay(
    glm::vec3 origin,
    glm::vec3 direction,
    glm::vec3 *hitPosition,
    glm::vec3 *hitNormal,
    material_t *hitMaterial) const
{
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    RTCRayHit rayHit;
    rayHit.ray.org_x = origin.x;
    rayHit.ray.org_y = origin.y;
    rayHit.ray.org_z = origin.z;
    rayHit.ray.dir_x = direction.x;
    rayHit.ray.dir_y = direction.y;
    rayHit.ray.dir_z = direction.z;
    rayHit.ray.tnear = 0.0001f;
    rayHit.ray.tfar = std::numeric_limits<float>::infinity();
    rayHit.ray.mask = 0;
    rayHit.ray.flags = 0;
    rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

    rtcIntersect1(embreeScene, &context, &rayHit);

    if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
    {
        *hitPosition = origin + direction * rayHit.ray.tfar;
        *hitNormal = glm::normalize(glm::vec3(rayHit.hit.Ng_x, rayHit.hit.Ng_y, rayHit.hit.Ng_z));
        if (rayHit.hit.instID[0] == RTC_INVALID_GEOMETRY_ID)
        {
            // hit a triangle
            *hitMaterial = triMaterials[rayHit.hit.primID];

            if (hitMaterial->triangleData.interpolate)
            {
                //std::cout << rayHit.hit.u << " " << rayHit.hit.v << std::endl;
                float u = rayHit.hit.u;
                float v = rayHit.hit.v;

                glm::vec3 a = hitMaterial->triangleData.normals[0];
                glm::vec3 b = hitMaterial->triangleData.normals[1];
                glm::vec3 c = hitMaterial->triangleData.normals[2];

                *hitNormal = (1 - u - v) * a + u * b + v * c;
            }
        }
        else
        {
            int sphereIndex = rayHit.hit.instID[0] - 1;
            *hitNormal = glm::normalize(sphereNormalTransforms[sphereIndex] * (*hitNormal));
            *hitMaterial = sphereMaterials[sphereIndex];
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool Scene::castOcclusionRay(glm::vec3 origin, glm::vec3 direction, float maxDistance) const
{
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    RTCRay ray;
    ray.org_x = origin.x;
    ray.org_y = origin.y;
    ray.org_z = origin.z;
    ray.dir_x = direction.x;
    ray.dir_y = direction.y;
    ray.dir_z = direction.z;
    ray.tnear = 0.0001f;
    ray.tfar = maxDistance - 0.0001f;
    ray.mask = 0;
    ray.flags = 0;

    rtcOccluded1(embreeScene, &context, &ray);

    return (ray.tfar < 0.0f);
}
