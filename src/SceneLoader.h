#pragma once

#include <glm/glm.hpp>
#include <embree3/rtcore.h>

#include "Scene.h"

void loadScene(
    const std::string& filePath,
    RTCDevice device,
    Scene** scene);
