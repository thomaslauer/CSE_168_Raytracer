#include "BRDF.h"

#include <glm/glm.hpp>

thread_local std::default_random_engine BRDF::rng = std::default_random_engine((std::random_device())());

BRDF::BRDF() {
    gen = std::uniform_real_distribution<float>(0.0f, 1.0f);
}
