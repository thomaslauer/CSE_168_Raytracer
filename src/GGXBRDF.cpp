#include "BRDF.h"

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>

#include "Constants.h"

glm::vec3 GGXBRDF::brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) {
    glm::vec3 halfVector = glm::normalize(glm::normalize(w_in) + glm::normalize(w_out));
    float halfAngle = glm::acos(glm::dot(halfVector, normal));

    float D = microfacetDistribution(halfAngle, material);
    float G = microfacetSelfShadowing(normal, w_in, material) * microfacetSelfShadowing(normal, w_out, material);
    glm::vec3 F = fresnel(w_in, halfVector, material);

    /*
    if (glm::isnan(D)) std::cout << "D IS NAN" << std::endl;
    if (glm::isnan(G)) std::cout << "G IS NAN" << std::endl;
    if (glm::any(glm::isnan(F))) std::cout << "F IS NAN" << std::endl;
    */

    /*
    if (glm::isnan(halfAngle)) {
        std::cout << "half angle IS NAN" << std::endl;
        std::cout << glm::to_string(w_out) << " " << glm::to_string(w_in) << std::endl;
    }
    */



    float normalization = 4 * glm::dot(w_in, normal) * glm::dot(w_out, normal);

    return material.diffuse / PI + D * G * F / normalization;
}

glm::vec3 GGXBRDF::importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material) {
    float epsilon1 = gen(rng);
    float epsilon2 = gen(rng);

    float theta = 0;
    float phi = 0;

    float k_s = averageVector(material.specular);
    float k_d = averageVector(material.diffuse);

    float t = glm::max(0.25f, k_s / (k_s + k_d));

    if (gen(rng) < t)
    {
        // ggx pdf
        theta = glm::atan(material.roughness * glm::sqrt(epsilon1), glm::sqrt(1 - epsilon1));
        phi = TWO_PI * epsilon2;

        glm::vec3 halfVector = sphereCoordsToVector(theta, phi, normal);
        glm::vec3 w_in = glm::normalize(2.0f * halfVector - w_out);

        return w_in;
    }
    else
    {
        // diffuse pdf
        theta = glm::acos(glm::sqrt(epsilon1));
        phi = TWO_PI * epsilon2;
        glm::vec3 w_in = sphereCoordsToVector(theta, phi, normal);
        return w_in;
    }
}

float GGXBRDF::pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) {
    float k_s = averageVector(material.specular);
    float k_d = averageVector(material.diffuse);
    float t = glm::max(0.25f, k_s / (k_s + k_d));
    glm::vec3 halfVector = glm::normalize(glm::normalize(w_in) + glm::normalize(w_out));
    float halfAngle = glm::acos(glm::dot(halfVector, normal));

    float diffuseTerm = (1-t) * glm::dot(normal, w_in) / PI;
    float ggxTerm = t * microfacetDistribution(halfAngle, material) * glm::dot(normal, halfVector) / (4 * glm::dot(halfVector, w_in));

    return glm::dot(normal, w_in) / (diffuseTerm + ggxTerm);
}

float GGXBRDF::microfacetDistribution(float halfAngle, material_t material) {
    float a_squared = material.roughness * material.roughness;

    float denominator = PI * glm::pow(glm::cos(halfAngle), 4) * glm::pow(a_squared + glm::pow(glm::tan(halfAngle), 2), 2);
    return a_squared / denominator;
}

float GGXBRDF::microfacetSelfShadowing(glm::vec3 normal, glm::vec3 view, material_t material) {
    if (glm::dot(view, normal) <= 0)
    {
        return 0;
    }

    float thetaV = glm::acos(glm::dot(view, normal));

    return 2.0f / (1.0f + glm::sqrt(1.0f + glm::pow(material.roughness, 2.0f) * glm::pow(glm::tan(thetaV), 2.0f)));
}

glm::vec3 GGXBRDF::fresnel(glm::vec3 w_in, glm::vec3 halfVector, material_t material) {
    return material.specular + (glm::vec3(1) - material.specular) * glm::pow(1.0f - glm::dot(w_in, halfVector), 5.0f);
}