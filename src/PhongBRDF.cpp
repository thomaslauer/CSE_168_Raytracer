#include "BRDF.h"

#include "Constants.h"
#include "MathUtils.h"

glm::vec3 PhongBRDF::brdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) {
    glm::vec3 reflection = 2 * glm::dot(normal, w_out) * normal - w_out;

    glm::vec3 diffuse = material.diffuse * INV_PI;
    glm::vec3 specular = material.specular * (material.shininess + 2) / (2 * PI) * glm::pow(glm::dot(reflection, w_in), material.shininess);

    return diffuse + specular;
}

glm::vec3 PhongBRDF::importanceSample(glm::vec3 normal, glm::vec3 w_out, material_t material, float& pdfNormalization) {
    float epsilon1 = gen(rng);
    float epsilon2 = gen(rng);

    float theta = 0;
    float phi = 0;

    float k_s = averageVector(material.specular);
    float k_d = averageVector(material.diffuse);

    float t = k_s / (k_s + k_d);

    glm::vec3 samplingSpaceCenter = normal;
    glm::vec3 reflection = (2 * glm::dot(normal, w_out) * normal - w_out);

    if (gen(rng) < t)
    {
        // specular pdf
        theta = glm::acos(glm::pow(epsilon1, 1.0 / (1.0 + material.shininess)));
        phi = TWO_PI * epsilon2;

        samplingSpaceCenter = reflection;
    }
    else
    {
        // diffuse pdf
        theta = glm::acos(glm::sqrt(epsilon1));
        phi = TWO_PI * epsilon2;
    }

    glm::vec3 w_in = sphereCoordsToVector(theta, phi, samplingSpaceCenter);

    // calculate pdf normalization
    float cosTerm = glm::max(0.0f, glm::dot(w_in, normal));

    float diffuse = (1 - t) * cosTerm / PI;
    float specular = t * (material.shininess + 1) / TWO_PI * glm::pow(glm::max(0.0f, glm::dot(reflection, w_in)), material.shininess);
    //return cosTerm / (diffuseNormalization + specularNormalization);

    pdfNormalization = diffuse + specular;


    return w_in;
}

/*
float PhongBRDF::pdf(glm::vec3 normal, glm::vec3 w_in, glm::vec3 w_out, material_t material) {

    float k_s = averageVector(material.specular);
    float k_d = averageVector(material.diffuse);

    float t = k_s / (k_s + k_d);
    glm::vec3 reflection = (2 * glm::dot(normal, w_out) * normal - w_out);
    //glm::vec3 reflection = glm::reflect(w_out, normal);

    float cosTerm = glm::max(0.0f, glm::dot(w_in, normal));

    float diffuse = (1 - t) * cosTerm / PI;
    float specular = t * (material.shininess + 1) / TWO_PI * glm::pow(glm::max(0.0f, glm::dot(reflection, w_in)), material.shininess);
    //return cosTerm / (diffuseNormalization + specularNormalization);

    return diffuse + specular;
}
*/