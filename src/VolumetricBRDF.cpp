#include "BRDF.h"

glm::vec3 VolumetricBRDF::brdf(
    glm::vec3 normal, 
    glm::vec3 w_in, 
    glm::vec3 w_out, 
    material_t material)
{
    return glm::vec3(0);
}

glm::vec3 VolumetricBRDF::importanceSample(
    glm::vec3 normal, 
    glm::vec3 w_out, 
    material_t material, 
    float& pdfNormalization)
{
}

float VolumetricBRDF::pdf(
    glm::vec3 normal, 
    glm::vec3 w_in, 
    glm::vec3 w_out, 
    material_t material)
{
}