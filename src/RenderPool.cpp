#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <iostream>

#include <glm/glm.hpp>

#include "Scene.h"
#include "Integrator.h"

#include "RenderPool.h"

thread_local std::default_random_engine RenderJob::rng = std::default_random_engine((std::random_device())());

RenderJob::RenderJob(glm::uvec2 startPixel, glm::uvec2 windowSize)
    : startPixel(startPixel),
      windowSize(windowSize),
      _result(windowSize.x * windowSize.y)
{
    gen = std::uniform_real_distribution<float>(-0.5f, 0.5f);
}

void RenderJob::render(Scene* scene, Integrator* integrator)
{
    for (size_t wy = 0; wy < windowSize.y; wy++) {
        size_t y = startPixel.y + wy;
        for (size_t wx = 0; wx < windowSize.x; wx++) {
            size_t x = startPixel.x + wx;

            _result[wy * windowSize.x + wx] = glm::vec3(0);

            for (int sampleNum = 0; sampleNum < scene->samplesPerPixel; sampleNum++) {

                float offsetX = 0;
                float offsetY = 0;

                if (sampleNum != 0) {
                    offsetX = gen(rng);
                    offsetY = gen(rng);
                }

                glm::vec3 target =
                    scene->camera.imagePlaneTopLeft
                    + (x + 0.5f + offsetX) * scene->camera.pixelRight
                    + (y + 0.5f + offsetY) * scene->camera.pixelDown;
                glm::vec3 direction = glm::normalize(target - scene->camera.origin);

                glm::vec3 currentResult = integrator->traceRay(scene->camera.origin, direction) / ((float) scene->samplesPerPixel);

                if (!glm::any(glm::isnan(currentResult))) _result[wy * windowSize.x + wx] += currentResult;
            }

            // adjust gamma
            _result[wy * windowSize.x + wx].x = glm::pow(_result[wy * windowSize.x + wx].x, 1.0f / scene->gamma);
            _result[wy * windowSize.x + wx].y = glm::pow(_result[wy * windowSize.x + wx].y, 1.0f / scene->gamma);
            _result[wy * windowSize.x + wx].z = glm::pow(_result[wy * windowSize.x + wx].z, 1.0f / scene->gamma);

        }
    }
}

std::vector<glm::vec3> RenderJob::getResult()
{
    return std::move(_result);
}

RenderPool::RenderPool(Scene* scene, Integrator* integrator, int numThreads, std::vector<RenderJob*>& jobs)
    : _scene(scene), _integrator(integrator), _nextJob(0), _jobQueue(jobs)
{
    for (int i = 0; i < numThreads; i++) {
        _threads.push_back(std::thread(threadMain, this));
    }
}

RenderPool::~RenderPool()
{
    for (std::thread& thread : _threads) {
        thread.join();
    }
}

void RenderPool::getCompletedJobs(std::vector<RenderJob*>& completedJobs)
{
    {
        std::unique_lock<std::mutex> lock(_mutex);

        _condition.wait(lock, [this]{ return _completedJobs.size() > 0; });
        completedJobs = std::move(_completedJobs);
    }
}

void RenderPool::threadMain(RenderPool* pool)
{
    while (true) {

        size_t jobIndex;
        {
            std::unique_lock<std::mutex> lock(pool->_mutex);

            if (pool->_nextJob >= pool->_jobQueue.size()) break;

            jobIndex = pool->_nextJob;
            pool->_nextJob++;
        }

        pool->_jobQueue[jobIndex]->render(pool->_scene, pool->_integrator);

        {
            std::unique_lock<std::mutex> lock(pool->_mutex);

            pool->_completedJobs.push_back(pool->_jobQueue[jobIndex]);
            pool->_condition.notify_all();
        }
    }
}
