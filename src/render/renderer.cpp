#include "renderer.h"
#include <algorithm>
#include <algorithm> // std::fill
#include <cmath>
#include <functional>
#include <glm/common.hpp>
#include <glm/gtx/component_wise.hpp>
#include <iostream>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tuple>
#include <math.h>

namespace render {

// The renderer is passed a pointer to the volume, gradinet volume, camera and an initial renderConfig.
// The camera being pointed to may change each frame (when the user interacts). When the renderConfig
// changes the setConfig function is called with the updated render config. This gives the Renderer an
// opportunity to resize the framebuffer.
Renderer::Renderer(
    const volume::Volume* pVolume,
    const volume::GradientVolume* pGradientVolume,
    const render::RayTraceCamera* pCamera,
    const RenderConfig& initialConfig)
    : m_pVolume(pVolume)
    , m_pGradientVolume(pGradientVolume)
    , m_pCamera(pCamera)
    , m_config(initialConfig)
{
    resizeImage(initialConfig.renderResolution);
}

// Set a new render config if the user changed the settings.
void Renderer::setConfig(const RenderConfig& config)
{
    if (config.renderResolution != m_config.renderResolution)
        resizeImage(config.renderResolution);

    m_config = config;
}

// Resize the framebuffer and fill it with black pixels.
void Renderer::resizeImage(const glm::ivec2& resolution)
{
    m_frameBuffer.resize(size_t(resolution.x) * size_t(resolution.y), glm::vec4(0.0f));
}

// Clear the framebuffer by setting all pixels to black.
void Renderer::resetImage()
{
    std::fill(std::begin(m_frameBuffer), std::end(m_frameBuffer), glm::vec4(0.0f));
}

// Return a VIEW into the framebuffer. This view is merely a reference to the m_frameBuffer member variable.
// This does NOT make a copy of the framebuffer.
gsl::span<const glm::vec4> Renderer::frameBuffer() const
{
    return m_frameBuffer;
}

// Main render function. It computes an image according to the current renderMode.
// Multithreading is enabled in Release/RelWithDebInfo modes. In Debug mode multithreading is disabled to make debugging easier.
void Renderer::render()
{
    resetImage();

    static constexpr float sampleStep = 1.0f;
    const glm::vec3 planeNormal = -glm::normalize(m_pCamera->forward());
    const glm::vec3 volumeCenter = glm::vec3(m_pVolume->dims()) / 2.0f;
    const Bounds bounds { glm::vec3(0.0f), glm::vec3(m_pVolume->dims() - glm::ivec3(1)) };

    // 0 = sequential (single-core), 1 = TBB (multi-core)
#ifdef NDEBUG
    // If NOT in debug mode then enable parallelism using the TBB library (Intel Threaded Building Blocks).
#define PARALLELISM 1
#else
    // Disable multi threading in debug mode.
#define PARALLELISM 0
#endif

#if PARALLELISM == 0
    // Regular (single threaded) for loops.
    for (int x = 0; x < m_config.renderResolution.x; x++) {
        for (int y = 0; y < m_config.renderResolution.y; y++) {
#else
    // Parallel for loop (in 2 dimensions) that subdivides the screen into tiles.
    const tbb::blocked_range2d<int> screenRange { 0, m_config.renderResolution.y, 0, m_config.renderResolution.x };
        tbb::parallel_for(screenRange, [&](tbb::blocked_range2d<int> localRange) {
        // Loop over the pixels in a tile. This function is called on multiple threads at the same time.
        for (int y = std::begin(localRange.rows()); y != std::end(localRange.rows()); y++) {
            for (int x = std::begin(localRange.cols()); x != std::end(localRange.cols()); x++) {
#endif
            // Compute a ray for the current pixel.
            const glm::vec2 pixelPos = glm::vec2(x, y) / glm::vec2(m_config.renderResolution);
            Ray ray = m_pCamera->generateRay(pixelPos * 2.0f - 1.0f);

            // Compute where the ray enters and exists the volume.
            // If the ray misses the volume then we continue to the next pixel.
            if (!instersectRayVolumeBounds(ray, bounds))
                continue;

            // Get a color for the current pixel according to the current render mode.
            glm::vec4 color {};
            switch (m_config.renderMode) {
            case RenderMode::RenderSlicer: {
                color = traceRaySlice(ray, volumeCenter, planeNormal);
                break;
            }
            case RenderMode::RenderMIP: {
                color = traceRayMIP(ray, sampleStep);
                break;
            }
            case RenderMode::RenderComposite: {
                color = traceRayComposite(ray, sampleStep);
                break;
            }
            case RenderMode::RenderIso: {
                color = traceRayISO(ray, sampleStep);
                break;
            }
            case RenderMode::RenderTF2D: {
                color = traceRayTF2D(ray, sampleStep);
                break;
            }
            };
            // Write the resulting color to the screen.
            fillColor(x, y, color);

#if PARALLELISM == 1
        }
    }
});
#else
            }
        }
#endif
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// This function generates a view alongside a plane perpendicular to the camera through the center of the volume
//  using the slicing technique.
glm::vec4 Renderer::traceRaySlice(const Ray& ray, const glm::vec3& volumeCenter, const glm::vec3& planeNormal) const
{
    const float t = glm::dot(volumeCenter - ray.origin, planeNormal) / glm::dot(ray.direction, planeNormal);
    const glm::vec3 samplePos = ray.origin + ray.direction * t;
    const float val = m_pVolume->getSampleInterpolate(samplePos);
    return glm::vec4(glm::vec3(std::max(val / m_pVolume->maximum(), 0.0f)), 1.f);
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// Function that implements maximum-intensity-projection (MIP) raycasting.
// It returns the color assigned to a ray/pixel given it's origin, direction and the distances
// at which it enters/exits the volume (ray.tmin & ray.tmax respectively).
// The ray must be sampled with a distance defined by the sampleStep
glm::vec4 Renderer::traceRayMIP(const Ray& ray, float sampleStep) const
{
    float maxVal = 0.0f;

    // Incrementing samplePos directly instead of recomputing it each frame gives a measureable speed-up.
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        maxVal = std::max(val, maxVal);
    }

    // Normalize the result to a range of [0 to mpVolume->maximum()].
    return glm::vec4(glm::vec3(maxVal) / m_pVolume->maximum(), 1.0f);
}

// ======= TODO: IMPLEMENT ========
// This function should find the position where the ray intersects with the volume's isosurface.
// If volume shading is DISABLED then simply return the isoColor.
// If volume shading is ENABLED then return the phong-shaded color at that location using the local gradient (from m_pGradientVolume).
//   Use the camera position (m_pCamera->position()) as the light position.
// Use the bisectionAccuracy function (to be implemented) to get a more precise isosurface location between two steps.
glm::vec4 Renderer::traceRayISO(const Ray& ray, float sampleStep) const
{
    static constexpr glm::vec3 isoColor { 0.8f, 0.8f, 0.2f };
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        if (val >= m_config.isoValue) {
            const float optimalT = bisectionAccuracy(ray, t - 1, t, m_config.isoValue);
            samplePos = ray.origin + optimalT * sampleStep * ray.direction;
            if (m_config.volumeShading) {
                const volume::GradientVoxel gradient = m_pGradientVolume->getGradientInterpolate(samplePos);
                const glm::vec3 normalizedCameraPosition = glm::normalize(m_pCamera->position());
                return glm::vec4(computePhongShading(glm::vec3(isoColor), gradient, normalizedCameraPosition, normalizedCameraPosition), 1);
            } else {
                return glm::vec4(isoColor, 1.0f);
            }
        }
    }
    return glm::vec4(glm::vec3(0.0f), 0.0f);
}

// ======= TODO: IMPLEMENT ========
// Given that the iso value lies somewhere between t0 and t1, find a t for which the value
// closely matches the iso value (less than 0.01 difference). Add a limit to the number of
// iterations such that it does not get stuck in degerate cases.
float Renderer::bisectionAccuracy(const Ray& ray, float t0, float t1, float isoValue) const
{
    // We use i to track the number of iterations and terminate the loop when i equals maxIterations
    int i = 0;
    int maxIterations = 100;
    
    float isoValueDifference = 100;
    float currentT = t1;
    while (isoValueDifference > 0.01 && i < maxIterations) {
        // every iteration step we cut the interval in half
        float newPoint = (t0 + t1) / 2;
        const glm::vec3 samplePos = ray.origin + newPoint * ray.direction;
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        isoValueDifference = abs(val - isoValue);
        // depending on where the value lies related to the isovalue we choose the left or right half
        if (val < isoValue) {
            t0 = newPoint;
        } else {
            t1 = newPoint;
        }
        currentT = newPoint;
        i++;
    }
    return currentT;
}

// ======= TODO: IMPLEMENT ========
// Compute Phong Shading given the voxel color (material color), the gradient, the light vector and view vector.
// You can find out more about the Phong shading model at:
// https://en.wikipedia.org/wiki/Phong_reflection_model
//
// Use the given color for the ambient/specular/diffuse (you are allowed to scale these constants by a scalar value).
// You are free to choose any specular power that you'd like.
glm::vec3 Renderer::computePhongShading(const glm::vec3& color, const volume::GradientVoxel& gradient, const glm::vec3& L, const glm::vec3& V)
{
    glm::vec3 normalizedGradient = glm::normalize(gradient.dir);
    glm::vec3 reflection = glm::normalize(L - 2 * glm::dot(L, normalizedGradient) * normalizedGradient);

    glm::vec3 whiteLight = glm::vec3(1.0f, 1.0f, 1.0f);

    // These are the Phong shading parameters givin in the exercise
    float ka = 0.1;
    float kd = 0.7;
    float ks = 0.2;
    float alpha = 100;

    // We can directly calculate the ambient term
    glm::vec3 ambientTerm = ka * whiteLight * color;

    glm::vec3 diffuseTerm = glm::vec3(0.0f);
    glm::vec3 specularTerm = glm::vec3(0.0f); 
    float diffuseCosine = glm::dot(-L, normalizedGradient);
    // The diffuse and specular terms are only nonzero in case the diffuseCosine is positive.
    if (diffuseCosine > 0) {
        diffuseTerm = kd * whiteLight * color * diffuseCosine;
        float specularCosine = glm::dot(-reflection, V);
        // The specular term is only nonzero when the specularCosine is also positive
        if (specularCosine > 0) {
            specularTerm = ks * whiteLight * color * pow(specularCosine, alpha);
        }
    }
   
    glm::vec3 res = ambientTerm + diffuseTerm + specularTerm; 

    return res;
}

// ======= TODO: IMPLEMENT ========
// In this function, implement 1D transfer function raycasting.
// Use getTFValue to compute the color for a given volume value according to the 1D transfer function.
glm::vec4 Renderer::traceRayComposite(const Ray& ray, float sampleStep) const
{
    //const float earlyTerminationDelta = 0.01;
    
    // we use these two variables to increment over the ray
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    // these two values are used to hold the accumulated values for the RayComposite
    glm::vec3 accumulatedColor = glm::vec3(0.0f);
    float accumulatedOpacity = 0.0f;

    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        const glm::vec4 tfValue = getTFValue(val);
        glm::vec3 color = glm::vec3(tfValue);

        // we apply Phong Shading when necessary
        if (m_config.volumeShading) {
            const glm::vec3 normalizedCameraPosition = glm::normalize(m_pCamera->position());
            const volume::GradientVoxel gradient = m_pGradientVolume->getGradientInterpolate(samplePos);
            color = computePhongShading(color, gradient, normalizedCameraPosition, normalizedCameraPosition);
        }

        // we increment the values using front-to-back compositing
        accumulatedColor += (1 - accumulatedOpacity) * color * tfValue[3];
        accumulatedOpacity += (1 - accumulatedOpacity) * tfValue[3];

        //if (1.0f - accumulatedOpacity < earlyTerminationDelta)
        //    break;
    }

    return glm::vec4(accumulatedColor, accumulatedOpacity);
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// Looks up the color+opacity corresponding to the given volume value from the 1D tranfer function LUT (m_config.tfColorMap).
// The value will initially range from (m_config.tfColorMapIndexStart) to (m_config.tfColorMapIndexStart + m_config.tfColorMapIndexRange) .
glm::vec4 Renderer::getTFValue(float val) const
{
    // Map value from [m_config.tfColorMapIndexStart, m_config.tfColorMapIndexStart + m_config.tfColorMapIndexRange) to [0, 1) .
    const float range01 = (val - m_config.tfColorMapIndexStart) / m_config.tfColorMapIndexRange;
    const size_t i = std::min(static_cast<size_t>(range01 * static_cast<float>(m_config.tfColorMap.size())), m_config.tfColorMap.size() - 1);
    return m_config.tfColorMap[i];
}

// ======= TODO: IMPLEMENT ========
// In this function, implement 2D transfer function raycasting.
// Use the getTF2DOpacity function that you implemented to compute the opacity according to the 2D transfer function.
glm::vec4 Renderer::traceRayTF2D(const Ray& ray, float sampleStep) const
{
    // this function is almost identical to traceRayComposite the only difference is the way we ge the tfValue
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    glm::vec3 accumulatedColor = glm::vec3(0.0f);
    float accumulatedOpacity = 0.0f;
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        const volume::GradientVoxel gradient = m_pGradientVolume->getGradientInterpolate(samplePos);

        // here we get the tfValue using the 2d transfer function
        const glm::vec4 tfValue = glm::vec4(glm::vec3(m_config.TF2DColor), m_config.TF2DColor[3]*getTF2DOpacity(val, gradient.magnitude));
        
        glm::vec3 color = glm::vec3(tfValue);
        if (m_config.volumeShading) {
            const glm::vec3 normalizedCameraPosition = glm::normalize(m_pCamera->position());
            color = computePhongShading(color, gradient, normalizedCameraPosition, normalizedCameraPosition);
        }

        accumulatedColor += (1 - accumulatedOpacity) * color * tfValue[3];
        accumulatedOpacity += (1 - accumulatedOpacity) * tfValue[3];
    }
    return glm::vec4(accumulatedColor, accumulatedOpacity);
}

// ======= TODO: IMPLEMENT ========
// This function should return an opacity value for the given intensity and gradient according to the 2D transfer function.
// Calculate whether the values are within the radius/intensity triangle defined in the 2D transfer function widget.
// If so: return a tent weighting as described in the assignment
// Otherwise: return 0.0f
//
// The 2D transfer function settings can be accessed through m_config.TF2DIntensity and m_config.TF2DRadius.
float Renderer::getTF2DOpacity(float intensity, float gradientMagnitude) const
{
    if (abs(intensity - m_config.TF2DIntensity) >= m_config.TF2DRadius){
        return 0.0f;
    }
    // Calculate the slope for the left side triangle
    float a1 = -m_pGradientVolume->maxMagnitude() / m_config.TF2DRadius;
    // Calculate the slope for the right side of the triangle
    float a2 = -a1;
    // Calculate the intercepts for the left and right side
    float b1 = -m_config.TF2DIntensity * a1;
    float b2 = -m_config.TF2DIntensity * a2;
    // Calculate the x coordinate of the intersection between the line y=gradientMagnitude and the two sides
    float x1 = (gradientMagnitude - b1) / a1;
    float x2 = (gradientMagnitude - b2) / a2;
    if (intensity < x1 || intensity > x2) {
        return 0.0f;
    }
    return 1 - abs(intensity - m_config.TF2DIntensity) / abs(x2 - m_config.TF2DIntensity);
}

// This function computes if a ray intersects with the axis-aligned bounding box around the volume.
// If the ray intersects then tmin/tmax are set to the distance at which the ray hits/exists the
// volume and true is returned. If the ray misses the volume the the function returns false.
//
// If you are interested you can learn about it at.
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
bool Renderer::instersectRayVolumeBounds(Ray& ray, const Bounds& bounds) const
{
    const glm::vec3 invDir = 1.0f / ray.direction;
    const glm::bvec3 sign = glm::lessThan(invDir, glm::vec3(0.0f));

    float tmin = (bounds.lowerUpper[sign[0]].x - ray.origin.x) * invDir.x;
    float tmax = (bounds.lowerUpper[!sign[0]].x - ray.origin.x) * invDir.x;
    const float tymin = (bounds.lowerUpper[sign[1]].y - ray.origin.y) * invDir.y;
    const float tymax = (bounds.lowerUpper[!sign[1]].y - ray.origin.y) * invDir.y;

    if ((tmin > tymax) || (tymin > tmax))
        return false;
    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    const float tzmin = (bounds.lowerUpper[sign[2]].z - ray.origin.z) * invDir.z;
    const float tzmax = (bounds.lowerUpper[!sign[2]].z - ray.origin.z) * invDir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    ray.tmin = std::max(tmin, tzmin);
    ray.tmax = std::min(tmax, tzmax);
    return true;
}

// This function inserts a color into the framebuffer at position x,y
void Renderer::fillColor(int x, int y, const glm::vec4& color)
{
    const size_t index = static_cast<size_t>(m_config.renderResolution.x * y + x);
    m_frameBuffer[index] = color;
}
}