#include <graphics/camera.hpp>

namespace Garuda {
Camera::Camera();

glm::vec3 Camera::getPosition() { return position; }

} // namespace Garuda