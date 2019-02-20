#include <iostream>

namespace Ramuh {
template <typename T> class Vector3 {

public:
  Vector3() {}
  Vector3(T val) : _x(val), _y(val), _z(val) {}
  Vector3(T x, T y, T z) : _x(x), _y(y), _z(z) {}

  T x() const { return _x; }
  T y() const { return _y; }
  T z() const { return _z; }

  void x(T newX) { _x = newX; }
  void y(T newY) { _y = newY; }
  void z(T newZ) { _z = newZ; }

  template <typename type>
  Vector3<double> operator/(const Vector3<type> &v) const {
    return Vector3<double>(_x / v.x(), _y / v.y(), _z / v.z());
  }

  friend std::ostream &operator<<(std::ostream &os, const Vector3<T> &v) {
    os << "[" << v.x() << ", " << v.y() << ", " << v.z() << "]";
    return os;
  }

protected:
  T _x, _y, _z;
};

#define Vector3d Vector3<double>
#define Vector3i Vector3<int>

} // namespace Ramuh
