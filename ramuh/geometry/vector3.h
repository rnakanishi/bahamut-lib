#ifndef __RAMUH_VECTOR3_H__
#define __RAMUH_VECTOR3_H__

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

  Vector3<T> operator+(Vector3<T> v) {
    return Vector3<T>(_x + v.x(), _y + v.y(), _z + v.z());
  }

  Vector3<T> operator-(Vector3<T> v) {
    return Vector3<T>(_x - v.x(), _y - v.y(), _z - v.z());
  }

  Vector3<T> operator-(T alpha) {
    return Vector3<T>(_x - alpha, _y - alpha, _z - alpha);
  }

  Vector3<double> operator/(double s) {
    return Vector3<double>(_x / s, _y / s, _z / s);
  }

  template <typename type>
  Vector3<double> operator*(const Vector3<type> &v) const {
    return Vector3<double>(_x * v.x(), _y * v.y(), _z * v.z());
  }

  Vector3<double> operator*(const double &s) const {
    return Vector3<double>(_x * s, _y * s, _z * s);
  }

  double dot(const Vector3<T> &v) {
    return _x * v.x() + _y * v.y() + _z * v.z();
  }

  void set(const Vector3<T> &v) {
    _x = v.x();
    _y = v.y();
    _z = v.z();
  }

  void set(const T &x, const T &y, const T &z) {
    _x = x;
    _y = y;
    _z = z;
  }

  friend std::ostream &operator<<(std::ostream &os, const Vector3<T> &v) {
    os << "[" << v.x() << ", " << v.y() << ", " << v.z() << "]";
    return os;
  }

protected:
  T _x, _y, _z;
};

typedef Vector3<int> Vector3i;
typedef Vector3<double> Vector3d;

} // namespace Ramuh

#endif