#ifndef __RAMUH_VECTOR2_H__
#define __RAMUH_VECTOR2_H__

#include <iostream>
#include <cmath>

namespace Ramuh {
template <typename T> class Vector2 {

public:
  Vector2() {}
  Vector2(T val) : _x(val), _y(val) {}
  Vector2(T x, T y) : _x(x), _y(y) {}

  T x() const { return _x; }
  T y() const { return _y; }

  void x(T newX) { _x = newX; }
  void y(T newY) { _y = newY; }

  template <typename type>
  Vector2<double> operator/(const Vector2<type> &v) const {
    return Vector2<double>(_x / v.x(), _y / v.y());
  }

  Vector2<T> operator+(Vector2<T> v) {
    return Vector2<T>(_x + v.x(), _y + v.y());
  }

  Vector2<T> operator-(Vector2<T> v) {
    return Vector2<T>(_x - v.x(), _y - v.y());
  }

  Vector2<T> operator-(T alpha) { return Vector2<T>(_x - alpha, _y - alpha); }

  Vector2<double> operator/(double s) {
    return Vector2<double>(_x / s, _y / s);
  }

  bool operator>(const Vector2<T> v) {
    if (_x <= v.x())
      return false;
    if (_y <= v.y())
      return false;
    return true;
  }

  bool operator>=(const Vector2<T> v) {
    if (_x < v.x())
      return false;
    if (_y < v.y())
      return false;
    return true;
  }

  bool operator<(const Vector2<T> &v) {
    if (_x >= v.x())
      return false;
    if (_y >= v.y())
      return false;
    return true;
  }

  bool operator<=(const Vector2<T> &v) {
    if (_x > v.x())
      return false;
    if (_y > v.y())
      return false;
    return true;
  }

  bool operator==(const Vector2<T> v) { return (_x == v.x() && _y == v.y()); }

  template <typename type>
  Vector2<double> operator*(const Vector2<type> &v) const {
    return Vector2<double>(_x * v.x(), _y * v.y());
  }

  Vector2<double> operator*(const double &s) const {
    return Vector2<double>(_x * s, _y * s);
  }

  double dot(const Vector2<T> &v) { return _x * v.x() + _y * v.y(); }

  double min() { return std::min(_x, _y); }

  double length() { return std::sqrt(_x * _x + _y * _y); }

  double sqrtLength() { return _x * _x + _y * _y; }

  Vector2<T> abs() { return Vector2<T>(std::abs(_x), std::abs(_y)); }

  void set(const Vector2<T> &v) {
    _x = v.x();
    _y = v.y();
  }

  void set(const T &x, const T &y) {
    _x = x;
    _y = y;
  }

  friend std::ostream &operator<<(std::ostream &os, const Vector2<T> &v) {
    os << "[" << v.x() << ", " << v.y() << "]";
    return os;
  }

protected:
  T _x, _y;
};

template <typename T> class Vector2CompareLess {
public:
  bool operator()(const Vector2<T> *lhs, const Vector2<T> *rhs) {
    return *lhs < *rhs;
  }
};
typedef Vector2<int> Vector2i;
typedef Vector2<double> Vector2d;

} // namespace Ramuh

#endif