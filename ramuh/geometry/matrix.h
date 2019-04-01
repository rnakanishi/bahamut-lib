#ifndef __RAMUH_MATRIX_H__
#define __RAMUH_MATRIX_H__

#include <geometry/vector3.h>
#include <vector>

namespace Ramuh {

template <typename T> class Matrix3 {
public:
  Matrix3();

  Matrix3(Vector3i size);

  ~Matrix3();

  void changeSize(Vector3i size);

  void changeSize(Vector3i size, T value);

  std::vector<std::vector<T>> &operator[](int i);

protected:
  std::vector<std::vector<std::vector<T>>> _data;
};

// #include <geometry/matrix.inl>
//*****************************************************************************
//*****************************************************************************

template <typename T> Matrix3<T>::Matrix3() {}

template <typename T> Matrix3<T>::~Matrix3() {}

template <typename T> Matrix3<T>::Matrix3(Vector3i size) {
  this->changeSize(size);
}

template <typename T> void Matrix3<T>::changeSize(Vector3i size) {
  _data.resize(size.x());
  for (auto &row : _data) {
    row.resize(size.y());
    for (auto &depth : row)
      depth.resize(size.z());
  }
}

template <typename T> void Matrix3<T>::changeSize(Vector3i size, T value) {
  _data.resize(size.x());
  for (auto &row : _data) {
    row.resize(size.y());
    for (auto &depth : row)
      depth.resize(size.z(), value);
  }
}

template <typename T>
std::vector<std::vector<T>> &Matrix3<T>::operator[](int i) {
  return _data[i];
}

} // namespace Ramuh
#endif
