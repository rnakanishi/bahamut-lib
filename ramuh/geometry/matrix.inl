
template <typename T> Matrix3<T>::Matrix3() {}

template <typename T> Matrix3<T>::Matrix3(Vector3i size) {

  _data.resize(size.x() + 1);
  for (auto &row : _data) {
    row.resize(size.y());
    for (auto &depth : row)
      depth.resize(size.z());
  }
}

template <typename T>
std::vector<std::vector<T>> Matrix3<T>::operator[](int i) {
  return _data[i];
}
