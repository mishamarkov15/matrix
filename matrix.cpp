//
// Created by Михаил Марков on 11/07/2023.
//

#include "matrix.h"

template<typename vt, size_t height, size_t width>
void matrix<vt, height, width>::allocate_memory() {
    data_ = new vt *[height];
    for (size_t i = 0; i < height; ++i)
        data_[i] = new vt[width];
}

template<typename vt, size_t height, size_t width>
void matrix<vt, height, width>::deallocate_memory() {
    for (size_t i = 0; i < height; ++i)
        delete[] data_[i];
    delete[] data_;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width>::matrix() {
    allocate_memory();
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] = 0;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width>::matrix(const matrix<vt, height, width> &other) {
    allocate_memory();
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] = other.data_[i][j];
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width>::matrix(const vt &num) {
    allocate_memory();
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] = num;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width>::~matrix() {
    deallocate_memory();
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> &matrix<vt, height, width>::operator=(const matrix<vt, height, width> &other) {
    if (this == &other)
        return *this;
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] = other.data_[i][j];
    return *this;
}

template<typename vt, size_t height, size_t width>
template<typename T, size_t H, size_t W>
bool matrix<vt, height, width>::operator==(const matrix<T, H, W> &other) const {
    if (height != H || width != W)
        return false;
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            if (data_[i][j] != other.data_[i][j])
                return false;
    return true;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> &matrix<vt, height, width>::operator+=(const matrix<vt, height, width> &other) {
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] += other.data_[i][j];
    return *this;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> matrix<vt, height, width>::operator+(const matrix<vt, height, width> &other) const {
    matrix<vt, height, width> tmp(*this);
    tmp += other;
    return tmp;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> &matrix<vt, height, width>::operator-=(const matrix<vt, height, width> &other) {
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] -= other.data_[i][j];
    return *this;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> matrix<vt, height, width>::operator-(const matrix<vt, height, width> &other) const {
    matrix<vt, height, width> tmp(*this);
    tmp -= other;
    return tmp;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> &matrix<vt, height, width>::operator*=(const long long int &value) {
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] *= value;
    return *this;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> matrix<vt, height, width>::operator*(const long long int &value) const {
    matrix<vt, height, width> tmp(*this);
    tmp *= value;
    return tmp;
}

template<typename vt, size_t height, size_t width>
std::ostream &operator<<(std::ostream &out, const matrix<vt, height, width> &matrix) {
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j)
            out << matrix.data_[i][j] << " ";
        out << "\n";
    }
    return out;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> &matrix<vt, height, width>::operator*=(const matrix<vt, width, width> &other) {
    matrix<vt, height, width> res;
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            for (size_t k = 0; k < width; ++k)
                res.at(i, j) += data_[i][k] * other.at(k, j);
    *this = res;
    return *this;
}

template<typename vt, size_t height, size_t width>
template<size_t W>
matrix<vt, height, W> matrix<vt, height, width>::operator*(const matrix<vt, width, W> &other) const {
    matrix<vt, height, W> res;
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < W; ++j)
            for (size_t k = 0; k < width; ++k)
                res.at(i, j) += data_[i][k] * other.at(k, j);
    return res;
}

template<typename vt, size_t height, size_t width>
matrix<vt, height, width> matrix<vt, height, width>::operator-() const {
    matrix<vt, height, width> res(*this);
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            res.data_[i][j] = -data_[i][j];
    return res;
}

template<typename vt, size_t height, size_t width>
vt &matrix<vt, height, width>::at(const size_t &i, const size_t &j) {
    return data_[i][j];
}

template<typename vt, size_t height, size_t width>
const vt &matrix<vt, height, width>::at(const size_t &i, const size_t &j) const {
    return data_[i][j];
}

template<typename vt, size_t height, size_t width>
matrix<vt, width, height> matrix<vt, height, width>::transposed() const {
    matrix<vt, width, height> res;
    for (size_t i = 0; i < width; ++i)
        for (size_t j = 0; j < height; ++j)
            res.at(i, j) = data_[j][i];
    return res;
}

template<typename vt, size_t height, size_t width>
vt matrix<vt, height, width>::trace() const {
    if (height != width)
        throw std::out_of_range("matrix height != width");
    vt sum = 0;
    for (size_t i = 0; i < height; ++i)
        sum += data_[i][i];
    return sum;
}

template<typename vt, size_t height, size_t width>
vt matrix<vt, height, width>::det() const {
    if (height != width)
        throw std::out_of_range("Can not calculate determinant when height != width");
    if (height == 1)
        return data_[0][0];
    if (height == 2)
        return data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];
    vt sum = 0;
    for (size_t i = 0; i < height; ++i) {
        sum += ((i + 2) % 2 == 0 ? 1 : -1) * create_minor(*this, 0, i).det() * data_[0][i];
    }

    return sum;
}

template <typename vt, size_t height, size_t width>
matrix<vt, height, width> &matrix<vt, height, width>::operator+=(const vt &value) {
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] += value;
    return *this;
}

template <typename vt, size_t height, size_t width>
matrix<vt, height, width> matrix<vt, height, width>::operator+(const vt &value) const {
    matrix<vt, height, width> res(*this);
    res += value;
    return res;
}

template <typename vt, size_t height, size_t width>
matrix<vt, height, width> &matrix<vt, height, width>::operator-=(const vt &value) {
    for (size_t i = 0; i < height; ++i)
        for (size_t j = 0; j < width; ++j)
            data_[i][j] -= value;
    return *this;
}

template <typename vt, size_t height, size_t width>
matrix<vt, height, width> matrix<vt, height, width>::operator-(const vt &value) const {
    matrix<vt, height, width> res(*this);
    res -= value;
    return res;
}

template <typename vt, size_t height, size_t width>
matrix<vt, height, width> operator+(const vt &value, const matrix<vt, height, width>& m) {
    return m + value;
}

template <typename vt, size_t height, size_t width>
matrix<vt, height, width> operator-(const vt &value, const matrix<vt, height, width>& m) {
    return m - value;
}

template <typename vt, size_t height, size_t width>
matrix<vt, height, width> operator*(const vt &value, const matrix<vt, height, width>& m) {
    return m * value;
}
