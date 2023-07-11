#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <iostream>
#include <stdexcept>

template<typename vt, size_t height, size_t width>
class matrix {
public:
    typedef matrix<vt, std::max<size_t>(1, height - 1), std::max<size_t>(1, width - 1)> matrix_minor;

    matrix_minor create_minor(const matrix<vt, height, width>& other, const size_t &i_del, const size_t &j_del) const {
        matrix_minor res;
        for (size_t i = 0, cur_i = 0; i < height; ++i) {
            if (i == i_del)
                continue;
            for (size_t j = 0, cur_j = 0; j < width; ++j) {
                if (j == j_del)
                    continue;
                res.at(cur_i, cur_j++) = data_[i][j];
            }
            ++cur_i;
        }
        return res;
    }

    explicit matrix();

    explicit matrix(const vt &num);

    matrix(const matrix<vt, height, width> &other);

    ~matrix();

    matrix &operator=(const matrix<vt, height, width> &other);

    matrix &operator+=(const matrix<vt, height, width> &other);

    matrix operator+(const matrix<vt, height, width> &other) const;

    matrix &operator+=(const vt &value);

    matrix operator+(const vt &value) const;

    matrix &operator-=(const matrix<vt, height, width> &other);

    matrix operator-(const matrix<vt, height, width> &other) const;

    matrix &operator-=(const vt &value);

    matrix operator-(const vt &value) const;

    matrix &operator*=(const long long int &value);

    [[ nodiscard ]] matrix operator*(const long long int &value) const;

    matrix &operator*=(const matrix<vt, width, width>& other);

    template <size_t W>
    matrix<vt, height, W> operator*(const matrix<vt, width, W>& other) const;

    matrix operator+() const { return *this; }

    matrix operator-() const;

    vt &at(const size_t &i, const size_t &j);

    [[nodiscard]] const vt &at(const size_t &i, const size_t &j) const;

    [[nodiscard]] matrix<vt, width, height> transposed() const;

    vt trace() const;

    vt det() const;

    template<typename T, size_t H, size_t W>
    bool operator==(const matrix<T, H, W> &other) const;

    template<typename T, size_t H, size_t W>
    bool operator!=(const matrix<T, H, W> &other) const { return !(*this == other); }

    template<typename T, size_t H, size_t W>
    friend std::ostream &operator<<(std::ostream &out, const matrix<T, H, W> &matrix);

private:
    vt **data_;

    void allocate_memory();

    void deallocate_memory();
};

#endif //MATRIX_MATRIX_H
