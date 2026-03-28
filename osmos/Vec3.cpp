#pragma once

#include <cmath>

struct Vec3 {
  double x, y, z;
  Vec3() : x(0.0), y(0.0), z(0.0) {}
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

  Vec3 operator+(const Vec3 &v) const { return {x + v.x, y + v.y, z + v.z}; }
  Vec3 operator-(const Vec3 &v) const { return {x - v.x, y - v.y, z - v.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

  Vec3 &operator+=(const Vec3 &v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  Vec3 &operator-=(const Vec3 &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  Vec3 &operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }

  Vec3 &operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  double dot(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }

  double norm_sq() const { return x * x + y * y + z * z; }

  double norm() const { return std::sqrt(norm_sq()); }
};