#pragma once

#include <cmath>
template <class T> struct Vect {
  T x, y;

  Vect operator*(T rhs) const { return {x * rhs, y * rhs}; }
  Vect operator/(T rhs) const { return {x / rhs, y / rhs}; }
  Vect operator+(Vect rhs) const { return {x + rhs.x, y + rhs.y}; }
  Vect operator-(Vect rhs) const { return {x - rhs.x, y - rhs.y}; }
  T dot_product(Vect rhs) const { return rhs.x * x + rhs.y * y; }

  T euclid_dist(Vect rhs) const {
    Vect new_rhs = (*this - rhs);
    return sqrt(new_rhs.dot_product(new_rhs));
  }
};
