#pragma once

template <class T> struct Vect {
  T x, y;

  Vect operator*(T rhs) const { return {x * rhs, y * rhs}; }
  Vect operator/(T rhs) const { return {x / rhs, y / rhs}; }
  Vect operator+(Vect rhs) const { return {x + rhs.x, y + rhs.y}; }
  T dot_product(Vect rhs) const { return rhs.x * x + rhs.y * y; }
};
