#pragma once

struct Vect {
  float x, y;

  Vect operator*(float rhs) const { return {x * rhs, y * rhs}; }
  Vect operator/(float rhs) const { return {x / rhs, y / rhs}; }
  Vect operator+(Vect rhs) const { return {x + rhs.x, y + rhs.y}; }
  float dot_product(Vect rhs) const { return rhs.x * x + rhs.y * y; }
};
