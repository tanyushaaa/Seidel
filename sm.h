#pragma once
#include "matrix.h"
#include "cmath"

float norm(float* vec1, float* vec2, int n);
std::pair<float**, int> new_mtrx(std::pair<float**, int> m, float* x);
float* iteration(std::pair<float**, int> m, float* x);
float* sm(const Matrix& mtrx);
