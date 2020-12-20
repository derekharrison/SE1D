/*
 * complex.cpp
 *
 *  Created on: Dec 19, 2020
 *      Author: d-w-h
 */

#include "complex.hpp"

Complex::Complex() {}

Complex::Complex(double a, double b) {
    this->a = a;
    this->b = b;
}

Complex Complex::operator+(const Complex& m) {
    Complex result_addition(0,0);
    result_addition.a = this->a + m.a;
    result_addition.b = this->b + m.b;

    return result_addition;
}

Complex Complex::operator*(const Complex& m) {
    Complex result_multiplication(0,0);
    result_multiplication.a = this->a*m.a - this->b*m.b;
    result_multiplication.b = this->a*m.b + this->b*m.a;

    return result_multiplication;
}

Complex Complex::operator-(const Complex& m) {
    Complex result_subtraction(0,0);
    result_subtraction.a = this->a - m.a;
    result_subtraction.b = this->b - m.b;

    return result_subtraction;
}
