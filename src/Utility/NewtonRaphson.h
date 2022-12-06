//
// Created by Roman Ellerbrock on 11/7/22.
//

#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H
#include <functional>
#include "Core/Matrix.h"

using GradientFunction = function<void(Vectord& grad, const Vectord& x)>;
using HessianFunction = function<void(Matrixd& H, const Vectord& grad, const Vectord& x)>;

void newtonRaphson(Vectord& x, Vectord& grad, Matrixd& hessian,
	GradientFunction& ddx, HessianFunction& dxdy,
	size_t n_iter, double eps) {

	Vectord xn = x;
	for (size_t i = 0; i < n_iter; ++i) {
		ddx(grad, x);
		dxdy(hessian, grad, x);
		auto hinv_x = inverse(diagonalize(hessian), eps);
		auto hinv = toMatrix(hinv_x);
		xn = hinv * grad;
		double delta = (x - xn).norm();
		if (delta < eps) { break; }
	}

}

#endif //NEWTONRAPHSON_H
