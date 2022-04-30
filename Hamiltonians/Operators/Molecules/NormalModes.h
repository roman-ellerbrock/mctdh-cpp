//
// Created by Roman Ellerbrock on 4/29/22.
//

#ifndef NORMALMODES_H
#define NORMALMODES_H
#include "TreeOperators/CoordinateTransformation.h"


class NormalModes : public CoordinateTransformation {
public:
	NormalModes(const string& file, size_t nmodes) {
		ifstream is(file);
		size_t dim1 = 1;
		size_t dim2 = 1;
		is >> dim1 >> dim2;
		Matrixd U(dim1, dim2);
		for (size_t i = 0; i < dim1; ++i) {
			for (size_t j = 0; j < dim2; ++j) {
				is >> U(j, i);
			}
		}
		U = Matrixd(dim1, nmodes);
		int shift = dim2 - nmodes;
		for (size_t i = 0; i < dim1; ++i) {
			for (size_t j = 0; j < nmodes; ++j) {
				U_(i, j) = U(i, shift + j);
			}
		}
	}

	Vectord transform(const Vectord& q) const override;
protected:
	Matrixd U_;
};


#endif //NORMALMODES_H
