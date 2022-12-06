//
// Created by Roman Ellerbrock on 10/31/22.
//

#ifndef QUADRATICSOP_H
#define QUADRATICSOP_H
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeClasses/MatrixTreeFunctions.h"

class QuadraticSOP {
public:
	QuadraticSOP() = default;
	~QuadraticSOP() = default;

	QuadraticSOP(const vector<MLOcd>& Ms, Matrixcd c,
		const Tree& tree) : c_(move(c)) {
		for (const auto& M : Ms) {
			hs_.push_back(M);
			mats_.emplace_back(MatrixTreecd(tree));
			hole_.emplace_back(MatrixTreecd(tree));
		}
	}

	void calculate(const TensorTreecd& Bra,
		const TensorTreecd& Ket, const Tree& tree) {
		/// build with O(N^2 chi^4)
		for (size_t i = 0; i < hs_.size(); ++i) {
			auto hKet = hs_[i].apply(Ket, tree);
			TreeFunctions::dotProduct(mats_[i], Bra, hKet, tree);
			TreeFunctions::contraction(hole_[i], Bra, Ket, mats_[i], tree);
		}
	}

	Matrixcd c_;

private:
	vector<MLOcd> hs_;
	vector<MatrixTreecd> mats_;
	vector<MatrixTreecd> hole_;
};

#endif //QUADRATICSOP_H
