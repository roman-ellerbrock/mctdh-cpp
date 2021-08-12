//
// Created by Roman Ellerbrock on 7/27/21.
//

#ifndef JORDANWIGNER_H
#define JORDANWIGNER_H
#include "TreeOperators/LeafMatrix.h"
#include "TreeOperators/SumOfProductsOperator.h"


namespace JordanWigner {
	typedef vector<tuple<int, int, double>> TwoIndex;

	typedef vector<tuple<int, int, int, int, double>> FourIndex;

	MLOd twoIndexOperator(size_t p, size_t q, double eps);

	MLOd fourIndexOperator(size_t p, size_t q, size_t r, size_t s, double eps);

	SOPd electronicHamiltonian(const TwoIndex& hpq, const FourIndex& Hpqrs);
}

#endif //JORDANWIGNER_H
