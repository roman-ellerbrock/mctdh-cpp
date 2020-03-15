//
// Created by Roman Ellerbrock on 3/13/20.
//

#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "Core/Vector.h"
#include "TreeOperators/MultiLeafOperator.h"

class Potential {
public:
	Potential() = default;
	virtual ~Potential() = default;

	virtual void Initialize(const Tree& tree) {}

	virtual vector<MLOcd> Build(const Tree& tree) {
		MLOcd v_mlo;
		PotentialOperator v(tree.nLeaves() - 1, 0);
		v_mlo.V() = v;
		return {v_mlo};
	}

	virtual double Evaluate(const Vectord& Xv, size_t part) const = 0;
};


#endif //POTENTIAL_H
