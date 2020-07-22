//
// Created by Roman Ellerbrock on 3/13/20.
//

#ifndef CDVR_H
#define CDVR_H
#include "DVR/TDDVR.h"
#include "TreeOperators/Potential.h"
#include "DVR/DeltaVTree.h"
#include "DVR/cdvr_functions.h"

class CDVR {
public:
	CDVR(const Tree& tree);
	CDVR(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree, size_t part = 0);
	~CDVR() = default;

	void Update(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree, size_t part = 0);

	void Update(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree,
		const Tree& smalltree, size_t part = 0);

	Tensorcd Apply(Tensorcd Phi, const Matrixcd& invsq_rho, const Node& node) const;

	TensorTreecd Apply(const Wavefunction& Psi, const Tree& tree) const;

	TDDVR tddvr_;

private:
	TensorTreecd Vnode_;
	MatrixTreed Vedge_;

	DeltaVTree deltaV_;

	ExplicitEdgeWavefunction Chi_;
	TensorTreecd Cdown_;
};

#endif //CDVR_H
