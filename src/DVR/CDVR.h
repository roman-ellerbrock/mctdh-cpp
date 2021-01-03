//
// Created by Roman Ellerbrock on 3/13/20.
//

#ifndef CDVR_H
#define CDVR_H
#include "DVR/TDDVR.h"
#include "TreeOperators/Potential.h"
#include "DVR/DeltaVTree.h"
#include "DVR/cdvr_functions.h"
#include "TreeClasses/TensorTreeFunctions.h"

class CDVR {
public:
	CDVR(const Tree& tree);
//	CDVR(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree, size_t part = 0);
	~CDVR() = default;

	void Update(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree,
		size_t part = 0, bool out = false, ostream& os = cout);

	void Update2(Wavefunction Psi, const PotentialOperator& V,
		const Tree& smalltree, size_t part = 0,
		bool out = false, ostream& os = cout);

//	Tensorcd Apply(Tensorcd Phi, const Matrixcd& invsq_rho, const Node& node) const;
	Tensorcd Apply(Tensorcd Phi, const SpectralDecompositioncd& rho_x, const Node& node) const;

	TDDVR tddvr_;

private:
	TensorTreecd Vnode_;
	MatrixTreed Vedge_;

	DeltaVTree deltaV_;

	MatrixTensorTree Chi_;
	TensorTreecd Cdown_;

	Tree ltree_;
};

#endif //CDVR_H
