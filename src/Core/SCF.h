//
// Created by Roman Ellerbrock on 11/14/22.
//

#ifndef SCF_H
#define SCF_H
#include "TreeShape/Tree.h"
#include "Core/Wavefunction.h"
#include "Core/Hamiltonian.h"


vector<const Node*> scf_sweep(const Tree& tree);

struct SCF_parameters {
	size_t nIter{20};
	size_t nKrylov{5};
	double beta{1.};

	Wavefunction* psi{nullptr};
	const Hamiltonian* h{nullptr};
	const Tree* tree{nullptr};
};

struct KrylovSpace {
	KrylovSpace(vector<Tensorcd> space, SpectralDecompositioncd spectrum) :
		space_(space), spectrum_(spectrum) {}
	vector<Tensorcd> space_;
	SpectralDecompositioncd spectrum_;
};


void scf(SCF_parameters& par);

#endif //SCF_H
