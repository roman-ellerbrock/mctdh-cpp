//
// Created by Roman Ellerbrock on 3/13/20.
//

#ifndef CDVR_H
#define CDVR_H
#include "DVR/TDDVR.h"
#include "DVR/Potential.h"
#include "DVR/DeltaVTree.h"
#include "DVR/cdvr_functions.h"

class CDVR {
public:
	CDVR(const Wavefunction& Psi, const Potential& V, const Tree& tree, size_t part = 0);
	~CDVR() = default;

	void Update(const Wavefunction& Psi, const Potential& V, const Tree& tree, size_t part = 0);

	TDDVR tddvr_;

private:
	TensorTreecd Vnode_;
	MatrixTreed Vedge_;

	DeltaVTree deltaV_;
};

#endif //CDVR_H
