//
// Created by Roman Ellerbrock on 3/18/20.
//

#ifndef CDVR_FUNCTIONS_H
#define CDVR_FUNCTIONS_H
#include "Core/Wavefunction.h"
#include "TreeClasses/MatrixTree.h"
#include "DVR/DeltaVTree.h"

namespace cdvr_functions {

	void Update(DeltaVTree& deltaVTree, const Wavefunction& Psi, const TensorTreecd& dvr,
		const MatrixTreed& edge_dvr, const Tree& tree);

}

#endif //CDVR_FUNCTIONS_H
