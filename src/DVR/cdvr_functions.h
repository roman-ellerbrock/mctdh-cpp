//
// Created by Roman Ellerbrock on 3/18/20.
//

#ifndef CDVR_FUNCTIONS_H
#define CDVR_FUNCTIONS_H
#include "Core/Wavefunction.h"
#include "TreeClasses/MatrixTree.h"
#include "DVR/DeltaVTree.h"
#include "ExplicitEdgeWavefunction.h"
#include "TreeGrids.h"

namespace cdvr_functions {

	void fillXNode(Vectord& X, vector<size_t> idx, const TreeGrids& grids, const TreeGrids& holegrids,
		const Node& node);

	void fillXEdge(Vectord& X, vector<size_t> idx, const TreeGrids& grids, const TreeGrids& holegrids,
		const Node& node);

	void CalculateDeltaVs(DeltaVTree& deltaVTree, const ExplicitEdgeWavefunction& Chi,
		const TensorTreecd& Vnodes, const MatrixTreed& Vedges, const Tree& tree);

	TensorTreecd Apply(const ExplicitEdgeWavefunction& Chi, const TensorTreecd& V,
		const DeltaVTree& DeltaVs, const Tree& tree);
}

#endif //CDVR_FUNCTIONS_H
