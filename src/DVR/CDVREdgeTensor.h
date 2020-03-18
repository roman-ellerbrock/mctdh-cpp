//
// Created by Roman Ellerbrock on 3/18/20.
//

#ifndef CDVREDGETENSOR_H
#define CDVREDGETENSOR_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeShape/Tree.h"
#include "Core/Wavefunction.h"
#include "DVR/TreeGrids.h"
#include "TreeClasses/MatrixTree.h"

class CDVREdgeTensor : NodeAttribute<Tensorcd>{
public:
	CDVREdgeTensor() = default;
	~CDVREdgeTensor() = default;

	void Initialize(const Tree& tree);

	void Calculate(const Wavefunction& Psi, const TensorTreecd& dvr, const MatrixTreed& edge_dvr,
		const Tree& tree);

};


#endif //CDVREDGETENSOR_H
