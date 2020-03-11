//
// Created by Roman Ellerbrock on 3/9/20.
//

#ifndef TDDVR_H
#define TDDVR_H
#include "DVR/TreeGrids.h"
#include "Core/Wavefunction.h"
#include "DVR/XMatrixTrees.h"
#include "TreeClasses/MatrixTree.h"

class TDDVR {
public:
	explicit TDDVR(const Tree& tree) : Xs_(tree), rho_(tree), trafo_(tree), grids_(tree){
	}
	~TDDVR() = default;

	void Calculate(const Wavefunction& Psi, const Tree& tree);

	TreeGrids grids_;
//	TreeGrids hole_grids_;
	MatrixTreecd trafo_;
private:
	XMatrixTrees Xs_;
	MatrixTreecd rho_;

};


#endif //TDDVR_H
