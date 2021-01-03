//
// Created by Roman Ellerbrock on 3/9/20.
//

#ifndef TDDVR_H
#define TDDVR_H
#include "DVR/TreeGrids.h"
#include "Core/Wavefunction.h"
#include "DVR/XMatrixTrees.h"
#include "TreeClasses/MatrixTree.h"
#include "DVR/MatrixTensorTree.h"

class TDDVR {
public:
	TDDVR(const Wavefunction& Psi, const Tree& tree)
		: TDDVR(tree) {
		Update(Psi, tree);
	}

	explicit TDDVR(const Tree& tree)
		: Xs_(tree), rho_(tree), trafo_(tree), hole_trafo_(tree),
		  grids_(tree), hole_grids_(tree, true) {
	}

	~TDDVR() = default;

	void GridTransformation(MatrixTensorTree& Psi, const Tree& tree, bool inverse = false) const;

	void Update(const Wavefunction& Psi, const Tree& tree);

	void print(const Tree& tree) const;

	void NodeTransformation(Tensorcd& Phi, const Node& node, bool inverse) const;

	TreeGrids grids_;
	MatrixTreecd trafo_;

	TreeGrids hole_grids_;
	MatrixTreecd hole_trafo_;

	XMatrixTrees Xs_;
private:
	MatrixTreecd rho_;

	void NodeTransformation(Wavefunction& Psi, const Tree& tree, bool inverse) const;

	void EdgeTransformation(Matrixcd& B_inv, const Edge& edge, bool inverse) const;
	void EdgeTransformation(MatrixTreecd& B_inv, const Tree& tree, bool inverse) const;

	void GridTransformation(Wavefunction& Psi, const Tree& tree, bool inverse = false) const;
	void GridTransformationLocal(Tensorcd& Phi, const Node& node, bool inverse = false) const;
};


#endif //TDDVR_H
