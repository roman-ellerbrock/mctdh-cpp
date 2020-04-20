//
// Created by Roman Ellerbrock on 3/10/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Parser/yaml_parser.h"
#include "TreeClasses/SOPMatrixTrees.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "DVR/XMatrixTrees.h"
#include "DVR/TDDVR.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "Util/QMConstants.h"

SUITE (TDDVR) {
	class Parser {
	public:
		Parser() {
			string yaml_filename("../examples/coupledho.yaml");
			state_ = parser::run(yaml_filename);
			Psi_ = state_.wavefunctions_["Psi"];
		}
		mctdh_state state_;
		Wavefunction Psi_;
	};

	TEST(XMat) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		SOPcd X_ops;
		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			X_ops.push_back(MLOcd(&LeafInterface::applyX, leaf.Mode()), 1.);
		}
		MatrixTreescd Xs(X_ops, tree);
		TreeFunctions::Represent(Xs, X_ops, Psi, Psi, tree);
			CHECK_EQUAL(4, Xs.matrices_.size());
			CHECK_EQUAL(4, Xs.contractions_.size());
	}

	TEST(Xs) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		XMatrixTrees Xs(tree);
		Xs.Update(Psi, tree);
			CHECK_EQUAL(4, Xs.holes_.size());
			CHECK_EQUAL(4, Xs.mats_.size());
	}

	TEST(TDDVR) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		TDDVR tddvr(Psi, tree);
		ExplicitEdgeWavefunction Chi(Psi, tree, true);
		tddvr.GridTransformation(Chi, tree, false);
		auto Xi = Chi.BottomUpNormalized(tree);
		XMatrixTrees Xs(tree);
		Xs.Update(Xi, tree);
		Xs.print();
//		tddvr.print(tree);
	}

}


