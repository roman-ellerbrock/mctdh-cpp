//
// Created by Roman Ellerbrock on 6/7/22.
//

#include "UnitTest++/UnitTest++.h"
#include "Parser/yaml_parser.h"
#include "DVR/TDDVR.h"

SUITE (TDDVR) {
	double eps = 1e-7;

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

	TEST (XMat) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		SOPcd X_ops;
		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.getLeaf(l);
			X_ops.push_back(MLOcd(&LeafInterface::applyX, leaf.mode()), 1.);
		}

		MatrixTreescd Xs(X_ops, tree);
		TreeFunctions::represent(Xs, X_ops, Psi, Psi, tree);
			CHECK_EQUAL(4, Xs.matrices_.size());
			CHECK_EQUAL(4, Xs.contractions_.size());

/*		cout << "Xs mat:\n";
		size_t i = 0;
		for (const auto& x : Xs.matrices_) {
			cout << "i: " << i++ << endl;
			x.print();
		}*/
	}



}

