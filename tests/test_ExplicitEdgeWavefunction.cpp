//
// Created by Roman Ellerbrock on 4/3/20.
//

#include "UnitTest++/UnitTest++.h"
#include "DVR/ExplicitEdgeWavefunction.h"
#include "Parser/yaml_parser.h"

SUITE (ExplicitEdgeWavefunction) {
	double eps = 1e-7;

	TEST (Init) {
		string yaml_filename("../examples/ho_sl.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		/// Transform to canonical representation

		ExplicitEdgeWavefunction Chi(Psi, tree, true);

		/// Check bottom-up
		{
			const TensorTreecd& Atilde = Chi.nodes();
			auto rho = TreeFunctions::Contraction(Psi, tree, true);
			for (const Node& node : tree) {
				if (!node.isToplayer()) {
					auto x = Contraction(Atilde[node], Atilde[node], node.nChildren());
					auto r = Residual(x, rho[node]);
						CHECK_CLOSE(0., r, eps);
				}
			}

			/// Check top-down
			for (const Node& node : tree) {
				if (!node.isBottomlayer()) {
					for (size_t k = 0; k < node.nChildren(); ++k) {
						const Node& child = node.child(k);
						auto x = Contraction(Atilde[node], Atilde[node], k);
						auto r = Residual(x, rho[child]);
							CHECK_CLOSE(0., r, eps);
					}
				}
			}
		}

			CHECK_EQUAL(true, IsWorking(Chi, tree, eps));
	}
}

