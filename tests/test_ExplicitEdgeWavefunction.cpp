//
// Created by Roman Ellerbrock on 4/3/20.
//

#include "UnitTest++/UnitTest++.h"
#include "DVR/ExplicitEdgeWavefunction.h"
#include "Parser/yaml_parser.h"

SUITE(ExplicitEdgeWavefunction) {
	TEST(Init) {

		string yaml_filename("../examples/ho_sl.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		ExplicitEdgeWavefunction Chi(Psi, tree);
	}
}