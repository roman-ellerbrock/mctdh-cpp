//
// Created by Roman Ellerbrock on 3/13/20.
//

#include "DVR/CDVR.h"
#include "Parser/yaml_parser.h"
#include "UnitTest++/UnitTest++.h"
#include "../Hamiltonians/PESs/PESs.h"

SUITE(CDVR) {
	TEST(coupledHO) {
		string yaml_filename("../examples/ho_sl.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		CDVRModelV V(4);

		CDVR cdvr(Psi, V, tree);


	}
}