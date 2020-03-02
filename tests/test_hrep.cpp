//
// Created by Roman Ellerbrock on 3/1/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Core/HamiltonianRepresentation.h"
#include "TreeShape/TreeFactory.h"
#include "Hamiltonians.h"
#include "Parser/yaml_parser.h"

SUITE(HamiltonianRepresentation) {
	TEST(Build) {
		mt19937 gen;
		Tree tree = TreeFactory::BalancedTree(4, 6, 4);
		Hamiltonian H;
		H = CoupledHO(tree);
		HamiltonianRepresentation hrep(H, tree);
		Wavefunction Psi(gen, tree);
		hrep.build(H, Psi, tree);
	}

	TEST(BuildTree) {
		string yaml_filename("../examples/tree.yaml");
		YAML::Node yaml = YAML::LoadFile(yaml_filename);
		Tree tree = parser::create_tree(yaml);
	}
}