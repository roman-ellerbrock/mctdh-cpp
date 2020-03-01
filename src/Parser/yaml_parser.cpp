//
// Created by Roman Ellerbrock on 2/27/20.
//

#include "Parser/yaml_parser.h"
#include "yaml-cpp/yaml.h"
#include "TreeShape/TreeFactory.h"
#include "Hamiltonians.h"

namespace parser {

	Tree read_tree(const YAML::Node& node) {
		auto type = node["type"].as<string>();
		if (type == "balanced") {
			auto num_leaves = node["number_leaves"].as<size_t>();
			auto dim_leaves = node["dimension_leaves"].as<size_t>();
			auto dim_nodes = node["dimension_nodes"].as<size_t>();
			return TreeFactory::BalancedTree(num_leaves, dim_leaves, dim_nodes);
		} else {
			cerr << "No valid tree type." << endl;
			cerr << "Choices: (balenced, unbalenced, manual)" << endl;
			exit(1);
		}
		return Tree();
	}

	shared_ptr<Hamiltonian> read_hamiltonian(const YAML::Node& node, const Tree& tree) {
		auto name = node["name"].as<string>();
		shared_ptr<Hamiltonian> H_ptr(new Hamiltonian);
		Hamiltonian& H = *H_ptr;

		if (name == "coupled_ho") {
			H = CoupledHO(tree);
		} else {
			cerr << "No valid Hamiltonian name." << endl;
			cerr << "Choices: (coupled_ho)" << endl;
			exit(1);
		}
		assert(H_ptr->size() > 0);
		return H_ptr;

	}

	void new_wavefunction(mctdh_state& state, const YAML::Node& node) {
		auto type = node["type"].as<string>();
		auto name = node["name"].as<string>();
		if (type == "read") {
			auto filename = node["filename"].as<string>();
			state.wavefunctions_[name] = Wavefunction(filename);
		} else if (type == "create") {
			state.wavefunctions_[name] = Wavefunction(state.rng_, state.tree_);
		} else {
			cerr << "No valid Wavefunction initialization type." << endl;
			cerr << "Choices: (read, create)" << endl;
			exit(1);
		}
	}

	mctdh_state read(const string& yaml_filename) {
		YAML::Node config = YAML::LoadFile(yaml_filename);
		mctdh_state job;
		job.tree_ = read_tree(config["tree"]);
		job.hamiltonian_ = read_hamiltonian(config["hamiltonian"], job.tree_);
		new_wavefunction(job, config["wavefunction"]);
		return job;
	}

	mctdh_state run(const string& yaml_filename) {
		YAML::Node config = YAML::LoadFile(yaml_filename);
		mctdh_state state;
		for (const auto& node : config["run"]) {
			const auto& name = node["job"].as<string>();
			if (name == "tree") {
				state.tree_ = read_tree(node);
			} else if (name  == "hamiltonian") {
				state.hamiltonian_ = read_hamiltonian(node, state.tree_);
			} else if (name == "wavefunction") {
				new_wavefunction(state, node);
			}
		}
		return state;
	}
}


