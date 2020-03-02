//
// Created by Roman Ellerbrock on 2/27/20.
//

#include "Parser/yaml_parser.h"
#include "yaml-cpp/yaml.h"
#include "TreeShape/TreeFactory.h"
#include "Hamiltonians.h"

namespace parser {
	Leaf create_leaf(const YAML::Node& yaml_node) {
		auto dim = yaml_node["dim"].as<size_t>();
		auto mode = yaml_node["mode"].as<size_t>();
		auto type = yaml_node["leaftype"].as<size_t>();
		PhysPar par;
		return Leaf(dim, mode, type, 0, par);
	}

	Node create_node(const YAML::Node& yaml_node) {
		Node node;
		auto dim = yaml_node["dim"].as<size_t>();
		vector<size_t> dims;
		for (const auto& child : yaml_node["children"]) {
			auto type = child["type"].as<string>();
			dims.push_back(child["dim"].as<size_t>());
			if (type == "leaf") {
				Leaf leaf = create_leaf(child);
				node = Node(leaf, dim);
			} else if (type == "node") {
				node.push_back(create_node(child));
			} else {
				cerr << "Invalid node type." << endl;
				cerr << "Choices: (leaf, node)" << endl;
				exit(1);
			}
		}

		dims.push_back(dim);
		node.shape() = TensorShape(dims);
		return node;
	}

	void read_leaf_parameters(Tree& tree, const YAML::Node& node) {
		for (const auto& child : node["leaves"]) {
			auto mode = child["mode"].as<size_t>();
			auto& leaf = tree.GetLeaf(mode);
			PhysPar par;
			auto r0 = child["r0"].as<double>();
			par.setR0(r0);
			auto wfr0 = child["wfr0"].as<double>();
			par.setWFR0(wfr0);
			auto omega = child["omega"].as<double>();
			par.setOmega(wfr0);
			auto wfomega = child["wfomega"].as<double>();
			par.setWFOmega(wfr0);
			leaf.SetPar(par);
		}
	}

	Tree create_tree(const YAML::Node& node) {
		Tree tree;
		Node root = create_node(node["tree"]);
		tree.SetRoot(root);
		tree.Update();
		read_leaf_parameters(tree, node);
		return tree;
	}

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

	Tree read_tree(const string& filename) {
		YAML::Node yaml = YAML::LoadFile(filename);
		return read_tree(yaml);
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


