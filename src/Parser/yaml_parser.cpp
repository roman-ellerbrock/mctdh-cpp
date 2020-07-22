//
// Created by Roman Ellerbrock on 2/27/20.
//

#include "Parser/yaml_parser.h"
#include "yaml-cpp/yaml.h"
#include "TreeShape/TreeFactory.h"
#include "Hamiltonians.h"
#include "Core/Eigenstates.h"

namespace parser {

	template <typename T>
	T evaluate(const YAML::Node& node, const string& key) {
		T val;
		if (auto par = node[key]) {
			return par.as<T>();
		} else {
			cerr << "Did not specify key '" << key << "' in yaml node " << node << endl;
			exit(3);
		}
	}

	template <typename T>
	T evaluate(const YAML::Node& node, const string& key, T default_) {
		T val;
		if (auto par = node[key]) {
			return par.as<T>();
		} else {
			return default_;
		}
	}

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
			auto mode = evaluate<size_t>(child, "mode");
			auto& leaf = tree.GetLeaf(mode);
			auto r0 = evaluate<double>(child, "r0");
			auto wfr0 = evaluate<double>(child, "wfr0");
			auto omega = evaluate<double>(child, "omega");
			auto wfomega = evaluate<double>(child, "wfomega");
			auto& grid = leaf.PrimitiveGrid();
			grid.Initialize(omega, r0, wfr0, wfomega);
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
		auto type = evaluate<string>(node, "type");
		if (type == "balanced") {
			auto num_leaves = evaluate<size_t>(node, "number_leaves");
			auto dim_leaves = evaluate<size_t>(node, "dimension_leaves");
			auto dim_nodes = evaluate<size_t>(node, "dimension_nodes");
			return TreeFactory::BalancedTree(num_leaves, dim_leaves, dim_nodes);
		} else if (type == "manual") {
			return create_tree(node);
		} else if (type == "compact") {
			auto tree_str = evaluate<string>(node, "tree");
			stringstream ss(tree_str);
			Tree tree(ss);
//			tree.info();
/*			cout << "checking tree.." << endl;
			if (!tree.IsWorking()) {
				cerr << "Failed to read tree with .yaml parser in compact format.\n";
				cerr << "Error caused by tree input string reading:\n";
				cerr << tree_str << endl;
				exit(2);
			}
			cout << "tree checked." << endl;*/
			return tree;
		} else {
			cerr << "No valid tree type." << endl;
			cerr << "Choices: (balanced, manual)" << endl;
			exit(1);
		}
		return Tree();
	}

	shared_ptr<Hamiltonian> read_hamiltonian(const YAML::Node& node, const Tree& tree) {
		auto name = evaluate<string>(node, "name");
		shared_ptr<Hamiltonian> H_ptr(new Hamiltonian);
		Hamiltonian& H = *H_ptr;

		if (name == "coupled_ho") {
			H = CoupledHO(tree);
		} else if (name == "kinetic_energy") {
			H = Operator::KineticEnergy(tree);
		} else if (name == "nocl") {
			bool V = evaluate<bool>(node, "V", true);
			H = Operator::NOCl_H(V);
		} else if (name == "ch3_meanfield") {
			H = Operator::CH3_meanfield();
		} else if (name == "ch3_quasiexact") {
			CH3_quasiexact Hch3(tree);
			H = Hch3;
			cout << "YAML H size: " << H.size() << endl;
		} else {
			cout << "No valid Hamiltonian name." << endl;
			cout << "Choices: (coupled_ho)" << endl;
			exit(1);
		}
		assert(H_ptr->size() > 0);
		return H_ptr;
	}

	PotentialOperator set_potential(const YAML::Node& node, const Tree& tree) {
		auto name = evaluate<string>(node, "name");
		if (name == "coupled_ho") {
			auto V = make_shared<CDVRModelV>(tree.nLeaves());
			return PotentialOperator(V, 0, 0);
		} else if (name == "nocl") {
			auto state = evaluate<string>(node, "state", "S1");
			bool S1 = (state == "S1");
			auto V = make_shared<NOClPotential>(S1);
			return PotentialOperator(V, 0, 0);
		} else if(name == "ch3") {
			auto V = make_shared<CH3Potential>();
			Vectord mass(4);
			mass(0) = 12.0;
			mass(1) = 1.007;
			mass(2) = 1.007;
			mass(3) = 1.007;
			PotentialOperator Vop(V, 0, 0);
			Vop.Q_ = make_shared<TrafoCH3Quasiexact>(mass);
			return Vop;
		} else {
			cerr << "Did not recognise potential energy operator name\n";
			exit(1);
		}
		return PotentialOperator();
	}

	void new_wavefunction(mctdh_state& state, const YAML::Node& node) {
		auto name = evaluate<string>(node, "name");
		auto type = evaluate<string>(node, "type");
		if (type == "read") {
			auto filename = evaluate<string>(node, "filename");
			Wavefunction Psi(state.tree_);
			ifstream is(filename);
			is >> Psi;
			state.wavefunctions_[name] = Psi;
		} else if (type == "create") {
			state.wavefunctions_[name] = Wavefunction(state.rng_, state.tree_);
		} else {
			cerr << "No valid Wavefunction initialization type." << endl;
			cerr << "Choices: (read, create)" << endl;
			exit(1);
		}
	}

	IntegratorVariables new_ivar(const YAML::Node& node, mctdh_state& state) {
		auto t_end = evaluate<double>(node, "t_end", 100*41.362);
		auto t = evaluate<double>(node, "t", 0.);
		auto out = evaluate<double>(node, "out", 41.362);
		auto dt = evaluate<double>(node, "dt", 1.);
		auto cmf = evaluate<double>(node, "eps_cmf", 1e-4);
		auto bs = evaluate<double>(node, "eps_bs", 1e-5);
		auto file_in = evaluate<string>(node, "file_in", "in.dat");
		auto file_out = evaluate<string>(node, "file_out", "out.dat");
		auto save = evaluate<bool>(node, "save_psi", true);
		IntegratorVariables ivar(t, t_end, dt, out, cmf, bs,
			state.wavefunctions_["Psi"], *state.hamiltonian_,
			state.tree_, file_out, file_in, false);
		return ivar;
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
			} else if (name == "potential") {
				PotentialOperator V = set_potential(node, state.tree_);
				state.hamiltonian_->V_ = V;
				state.hamiltonian_->hasV = true;
			} else if (name == "wavefunction") {
				new_wavefunction(state, node);
			} else if (name == "eigenstates") {
				auto ivar = new_ivar(node, state);
				Eigenstates(ivar);
			} else if (name == "cmf") {
				auto ivar = new_ivar(node, state);
				const Hamiltonian& H = *ivar.h;
				const Tree& tree = *ivar.tree;
				CMFIntegrator cmf(H, tree, 1.);
				cmf.Integrate(ivar);
			}
		}
		return state;
	}
}


