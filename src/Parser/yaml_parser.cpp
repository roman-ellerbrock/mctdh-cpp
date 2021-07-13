//
// Created by Roman Ellerbrock on 2/27/20.
//

#include <string>
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
			auto& leaf = tree.getLeaf(mode);
			auto r0 = evaluate<double>(child, "r0");
			auto wfr0 = evaluate<double>(child, "wfr0");
			auto omega = evaluate<double>(child, "omega");
			auto wfomega = evaluate<double>(child, "wfomega");
			auto& grid = leaf.interface();
			grid.initialize(omega, r0, wfr0, wfomega);
		}
	}

	Tree create_tree(const YAML::Node& node) {
		Tree tree;
		Node root = create_node(node["tree"]);
		tree.setRoot(root);
		tree.update();
		read_leaf_parameters(tree, node);
		return tree;
	}

	Tree read_tree(const YAML::Node& node) {
		auto type = evaluate<string>(node, "type");
		if (type == "balanced") {
			auto num_leaves = evaluate<size_t>(node, "number_leaves");
			auto dim_leaves = evaluate<size_t>(node, "dimension_leaves");
			auto dim_nodes = evaluate<size_t>(node, "dimension_nodes");
			return TreeFactory::balancedTree(num_leaves, dim_leaves, dim_nodes);
		} else if (type == "manual") {
			return create_tree(node);
		} else if (type == "compact") {
			auto tree_str = evaluate<string>(node, "tree");
			stringstream ss(tree_str);
			Tree tree(ss);
			tree.info();
			cout << "checking tree.." << endl;
			if (!tree.isWorking()) {
				cerr << "Failed to read tree with .yaml parser in compact format.\n";
				cerr << "Error caused by tree input string reading:\n";
				cerr << tree_str << endl;
				exit(2);
			}
			cout << "tree checked." << endl;
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
		} else if (name == "exciton") {
			H = Operator::Exciton("matrix.out", tree);
		} else if (name == "ch3_quasiexact") {
            CH3_quasiexact Hch3(tree);
            H = Hch3;
            cout << "YAML H size: " << H.size() << endl;
        }else if (name == "schaepers") {
		    // find the masses supplied for this hamiltonian
            auto masses = evaluate<string>(node, "masses");
            stringstream masses_ss(masses);
            vector<double> massvec;
            while(masses_ss.good()){
                string substr;
                getline(masses_ss, substr, ',');
                massvec.push_back(stod(substr));
            }

            // find the coupling construction for this hamiltonian
            auto coupling = evaluate<string>(node, "coupling");
            stringstream coupling_ss(coupling);
            vector<int> couplingvec;
            while(coupling_ss.good()){
                string substr;
                getline(coupling_ss, substr, ',');
                couplingvec.push_back(stoi(substr));
            }

            // init schaepers vector
            H = Operator::schaepers(tree, massvec, couplingvec);

		} else {
			cout << "No valid Hamiltonian name." << endl;
			cout << "Chosen name: " << name << endl;
			cout << "See Parser/yaml_parser.cpp function read_hamiltonian for a list of choices.\n" << endl;
			exit(1);
		}
		assert(H_ptr->size() > 0);
		return H_ptr;
	}

	PotentialOperator set_potential(const YAML::Node& node, const Tree& tree) {
		auto name = evaluate<string>(node, "name");
		if (name == "coupled_ho") {
			auto coupling = evaluate<bool>(node, "coupling", "true");
			auto V = make_shared<CDVRModelV>(tree.nLeaves(), coupling);
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
		} else if(name == "liuch4cl") {
            // find the masses supplied for this potential
            auto masses = evaluate<string>(node, "masses");
            stringstream masses_ss(masses);
            vector<double> massvec;
            while(masses_ss.good()){
                string substr;
                getline(masses_ss, substr, ',');
                massvec.push_back(stod(substr));
            }
            // find the coupling construction for this potential
            auto coupling = evaluate<string>(node, "coupling");
            stringstream coupling_ss(coupling);
            vector<int> couplingvec;
            while(coupling_ss.good()){
                string substr;
                getline(coupling_ss, substr, ',');
                couplingvec.push_back(stoi(substr));
            }

            auto V = make_shared<liuch4cl>(massvec,couplingvec);
            PotentialOperator Vop(V, 0, 0);
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
            if(evaluate<string>(node, "filename").empty()){
                cerr << "supply filename to save and read directive" << endl;
                exit(1);
            }
            auto filename = evaluate<string>(node, "filename");
			Wavefunction Psi(state.tree_);
			ifstream is(filename);
			Psi.read(is);
			state.wavefunctions_[name] = Psi;
			is.close();
		} else if (type == "create") {
			state.wavefunctions_[name] = Wavefunction(state.rng_, state.tree_);
		} else if (type == "save") {
		    if(evaluate<string>(node, "filename").empty()){
		        cerr << "supply filename to save and read directive" << endl;
		        exit(1);
		    }
		    auto filename = evaluate<string>(node, "filename");
		    Wavefunction Psi(state.tree_);
		    Psi = state.wavefunctions_[name];
		    ofstream of(filename);
		    Psi.write(of);
		    of.close();
		} else {
			cerr << "No valid Wavefunction initialization type." << endl;
			cerr << "Choices: (read, create, save)" << endl;
			exit(1);
		}
	}

	IntegratorVariables new_ivar(const YAML::Node& node, mctdh_state& state) {
	    // TODO: this routine does not use the 'save'-directive
	    // TODO: also, it does not save the wavefunction and only works with a wavefunction called "Psi"
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
			state.tree_, state.cdvrtree_, file_out, file_in, false);
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
			const auto& job = node["job"].as<string>();
			if (job == "tree") {
				state.tree_ = read_tree(node);
				state.cdvrtree_ = state.tree_;
//				if (state.cdvrtree_.nNodes() == 0) { state.cdvrtree_ = state.tree_; }
			} else if (job == "hamiltonian") {
				state.hamiltonian_ = read_hamiltonian(node, state.tree_);
			} else if (job == "potential") {
				PotentialOperator V = set_potential(node, state.tree_);
				state.hamiltonian_->V_ = V;
				state.hamiltonian_->hasV = true;
			} else if (job == "wavefunction") {
				new_wavefunction(state, node);
			} else if (job == "eigenstates") {
				auto ivar = new_ivar(node, state);
				Eigenstates(ivar);
			} else if (job == "cmf") {
				auto ivar = new_ivar(node, state);
				const Hamiltonian& H = *ivar.h;
				const Tree& tree = *ivar.tree;
				const Tree& cdvrtree = *ivar.cdvrtree;
				CMFIntegrator cmf(H, tree, cdvrtree, 1.);
				cmf.Integrate(ivar);
			} else if (job == "cdvrtree") {
				state.cdvrtree_ = read_tree(node);
			}
		}
		return state;
	}
}


