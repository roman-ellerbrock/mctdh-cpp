//
// Created by Roman Ellerbrock on 2/27/20.
//

#ifndef TMCTDH_STATE_H
#define TMCTDH_STATE_H
#include "TreeShape/Tree.h"
#include "Core/Hamiltonian.h"
#include "Core/Wavefunction.h"

template <typename T>
struct tmctdh_state {
	tmctdh_state():rng_(0) {}
//	mctdh_state():rng_(time(nullptr)) {}
	~tmctdh_state() = default;

	Tree tree_;
	Tree cdvrtree_;
	shared_ptr<Hamiltonian> hamiltonian_;
	shared_ptr<SOP<T>> sop_;
	map<string, TensorTree<T>> wavefunctions_;

	mt19937 rng_;

	void print(ostream& os = cout) const {
		os << "Tree: " << endl;
		tree_.print(os);
		os << "Hamiltonian:" << endl;
		hamiltonian_->print(os);
		os << "Wavefunctions:" << endl;
		for (const auto& pair : wavefunctions_) {
			os << pair.first << endl;
			pair.second.print(tree_, os);
		}
	}

};


#endif //TMCTDH_STATE_H
