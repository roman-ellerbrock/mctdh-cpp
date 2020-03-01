//
// Created by Roman Ellerbrock on 2/27/20.
//

#ifndef MCTDH_STATE_H
#define MCTDH_STATE_H
#include "TreeShape/Tree.h"
#include "Core/Hamiltonian.h"
#include "Core/Wavefunction.h"

struct mctdh_state {
	mctdh_state():rng_(time(nullptr)) {}
	~mctdh_state() = default;

	Tree tree_;
	shared_ptr<Hamiltonian> hamiltonian_;
	map<string, Wavefunction> wavefunctions_;

	mt19937 rng_;

};


#endif //MCTDH_STATE_H
