//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef INTEGRATORVARIABLES_H
#define INTEGRATORVARIABLES_H
#include "Core/Wavefunction.h"

struct IntegratorVariables {
	IntegratorVariables(double time_now_, double time_end_, double dt_, double out_,
		double accuracy_root_, double accuracy_leaf_, Wavefunction& Psi_,
		const Hamiltonian& h_, const Tree& tree_, string ofname_, string ifname_,
		bool save_psi_):time_now(time_now_), time_end(time_end_), dt(dt_), out(out_),
		accuracy_root(accuracy_root_), accuracy_leaf(accuracy_leaf_), psi(&Psi_),
		h(&h_), tree(&tree_), ofname(move(ofname_)), ifname(move(ifname_)), save_psi(save_psi_) {}

	double time_now;
	double time_end;
	double dt;
	double out;

	double accuracy_root;
	double accuracy_leaf;

	Wavefunction* psi;
	const Hamiltonian* h;
	const Tree* tree;

	string ofname;
	string ifname;

	bool save_psi;

};


#endif //INTEGRATORVARIABLES_H
