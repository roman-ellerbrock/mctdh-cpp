//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef TINTEGRATORVARIABLES_H
#define TINTEGRATORVARIABLES_H
#include "Core/Wavefunction.h"

template <typename T>
class tIntegratorVariables {
public:
	tIntegratorVariables(double time_now_, double time_end_, double dt_, double out_,
		double accuracy_root_, double accuracy_leaf_, TensorTree<T>& Psi_,
		const Hamiltonian& h_, const SOP<T>& hsop_, const Tree& tree_, const Tree& cdvrtree_,
		string ofname_, string ifname_, bool save_psi_,
		bool cmf_bottom = true, bool cmf_upper = true, bool cmf_top = true)
		: time_now(time_now_), time_end(time_end_), dt(dt_), out(out_),
		accuracy_root(accuracy_root_), accuracy_leaf(accuracy_leaf_), psi(&Psi_),
		h(&h_), sop(&hsop_), tree(&tree_), cdvrtree(&cdvrtree_), ofname(move(ofname_)),
		ifname(move(ifname_)), save_psi(save_psi_),
		cmf_bottom_(cmf_bottom), cmf_upper_(cmf_upper), cmf_top_(cmf_top)
		{}

	double time_now;
	double time_end;
	double dt;
	double out;

	double accuracy_root;
	double accuracy_leaf;

	TensorTree<T>* psi;
	const Hamiltonian* h;
	const SOP<T>* sop;
	const Tree* tree;
	const Tree* cdvrtree;

	string ofname;
	string ifname;

	bool save_psi;

	bool cmf_bottom_;
	bool cmf_upper_;
	bool cmf_top_;
};


#endif //TINTEGRATORVARIABLES_H
